from sklearn.model_selection import train_test_split
import tempfile
import h2o
import matplotlib.pyplot as plt
import numpy as np
from h2o.automl import H2OAutoML
import pandas as pd 
import os 

def display_h2o_out(aml, test):
    lb_df = aml.leaderboard.as_data_frame()
    exa = aml.explain(test)
    return lb_df, exa

def plot_corr(aml, test, y):

    # Get the best model automatically
    best_model = aml.leader
    print(f"Using model: {best_model.model_id}")

    # Check CV metrics (these are available)
    print("\nCross-Validation Metrics:")
    print(f"CV RMSE: {best_model.rmse(xval=True):.4f}")
    print(f"CV R²: {best_model.r2(xval=True):.4f}")
    print(f"CV MAE: {best_model.mae(xval=True):.4f}")

    # Get test set predictions (this always works)
    preds = best_model.predict(test).as_data_frame()
    actual = test[y].as_data_frame()

    # Plot test predictions
    plt.figure(figsize=(10, 5))

    # Subplot 1: Predicted vs Actual (predicted on x-axis, actual on y-axis)
    plt.subplot(1, 2, 1)
    plt.scatter(preds['predict'], actual, alpha=0.6)
    plt.plot([preds['predict'].min(), preds['predict'].max()], 
            [preds['predict'].min(), preds['predict'].max()], 
            'r--', label='Perfect predictions', linewidth=2)
    plt.xlabel('Predicted IC50', fontsize=12)
    plt.ylabel('Actual IC50', fontsize=12)
    plt.title('Test Set Predictions', fontsize=14)
    plt.legend()
    plt.grid(True, alpha=0.3)

    # Calculate metrics
    test_corr = np.corrcoef(actual.values.flatten(), preds['predict'].values)[0, 1]
    test_rmse = np.sqrt(np.mean((actual.values.flatten() - preds['predict'].values)**2))

    # Add text box with metrics
    textstr = f'Test R² = {test_corr**2:.3f}\nTest RMSE = {test_rmse:.3f}\n\nCV R² = {best_model.r2(xval=True):.3f}\nCV RMSE = {best_model.rmse(xval=True):.3f}'
    plt.text(0.05, 0.95, textstr, transform=plt.gca().transAxes, 
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Subplot 2: Residuals
    plt.subplot(1, 2, 2)
    residuals = actual.values.flatten() - preds['predict'].values
    plt.scatter(preds['predict'], residuals, alpha=0.6)
    plt.axhline(y=0, color='r', linestyle='--', linewidth=2)
    plt.xlabel('Predicted IC50', fontsize=12)
    plt.ylabel('Residuals', fontsize=12)
    plt.title('Residual Plot', fontsize=14)
    plt.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()

    print(f"\nTest Set Performance:")
    print(f"  Correlation: {test_corr:.3f}")
    print(f"  R²: {test_corr**2:.3f}")
    print(f"  RMSE: {test_rmse:.3f}")

    print(f"\nCross-Validation Performance:")
    print(f"  CV R²: {best_model.r2(xval=True):.3f}")
    print(f"  CV RMSE: {best_model.rmse(xval=True):.3f}")

def run_h2o(df, train_size=0.8, num_models=20):
    df = df.sample(frac=1, random_state=42).reset_index(drop=True)

    train_pd, test_pd = train_test_split(df, train_size=train_size, random_state=42)

    tmp1 = tempfile.NamedTemporaryFile(suffix='.csv', delete=False)
    train_pd.to_csv(tmp1.name, index=False)

    tmp2 = tempfile.NamedTemporaryFile(suffix='.csv', delete=False)
    test_pd.to_csv(tmp2.name, index=False)

    h2o.init(max_mem_size="50G")
    train = h2o.import_file(tmp1.name)
    test = h2o.import_file(tmp2.name)

    y = "mean_tgcpl_log_ic50 (uM)"
    x = [col for col in train.columns if col != y]

    aml = H2OAutoML(max_models=num_models, seed=42, project_name="h2o_tgcpl", nfolds=5)
    aml.train(x=x, y=y, training_frame=train, leaderboard_frame=test)

    lb_df, exa = display_h2o_out(aml, test)
    plot_corr(aml, test, y)

    return lb_df, exa

def compile_df(records_path, analysis_output_path, features=[]):
    metadata = pd.read_csv(os.path.join(records_path, 'metadata.csv'))
    exp_reads = pd.read_csv(os.path.join(records_path, 'experiment_readouts.csv'))

    if not features:
        df = exp_reads.drop(columns=['mean_hscpl_log_ic50 (uM)'])
        return df
    
    # boltz_confidence requires boltz_preds — must happen BEFORE the merge blocks
    if 'boltz_confidence' in features and 'boltz_preds' not in features:
        features.append('boltz_preds')
    
    df = pd.merge(exp_reads, metadata, on='substance_id', how='left')
    
    if 'interaction_fps' in features:
        fps = pd.read_csv(os.path.join(analysis_output_path, 'tgcpl_covalent/fingerprints/interaction_fingerprints.csv'))
        df = pd.merge(df, fps, on='substance_id', how='left')
    
    if 'boltz_preds' in features:
        mean_pred = pd.read_csv(os.path.join(analysis_output_path, 'tgcpl_covalent/predictions_averaged.csv'))
        mean_pred = mean_pred[~mean_pred['protein'].str.contains('5MAJ', na=False)]
        df = pd.merge(df, mean_pred, on='substance_id', how='left')
        df = df.dropna(subset=['pred_log10ic50'])

    # --- Column cleanup ---
    drop_patterns = ['std', 'pic50']
    
    drop_exact = [
    'mean_hscpl_log_ic50 (uM)', 'vault_mol_id_y', 'inchi_key_y', 'inchi_y', 'smiles_y',
    'vault_mol_id_x', 'inchi_key_x', 'inchi_x', 'synonyms',
    'smiles_x', 'smiles', 'name', 'protein', 'n_replicates',
    'substance_id', 'inchi', 'binding_probability',
    ]

    # drop ligand features if not requested
    ligand_feature_cols = ['molecular_weight', 'log_p', 'log_d', 'log_s',
       'num_aromatic_rings', 'num_h_bond_donors', 'num_h_bond_acceptors',
       'num_rule_of_5_violations', 'p_k_a', 'p_k_a_basic', 'heavy_atom_count',
       'topological_polar_surface_area', 'num_rotatable_bonds',
       'cns_mpo_score', 'bbb2_score', 'fsp3', 'p_k_a_acidic']
    if 'ligand_features' not in features:
        drop_exact += ligand_feature_cols

    # drop confidence metrics if not requested
    boltz_confidence_cols = [
        'binding_probability', 'confidence_score', 'complex_plddt',
        'ptm', 'iptm', 'ligand_iptm', 'protein_iptm',
        'complex_iplddt', 'complex_pde', 'complex_ipde', 
    ]
    if 'boltz_confidence' not in features:
        drop_exact += boltz_confidence_cols
    
    cols_to_drop = [c for c in df.columns if any(p in c for p in drop_patterns)]
    cols_to_drop += [c for c in drop_exact if c in df.columns]
    df = df.drop(columns=cols_to_drop)
    
    df = df.dropna(subset=['mean_tgcpl_log_ic50 (uM)'])
    
    return df