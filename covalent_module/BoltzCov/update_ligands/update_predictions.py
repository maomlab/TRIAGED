# use boltz_analysis_env 

import pandas as pd
import os 
import glob     
import shutil
from BoltzCov.analysis import analysis_utils

def fetch_new(pred_df, metadata_df, protein, min_replicates=3):
    '''
    Checks if predictions of compounds are in predictions.csv records.
    Works with long-format predictions.csv (with 'protein' column).
    Does NOT check if .cif file of final prediction is present.
    
    pred_df: df with names of molecules, confidence values and predictions (long format)
    metadata_df: df with all compound information
    protein: any name for the protein. taken from pdb file name.
    min_replicates: minimum number of non-NaN predictions required (default: 3)
    Returns DataFrame of new compounds that have no predictions yet OR have fewer than min_replicates. 
    '''
    
    # Filter predictions for this specific protein
    prot_preds = pred_df[pred_df['protein'] == protein].copy()
    pred_col = 'pred_log10ic50'  # Column name without protein suffix
    
    # If pred_col doesn't exist, return all compounds as new
    if pred_col not in prot_preds.columns:
        return metadata_df[['substance_id', 'inchi_key', 'smiles']].copy()
    
    # Count non-NaN predictions per compound for this protein
    pred_counts = prot_preds.groupby(['substance_id', 'inchi_key'])[pred_col].apply(
        lambda x: x.notna().sum()
    ).reset_index(name='pred_count')
    
    # Merge metadata with prediction counts
    merged = metadata_df.merge(
        pred_counts,
        on=['substance_id', 'inchi_key'], 
        how='left'
    )
    
    # Fill NaN counts with 0 (compounds with no predictions)
    merged['pred_count'] = merged['pred_count'].fillna(0)
    
    # Keep compounds with no predictions OR fewer than min_replicates
    no_pred_prot = merged[merged['pred_count'] < min_replicates][['substance_id', 'inchi_key', 'smiles']]
    
    return no_pred_prot

def check_attempted(error_df, dock_compounds):
    '''
    Checks errored.csv records to skip attempted ligands that failed in past runs.
    Returns updated DataFrame with errored compounds removed.
    '''
    # Anti-join: keep only rows NOT in error_df
    dock_compounds = dock_compounds.merge(
        error_df[['substance_id', 'inchi_key']], 
        on=['substance_id', 'inchi_key'], 
        how='left', 
        indicator=True
    )
    dock_compounds = dock_compounds[dock_compounds['_merge'] == 'left_only']
    dock_compounds = dock_compounds.drop('_merge', axis=1)
    
    return dock_compounds

def get_identifier(pkl_id, boltz_cache_prot, record_path):
    '''
    Gets unique identifier information about molecule.
    
    :param pkl_id: ID generated during yaml and pkl creation before docking. 
    :param boltz_cache_prot: Path to raw predictions for a given protein.
    :param record_path: Path to metadata.csv record.
    '''

    metadata_record = os.path.join(record_path, 'metadata.csv')
    metadata_df = pd.read_csv(metadata_record)
    csv_files = glob.glob(os.path.join(boltz_cache_prot, '*.csv'))
    if csv_files: ligand_cache_csv = csv_files[0] # temp csv that was generating before boltz run to make yamls 
    ligand_csv_df = pd.read_csv(ligand_cache_csv)
    # Get compound info
    substance_id = ligand_csv_df[ligand_csv_df['pkl_id'] == pkl_id]['substance_id'].iloc[0]
    meta_row = metadata_df[metadata_df['substance_id'] == substance_id].iloc[0]
    inchi_key = meta_row['inchi_key']
    vault_mol_id = meta_row['vault_mol_id']

    return substance_id, inchi_key, vault_mol_id

def update_preds(pkl_id, boltz_cache_prot, record_path, protein_name):
    ''' 
    Update predictions in predictions.csv. Each run creates a unique row (no averaging).
    Now saves data in long format with a 'protein' column.
    '''
    boltz_preds = analysis_utils.read_boltz_predictions(predictions_dir=boltz_cache_prot, reps=False)
    
    new_pred = boltz_preds[boltz_preds['compound_id'] == pkl_id]
    if new_pred.empty:
        return
    
    new_pred = new_pred.iloc[0]
    
    substance_id, inchi_key, vault_mol_id = get_identifier(pkl_id, boltz_cache_prot, record_path)
    
    pred_path = glob.glob(f'{boltz_cache_prot}/*/boltz_results_*/predictions/{pkl_id}*')
    if pred_path:
        boltz_runID = pred_path[0].split('/')[-4]
    else:
        boltz_runID = None
    
    predictions_csv = os.path.join(record_path, 'predictions.csv')
    
    numeric_cols = ['Pred log10(IC50)', 'Pred pIC50', 'Binding Probability', 'Confidence Score',
                    'PTM', 'IPTM', 'Ligand IPTM', 'Protein IPTM', 'Complex pLDDT', 
                    'Complex iPLDDT', 'Complex PDE', 'Complex iPDE']
    
    # Create new row WITHOUT protein suffix on metric columns
    # Add 'protein' as a separate column instead
    prediction_data = pd.DataFrame([{
        'substance_id': substance_id,
        'inchi_key': inchi_key,
        'vault_mol_id': vault_mol_id,
        'boltz_runID': boltz_runID,
        'protein': protein_name,  # Add protein as a separate column
        **{col.lower().replace(" ", "_").replace("(", "").replace(")", ""): new_pred[col] 
           for col in numeric_cols if col in new_pred.index}
    }])
    
    # Append new row (no averaging, boltz_runID makes each row unique)
    if os.path.exists(predictions_csv):
        old_preds = pd.read_csv(predictions_csv)
        updated_preds = pd.concat([old_preds, prediction_data], ignore_index=True)
        updated_preds.to_csv(predictions_csv, index=False)
    else:
        prediction_data.to_csv(predictions_csv, index=False)

def check_pred(run_cache_prot, record_path, VERBOSE, protein_name):
    '''
    Given the path to raw predictions, check if .cif was produced or not. 
    Updates errored.csv and predictions.csv.
    '''
    error_record = os.path.join(record_path, 'errored.csv')

    csv_files = glob.glob(os.path.join(run_cache_prot, '*.csv'))
    if csv_files: 
        ligand_cache_csv = csv_files[0]
    else:
        print(f"No ligand CSV found in {run_cache_prot}")
        return
        
    ligand_csv_df = pd.read_csv(ligand_cache_csv)

    for pkl_id in ligand_csv_df['pkl_id']:
        cif_files = glob.glob(f'{run_cache_prot}/*/boltz_results_*/predictions/{pkl_id}*/*.cif')
        json_files = glob.glob(f'{run_cache_prot}/*/boltz_results_*/predictions/{pkl_id}*/*.json')
        has_both = len(cif_files) > 0 and len(json_files) > 0

        # Get compound info
        substance_id, inchi_key, vault_mol_id = get_identifier(pkl_id, run_cache_prot, record_path)

        if has_both:
            update_preds(pkl_id, run_cache_prot, record_path, protein_name)  
            if VERBOSE: print(f"Updated predictions for {substance_id}")
        else: 
            error_data = pd.DataFrame([{
                'substance_id': substance_id,
                'inchi_key': inchi_key,
                'vault_mol_id': vault_mol_id
            }])
            
            if os.path.exists(error_record): 
                error_data.to_csv(error_record, mode='a', header=False, index=False)
            else: 
                error_data.to_csv(error_record, mode='w', header=True, index=False)
            if VERBOSE: print(f"Added {substance_id} to errored.csv")

    print("Record updates complete.")


def reorg_preds(run_cache_prot, record_path, output_dir, VERBOSE):
    '''
    Copy .cif and .yaml files to output directory for successful predictions.
    Files are renamed to boltz_runID_pklID format.
    Also copies hparams.yaml from one of the boltz_results directories to output_dir.
    '''
    predictions_csv = os.path.join(record_path, 'predictions.csv')
    if not os.path.exists(predictions_csv):
        print("[ERROR] No predictions.csv found")
        return
    
    pred_df = pd.read_csv(predictions_csv)
    
    # Create output subdirectories
    cif_dir = os.path.join(output_dir, 'cifs')
    yaml_dir = os.path.join(output_dir, 'yamls')
    os.makedirs(cif_dir, exist_ok=True)
    os.makedirs(yaml_dir, exist_ok=True)
    
    # Flag to track if hparams.yaml has been copied
    hparams_copied = False
    
    for _, row in pred_df.iterrows():
        substance_id = row['substance_id']
        boltz_runID = row['boltz_runID']
        
        # Find prediction files
        pred_files = glob.glob(f'{run_cache_prot}/*/boltz_results_*/predictions/*/*')
        
        # Find the .cif file
        cif_files = [f for f in pred_files if f.endswith('.cif')]
        
        if cif_files:
            cif_file = cif_files[0]
            
            # Copy .cif file with new name
            cif_dest = os.path.join(cif_dir, f"{boltz_runID}_{substance_id}.cif")
            shutil.copy(cif_file, cif_dest)
            if VERBOSE: print(f"Copied {substance_id}.cif for {boltz_runID}")
            
            # Find and copy yaml file
            # Extract parent directory: run_cache_prot/YYMMDD_HHMMSS_pklID/
            parent_dir = cif_file.split('/boltz_results_')[0]
            
            yaml_files = glob.glob(f'{parent_dir}/*.yaml')
            if yaml_files:
                yaml_file = yaml_files[0]
                yaml_dest = os.path.join(yaml_dir, f"{boltz_runID}_{substance_id}.yaml")
                shutil.copy(yaml_file, yaml_dest)
                if VERBOSE: print(f"Copied {substance_id}.yaml for {boltz_runID}")
            
            # Copy hparams.yaml only once (from first successful prediction)
            if not hparams_copied:
                # Extract boltz_results directory path
                boltz_results_dir = cif_file.split('/predictions/')[0]
                hparams_path = os.path.join(boltz_results_dir, 'lightning_logs', 'version_0', 'hparams.yaml')
                
                if os.path.exists(hparams_path):
                    hparams_dest = os.path.join(output_dir, 'hparams.yaml')
                    shutil.copy(hparams_path, hparams_dest)
                    if VERBOSE: print(f"Copied hparams.yaml to {output_dir}")
                    hparams_copied = True
        else:
            print(f"[WARNING] No .cif file found for {substance_id}")
    
    if not hparams_copied:
        print("[WARNING] Could not find hparams.yaml in any boltz_results directory")
    
    if VERBOSE: print(f"Files copied to {output_dir}")

def remove_pkls(boltz_cache, run_cache_prot, VERBOSE):

    run_cache_csv = glob.glob(f'{run_cache_prot}/*.csv')

    run_cache_csv_df = pd.read_csv(run_cache_csv[0])

    for pkl_id in run_cache_csv_df['pkl_id']:
        pkl_file = os.path.join(boltz_cache, 'mols', f'{pkl_id}.pkl')
        if os.path.exists(pkl_file):
            os.remove(pkl_file)
    
    if VERBOSE: print(f"Deleted all generated pkls.")


def get_protein_stats(pred_df, protein=None):
    """
    Get statistics about predictions for each protein.
    Useful for monitoring progress across different proteins.
    
    Args:
        pred_df: predictions DataFrame (long format with 'protein' column)
        protein: specific protein to analyze (None = all proteins)
    
    Returns:
        DataFrame with stats per protein
    """
    if 'protein' not in pred_df.columns:
        print("Warning: 'protein' column not found. Cannot generate protein stats.")
        return None
    
    if protein:
        pred_df = pred_df[pred_df['protein'] == protein]
    
    stats = pred_df.groupby('protein').agg({
        'substance_id': 'count',  # Total rows
        'pred_log10ic50': lambda x: x.notna().sum(),  # Non-null predictions
        'boltz_runID': 'nunique'  # Unique runs
    }).rename(columns={
        'substance_id': 'total_rows',
        'pred_log10ic50': 'valid_predictions',
        'boltz_runID': 'unique_runs'
    })
    
    return stats