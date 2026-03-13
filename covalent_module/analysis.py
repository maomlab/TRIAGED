"""
Simple Boltz Analysis 

Inputs:
1. RECORDS directory - Contains experimental_readouts.csv and predictions.csv
2. CIF directory - Structure files

Outputs: metrics, plots, PLIP fingerprints

Usage: python simple_analysis.py --json_file config.json

JSON format:
{
    "RECORDS": "/path/to/records_dir",
    "CIF_DIR": "main_preds",
    "OUTPUT_DIR": "results",
    "PROTEIN_NAME": "tgcpl",
    "SCORE_COL": "Binding Probability",
    "TOP_N": 0.1
}
"""

import os
import json
import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from BoltzCov.analysis.analysis_utils import (
    process_invitro,
    calculate_metrics,
    affinity_scatter,
    plot_curves
)
from BoltzCov.analysis.analyze_boltz_preds import topN_affinity_scatter 

def read_json_config(json_file):
    with open(json_file, 'r') as f:
        config = json.load(f)
    
    for key in config:
        if isinstance(config[key], str):
            config[key] = os.path.expandvars(config[key])
    
    return config


def validate_config(config):
    required = ['RECORDS', 'OUTPUT_DIR']
    missing = [k for k in required if k not in config]
    
    if missing:
        raise ValueError(f"Missing: {', '.join(missing)}")
    
    records_dir = config['RECORDS']
    if not os.path.isdir(records_dir):
        raise FileNotFoundError(f"Records directory not found: {records_dir}")
    
    experimental_csv = os.path.join(records_dir, 'experiment_readouts.csv')
    predictions_csv = os.path.join(records_dir, 'predictions.csv')
    
    if not os.path.exists(experimental_csv):
        raise FileNotFoundError(f"Not found: {experimental_csv}")
    if not os.path.exists(predictions_csv):
        raise FileNotFoundError(f"Not found: {predictions_csv}")
    if 'CIF_DIR' in config and not os.path.isdir(config['CIF_DIR']):
        raise FileNotFoundError(f"Not found: {config['CIF_DIR']}")

def average_replicates(df_predictions, protein):
    """Average replicate predictions per compound for a specific protein."""
    if protein == 'tgcpl':
        protein = '3F75'
    elif protein == 'hscpl':
        protein = '5MAJ'
    # Filter for specific protein
    df_predictions = df_predictions[df_predictions['protein'] == protein].copy()
    if df_predictions.empty:
        print(f"[WARNING] No predictions found for protein: {protein}")
        return df_predictions

    score_cols = [
        'pred_log10ic50',
        'pred_pic50',
        'binding_probability',
        'confidence_score',
        'ptm',
        'iptm',
        'ligand_iptm',
        'protein_iptm',
        'complex_plddt',
        'complex_iplddt',
        'complex_pde',
        'complex_ipde'
    ]
    
    metadata_cols = ['inchi_key', 'vault_mol_id', 'boltz_runID']
    
    cols_to_average = [c for c in score_cols if c in df_predictions.columns]
    group_cols = ['substance_id', 'protein']

    if not cols_to_average:
        return df_predictions

    n_reps = df_predictions.groupby(group_cols).size()

    print(
        f"[{protein}] Averaging replicates: {len(n_reps)} compound–protein pairs, "
        f"{n_reps.mean():.1f} reps/pair"
    )

    agg_dict = {c: ['mean', 'std'] for c in cols_to_average}
    
    for col in metadata_cols:
        if col in df_predictions.columns:
            agg_dict[col] = 'first'

    df_averaged = (
        df_predictions
        .groupby(group_cols)[list(agg_dict.keys())]
        .agg(agg_dict)
        .reset_index()
    )

    df_averaged.columns = [
        f'{col}_{stat}' if stat else col
        for col, stat in df_averaged.columns
    ]

    df_averaged.rename(
        columns={f'{c}_mean': c for c in cols_to_average},
        inplace=True
    )

    df_averaged['n_replicates'] = df_averaged.set_index(group_cols).index.map(n_reps)

    return df_averaged

def merge_with_experimental(df_predictions, df_experimental, exp_col):
    """Merge predictions with experimental data and label actives/inactives."""
    
    id_col = None
    for col in ['substance_id', 'vault_mol_id']:
        if col in df_predictions.columns and col in df_experimental.columns:
            id_col = col
            break
    
    if id_col is None:
        raise ValueError("No common ID column found")
    
    print(f"Using ID column: {id_col}, Experimental column: {exp_col}")
    
    df_exp_labeled = process_invitro(invitro_df=df_experimental, exp_col=exp_col)
    df_exp_labeled['label'] = df_exp_labeled['is_binder'].astype(int)
    
    n_actives = df_exp_labeled['is_binder'].sum()
    print(f"Actives: {n_actives}, Inactives: {len(df_exp_labeled) - n_actives}")
    
    df_merged = pd.merge(df_exp_labeled, df_predictions, on=id_col, how='inner')
    print(f"Merged: {len(df_merged)} compounds")
    
    return df_merged

def compute_metrics_and_plots(df_merged, score_col, exp_col, topN, run_name, output_dir):
    """Compute metrics and generate plots."""

    print("Computing metrics...")
    metrics, curves = calculate_metrics(
        topN=topN,
        df_truth_pred=df_merged,
        score_col=score_col
    )
    
    for key in ['ROC AUC', 'PR-AUC', 'LogAUC', 'BEDROC', 'Spearman r']:
        if key in metrics:
            print(f"  {key}: {metrics[key]:.4f}")
    
    plots_dir = os.path.join(output_dir, 'plots')
    os.makedirs(plots_dir, exist_ok=True)
    
    try:
        scatter_dict = affinity_scatter(
        df_truth_pred=df_merged,
        run_name=run_name,
        score_col=score_col,
        exp_col=exp_col
    )

        analysis_dict = {
            'plots': {
                'scatter': scatter_dict
            }
        }
        scatter_fig = list(scatter_dict.keys())[0]
        scatter_fig.savefig(
            os.path.join(plots_dir, 'scatter.png'),
            dpi=300,
            bbox_inches='tight'
        )

        topN_affinity_scatter(
            truth_pred_df=df_merged,
            analysis_dict=analysis_dict,
            score_col=score_col,
            topN=topN,
            exp_col=exp_col,
            write_output=plots_dir,
            run_name=run_name
        )   

        plt.close(scatter_fig)
        print("  Saved scatter + topN affinity plot")

    except Exception as e:
        print(f"  Scatter plot failed: {e}")
    
    try:
        curve_plots = plot_curves(run_name=run_name, curves=curves, metrics=metrics)
        for curve_name, (fig, ax) in curve_plots.items():
            fig.savefig(os.path.join(plots_dir, f'{curve_name}.png'), dpi=300, bbox_inches='tight')
            plt.close(fig)
        print(f"  Saved {len(curve_plots)} curve plots")
    except Exception as e:
        print(f"  Curve plots failed: {e}")
    
    return metrics, curves


def run_plip_analysis(cif_dir, output_dir, protein_name, records, COVALENT, receptor_type='protein'):
    """Run PLIP fingerprinting."""

    if protein_name == 'tgcpl':
        protein_name='3F75'
    elif protein_name == 'hscpl':
        protein_name='5MAJ'

    print("Running PLIP analysis...")
    
    try:
        from BoltzCov.analysis import run_plip
    except ImportError:
        print("run_plip_batch_replicates.py not found, skipping PLIP")
        return
    
    fp_dir = os.path.join(output_dir, 'fingerprints')
    os.makedirs(fp_dir, exist_ok=True)
    
    plip_args = argparse.Namespace(
        directory=cif_dir,
        outdir=fp_dir,
        receptor_type=receptor_type,
        verbose=False,
        csv_name=f"{protein_name}_ifps",
        selection_method='first', 
        protein_name=protein_name, 
        records=records, 
        COVALENT=COVALENT
    )

    try:
        errors = run_plip.main(plip_args)
        if errors:
            print(f"  PLIP completed with {len(errors)} errors")
        else:
            print("  PLIP completed successfully")
    except Exception as e:
        print(f"  PLIP failed: {e}")


def save_outputs(df_averaged, df_merged, metrics, output_dir):
    """Save CSV outputs and metrics JSON."""
    
    df_averaged.to_csv(os.path.join(output_dir, 'predictions_averaged.csv'), index=False)
    df_merged.to_csv(os.path.join(output_dir, 'predictions_merged.csv'), index=False)
    
    with open(os.path.join(output_dir, 'metrics.json'), 'w') as f:
        json.dump(metrics, f, indent=2)
    
    with open(os.path.join(output_dir, 'summary.txt'), 'w') as f:
        f.write("Boltz Analysis Summary\n")
        f.write("=" * 60 + "\n\n")
        f.write("Performance Metrics:\n")
        f.write("-" * 40 + "\n")
        for key, value in metrics.items():
            if isinstance(value, (int, float)):
                f.write(f"{key:30s}: {value:.4f}\n")
    
    print(f"Saved outputs to {output_dir}")


def main():
    parser = argparse.ArgumentParser(description="Simple Boltz analysis with 3 inputs")
    parser.add_argument('--json_file', required=True, help='JSON config file')
    args = parser.parse_args()
    
    config = read_json_config(args.json_file)
    validate_config(config)
    
    records_dir = config['RECORDS']
    COVALENT = config['COVALENT']
    experimental_csv = os.path.join(records_dir, 'experiment_readouts.csv')
    predictions_csv = os.path.join(records_dir, 'predictions.csv')
    if not COVALENT: 
        predictions_csv = os.path.join(records_dir, 'noncov_predictions.csv')
    cif_dir = config.get('CIF_DIR')
    output_dir = config['OUTPUT_DIR']
    
    protein_name = config.get('PROTEIN_NAME')
    run_name = config.get('RUN_NAME', 'analysis')
    score_col = config.get('SCORE_COL', 'Pred log10(IC50)')
    topN = config.get('TOP_N', 0.1)
    
    run_output_dir = os.path.join(output_dir, run_name)
    os.makedirs(run_output_dir, exist_ok=True)

    print("=" * 80)
    print("Boltz Analysis")
    print("=" * 80)
    print(f"Run: {run_name}")
    print(f"Records: {records_dir}")
    print(f"Score: {score_col}")
    print(f"Output: {run_output_dir}\n")
    
    print("Loading predictions...")
    df_predictions = pd.read_csv(predictions_csv)
    
    df_averaged = average_replicates(df_predictions, protein_name)
    df_predictions = df_predictions.dropna(subset=[score_col])

    print("\nLoading experimental data...")
    df_experimental = pd.read_csv(experimental_csv)
    exp_col = f"mean_{protein_name.lower()}_log_ic50 (uM)"
    
    if exp_col not in df_experimental.columns:
        raise KeyError(f"IC50 column not found. Available: {list(df_experimental.columns)}")
    
    df_experimental = df_experimental.dropna(subset=[exp_col])
    df_merged = merge_with_experimental(df_averaged, df_experimental, exp_col)
    
    print("\n" + "=" * 80)
   
    metrics, _ = compute_metrics_and_plots(
        df_merged, score_col, exp_col, topN, run_name, run_output_dir
    )
    
    if cif_dir:
        print("\n" + "=" * 80)
        run_plip_analysis(cif_dir=cif_dir, output_dir=run_output_dir, receptor_type='protein', protein_name=protein_name, records=records_dir, COVALENT=COVALENT)
    
    print("\n" + "=" * 80)
    save_outputs(df_averaged, df_merged, metrics, run_output_dir)
    
    print("\n" + "=" * 80)
    print("Analysis complete")
    print("=" * 80)


if __name__ == "__main__":
    main()