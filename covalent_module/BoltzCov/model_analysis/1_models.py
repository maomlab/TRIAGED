"""
ablation_study.py

Runs H2O AutoML models on each feature subset CSV produced by feature_selector.py.
For each CSV found in the input directories:
  - Trains models for each target column
  - Computes Spearman rank correlation for ALL models
  - Re-ranks leaderboard by Spearman rho
  - Saves performance metrics and plots to a 'model_results/' subdir
    within the same timestamped directory as the feature CSV

Usage:
------
    python ablation_study.py --json_file ablation_config.json

JSON Config:
------------
{
    "FEATURE_SUBSET_BASE_DIR": "/path/to/feature_subsets",
    "TARGETS": ["mean_tgcpl_log_ic50 (uM)", "mean_hscpl_log_ic50 (uM)", "selectivity"],
    "MAX_RUNTIME_SECS": 120,
    "MAX_MODELS": 10,
    "NFOLDS": 5,
    "SEED": 42
}
"""

import os
import json
import argparse
import glob
from datetime import datetime
from scipy.stats import spearmanr

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import h2o
from h2o.automl import H2OAutoML

from sklearn.metrics import average_precision_score

def compute_enrichment_factor(y_actual, y_pred, threshold=0.0, fraction=0.1):
    """
    y_actual : log IC50 in uM (lower = more active)
    y_pred   : predicted log IC50 in uM (lower = more active)
    threshold: log IC50 cutoff for active (default 0.0 = 1uM)
    fraction : top fraction to consider for EF (default 0.1 = top 10%)
    """
    # active = 1 if log_ic50 < threshold (more potent)
    actives = (y_actual < threshold).astype(int)

    if actives.sum() == 0:
        print("  Warning: no actives found at threshold, skipping EF")
        return float('nan'), float('nan')

    # negate predictions so lower predicted IC50 ranks higher
    scores = -y_pred

    n_total  = len(y_actual)
    n_active = actives.sum()
    n_top    = max(1, int(np.ceil(fraction * n_total)))

    # get indices of top fraction ranked by score descending
    top_idx     = np.argsort(scores)[::-1][:n_top]
    n_active_top = actives[top_idx].sum()

    # EF = (actives in top fraction / top fraction size) / (total actives / total)
    ef = (n_active_top / n_top) / (n_active / n_total)

    # PR AUC
    pr_auc = average_precision_score(actives, scores)

    return ef, pr_auc

##########################
# Plotting
##########################

def plot_performance(results_df, target, output_dir, csv_name):
    """
    Saves two plots per target per CSV:
    1. Model leaderboard bar chart ranked by Spearman rho
    2. Predicted vs actual scatter for best model with R2 and Spearman annotated
    """
    fig = plt.figure(figsize=(14, 6))
    fig.patch.set_facecolor('#0f0f0f')
    gs = gridspec.GridSpec(1, 2, figure=fig, wspace=0.35)

    # --- leaderboard bar chart ranked by spearman ---
    ax1 = fig.add_subplot(gs[0])
    ax1.set_facecolor('#1a1a1a')

    models = results_df['model_id'].str.split('_').str[:3].str.join('_')

    if 'spearman_rho' in results_df.columns:
        metric_vals  = results_df['spearman_rho'].fillna(0)
        metric_label = 'Spearman ρ'
    else:
        metric_vals  = results_df['r2'].fillna(0)
        metric_label = 'R²'

    colors = ['#00e5ff' if v == metric_vals.max() else '#2979ff' for v in metric_vals]

    bars = ax1.barh(models, metric_vals, color=colors, edgecolor='none', height=0.6)
    ax1.set_xlabel(metric_label, color='#aaaaaa', fontsize=10)
    ax1.set_title(f'Model Leaderboard\n{target}', color='white', fontsize=11, pad=10)
    ax1.tick_params(colors='#aaaaaa', labelsize=8)
    ax1.spines[:].set_color('#333333')
    ax1.set_xlim(-1, 1)
    ax1.axvline(0, color='#333333', linewidth=0.5)
    for bar, val in zip(bars, metric_vals):
        ax1.text(val + 0.01, bar.get_y() + bar.get_height() / 2,
                 f'{val:.3f}', va='center', color='white', fontsize=8)

    # --- predicted vs actual scatter for spearman-best model ---
    ax2 = fig.add_subplot(gs[1])
    ax2.set_facecolor('#1a1a1a')

    best     = results_df.iloc[0]
    y_pred   = best.get('y_pred')
    y_actual = best.get('y_actual')

    if y_pred is not None and y_actual is not None:
        ax2.scatter(y_actual, y_pred, color='#00e5ff', alpha=0.6, s=30, edgecolors='none')

        lim = [min(min(y_actual), min(y_pred)) - 0.2,
               max(max(y_actual), max(y_pred)) + 0.2]
        ax2.plot(lim, lim, color='#ff4081', linewidth=1, linestyle='--', label='y=x')
        ax2.set_xlim(lim)
        ax2.set_ylim(lim)
        ax2.set_xlabel('Actual', color='#aaaaaa', fontsize=10)
        ax2.set_ylabel('Predicted', color='#aaaaaa', fontsize=10)

        rho  = best.get('spearman_rho', float('nan'))
        pval = best.get('spearman_pval', float('nan'))
        r2   = best.get('r2', float('nan'))
        rmse = best.get('rmse', float('nan'))

        ax2.set_title(
            f'Best Model: Pred vs Actual\n'
            f'ρ={rho:.3f} (p={pval:.2e})  R²={r2:.3f}  RMSE={rmse:.3f}',
            color='white', fontsize=10, pad=10
        )
        ax2.tick_params(colors='#aaaaaa', labelsize=8)
        ax2.spines[:].set_color('#333333')
        ax2.legend(fontsize=8, labelcolor='white',
                   facecolor='#1a1a1a', edgecolor='#333333')
    else:
        ax2.text(0.5, 0.5, 'No prediction data', ha='center', va='center',
                 color='#aaaaaa', transform=ax2.transAxes)

    target_slug = target.replace(' ', '_').replace('(', '').replace(')', '').replace('/', '')
    fname = f"{csv_name}_{target_slug}_performance.png"
    fpath = os.path.join(output_dir, fname)
    plt.savefig(fpath, dpi=150, bbox_inches='tight', facecolor='#0f0f0f')
    plt.close()
    print(f"  Plot saved: {fpath}")


def plot_cross_csv_summary(all_results, target, output_dir):
    """
    Summary bar chart comparing best Spearman rho and R2 across all CSVs.
    Two side-by-side panels. Saved in the base output dir.
    """
    fig, axes = plt.subplots(1, 2, figsize=(16, 5))
    fig.patch.set_facecolor('#0f0f0f')

    labels   = [r['csv_name'] for r in all_results]
    rho_vals = [r['best_spearman'] for r in all_results]
    r2_vals  = [r['best_r2'] for r in all_results]

    for ax, vals, ylabel, title_metric in zip(
        axes,
        [rho_vals, r2_vals],
        ['Best Spearman ρ', 'Best R²'],
        ['Spearman ρ', 'R²']
    ):
        ax.set_facecolor('#1a1a1a')
        colors = ['#00e5ff' if v == max(vals) else '#2979ff' for v in vals]
        bars = ax.bar(range(len(labels)), vals, color=colors, edgecolor='none', width=0.6)
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8, color='#aaaaaa')
        ax.set_ylabel(ylabel, color='#aaaaaa', fontsize=10)
        ax.set_title(f'{title_metric} per Feature Subset\n{target}',
                     color='white', fontsize=11, pad=10)
        ax.tick_params(colors='#aaaaaa')
        ax.spines[:].set_color('#333333')
        ax.set_ylim(-1 if 'ρ' in ylabel else 0, 1)
        ax.axhline(0, color='#333333', linewidth=0.5)
        for bar, val in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2, val + 0.01,
                    f'{val:.3f}', ha='center', color='white', fontsize=8)

    target_slug = target.replace(' ', '_').replace('(', '').replace(')', '').replace('/', '')
    fname = f"summary_{target_slug}.png"
    fpath = os.path.join(output_dir, fname)
    plt.tight_layout()
    plt.savefig(fpath, dpi=150, bbox_inches='tight', facecolor='#0f0f0f')
    plt.close()
    print(f"Summary plot saved: {fpath}")


##########################
# H2O Model Training
##########################

def train_h2o(df, target, max_runtime_secs, max_models, nfolds, seed):
    # drop other targets to avoid leakage
    all_targets = ['mean_tgcpl_log_ic50 (uM)', 'mean_hscpl_log_ic50 (uM)', 'selectivity']
    drop_cols   = [c for c in all_targets if c != target and c in df.columns]
    df_model    = df.drop(columns=drop_cols).copy()

    # drop rows where target is missing
    df_model = df_model[df_model[target].notna()].reset_index(drop=True)

    if df_model.shape[0] < 10:
        print(f"  Skipping {target} — too few samples ({df_model.shape[0]})")
        return None, None, None

    if 'split' in df_model.columns:
        fold_map = {'train': 0, 'val': 0, 'test': -1}  # train+val combined = 0, test = -1
        df_model['fold_column'] = df_model['split'].map(fold_map)
        df_model = df_model.drop(columns=['split', 'cluster_id'], errors='ignore')

        hf = h2o.H2OFrame(df_model)
        hf[target] = hf[target].asnumeric()
        features = [c for c in hf.columns if c not in [target, 'fold_column']]

        train_val_hf = hf[hf['fold_column'] == 0]
        test_hf      = hf[hf['fold_column'] == -1]

        print(f"  train+val={train_val_hf.shape[0]}, test={test_hf.shape[0]}")

        actuals = df_model.loc[df_model['fold_column'] == -1, target].values

        aml = H2OAutoML(
            max_runtime_secs = max_runtime_secs,
            max_models       = max_models,
            nfolds           = nfolds,  # CV on train+val
            seed             = seed,
            sort_metric      = 'RMSE',
            verbosity        = None,
        )
        aml.train(
            x                 = features,
            y                 = target,
            training_frame    = train_val_hf,
            leaderboard_frame = test_hf,
        )
        eval_frame = test_hf

    else:
        df_model = df_model.drop(columns=['split', 'cluster_id'], errors='ignore')
        hf = h2o.H2OFrame(df_model)
        hf[target] = hf[target].asnumeric()
        features = [c for c in hf.columns if c not in [target, 'fold_column']]
        actuals  = df_model[target].values

        aml = H2OAutoML(
            max_runtime_secs = max_runtime_secs,
            max_models       = max_models,
            nfolds           = nfolds,
            seed             = seed,
            sort_metric      = 'RMSE',
            verbosity        = None,
        )
        aml.train(x=features, y=target, training_frame=hf)
        eval_frame = hf

    lb = aml.leaderboard.as_data_frame()

    rho_list, pval_list, pred_dict = [], [], {}
    print(f"  Computing Spearman for {len(lb)} models...")
    for model_id in lb['model_id']:
        try:
            model = h2o.get_model(model_id)
            preds = model.predict(eval_frame).as_data_frame()['predict'].values
            rho, pval = spearmanr(actuals, preds)
            rho_list.append(rho)
            pval_list.append(pval)
            pred_dict[model_id] = preds
        except Exception as e:
            print(f"  Warning: could not get predictions for {model_id}: {e}")
            rho_list.append(float('nan'))
            pval_list.append(float('nan'))

    lb['spearman_rho']  = rho_list
    lb['spearman_pval'] = pval_list
    lb = lb.sort_values('spearman_rho', ascending=False).reset_index(drop=True)

    best_model_id = lb['model_id'].iloc[0]
    y_pred        = pred_dict.get(best_model_id)

    if y_pred is not None:
        ss_res  = np.sum((actuals - y_pred) ** 2)
        ss_tot  = np.sum((actuals - actuals.mean()) ** 2)
        r2_best = 1 - ss_res / ss_tot if ss_tot > 0 else float('nan')
    else:
        r2_best = float('nan')

    if 'r2' not in lb.columns:
        lb['r2'] = float('nan')
    lb.loc[0, 'r2'] = r2_best

    if y_pred is not None:
        ss_res  = np.sum((actuals - y_pred) ** 2)
        ss_tot  = np.sum((actuals - actuals.mean()) ** 2)
        r2_best = 1 - ss_res / ss_tot if ss_tot > 0 else float('nan')

        # EF and PR AUC on test set
        ef, pr_auc = compute_enrichment_factor(actuals, y_pred, threshold=0.0, fraction=0.1)
        print(f"  EF10%={ef:.3f}  PR-AUC={pr_auc:.3f}")
    else:
        r2_best = float('nan')
        ef, pr_auc = float('nan'), float('nan')

    if 'r2' not in lb.columns:
        lb['r2'] = float('nan')
    lb.loc[0, 'r2']     = r2_best
    lb.loc[0, 'ef10']   = ef
    lb.loc[0, 'pr_auc'] = pr_auc

    return lb, y_pred, actuals

##########################
# Per-CSV Runner
##########################

def run_on_csv(csv_path, targets, max_runtime_secs, max_models, nfolds, seed, split_csv=None):
    """
    Runs H2O AutoML for all targets on a single feature subset CSV.
    Saves results and plots to model_results/ subdir.
    """
    csv_dir     = os.path.dirname(csv_path)
    csv_name    = os.path.splitext(os.path.basename(csv_path))[0]
    results_dir = os.path.join(csv_dir, 'model_results')
    os.makedirs(results_dir, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"CSV: {csv_name}")
    print(f"Results dir: {results_dir}")
    print(f"{'='*60}")

    df = pd.read_csv(csv_path, index_col='substance_id')
    df = df.fillna(0)
    
    if split_csv is not None:
        splits = pd.read_csv(split_csv, index_col='substance_id')[['split', 'cluster_id']]
        df = df.join(splits, how='left')
        print(f"  Split labels merged: {df['split'].value_counts().to_dict()}")

    all_metrics = []
    csv_summary = []

    for target in targets:
        if target not in df.columns:
            print(f"  Target '{target}' not in CSV — skipping")
            continue

        print(f"\n  Training for target: {target}")
        lb, y_pred, y_actual = train_h2o(
            df, target, max_runtime_secs, max_models, nfolds, seed
        )

        if lb is None:
            continue

        lb['target']   = target
        lb['csv_name'] = csv_name

        # attach pred/actual to top row for scatter plot
        lb_plot = lb.copy()
        lb_plot['y_pred']   = [y_pred.tolist() if y_pred is not None else None] \
                              + [None] * (len(lb_plot) - 1)
        lb_plot['y_actual'] = [y_actual.tolist()] + [None] * (len(lb_plot) - 1)

        all_metrics.append(lb)

        best_rho  = lb['spearman_rho'].iloc[0]
        best_pval = lb['spearman_pval'].iloc[0]
        best_r2   = lb['r2'].iloc[0]

        csv_summary.append({
            'target':        target,
            'best_r2':       best_r2,
            'best_spearman': best_rho,
            'spearman_pval': best_pval,
        })

        print(f"  Best model (by ρ) — ρ={best_rho:.3f} (p={best_pval:.2e})  R²={best_r2:.3f}")

        # plot
        plot_performance(lb_plot, target, results_dir, csv_name)

        # save leaderboard
        target_slug = target.replace(' ', '_').replace('(', '').replace(')', '').replace('/', '')
        lb_path     = os.path.join(results_dir, f"{csv_name}_{target_slug}_leaderboard.csv")
        save_cols = ['model_id', 'rmse', 'mse', 'mae', 'r2',
             'spearman_rho', 'spearman_pval', 'ef10', 'pr_auc', 'target', 'csv_name']
        lb[[c for c in save_cols if c in lb.columns]].to_csv(lb_path, index=False)
        print(f"  Leaderboard saved: {lb_path}")

    # save combined metrics for this CSV
    if all_metrics:
        combined      = pd.concat(all_metrics, ignore_index=True)
        combined_path = os.path.join(results_dir, f"{csv_name}_all_metrics.csv")
        combined.to_csv(combined_path, index=False)

    return csv_summary


##########################
# Main
##########################

def main(args):

    with open(args.json_file) as f:
        cfg = json.load(f)

    base_dir         = os.path.expandvars(cfg['FEATURE_SUBSET_BASE_DIR'])
    targets          = cfg['TARGETS']
    max_runtime_secs = cfg.get('MAX_RUNTIME_SECS', 120)
    max_models       = cfg.get('MAX_MODELS', 10)
    nfolds           = cfg.get('NFOLDS', 5)
    seed             = cfg.get('SEED', 42)
    split_csv = cfg.get('SPLIT_CSV', None)

    # find all feature subset CSVs inside timestamped subdirs
    csv_paths = [
        p for p in glob.glob(os.path.join(base_dir, '**', '*.csv'), recursive=True)
        if 'model_results'    not in p
        and 'manifest'        not in p
        and 'ablation_summary' not in p
    ]

    if not csv_paths:
        raise FileNotFoundError(f"No feature subset CSVs found under: {base_dir}")

    print(f"Found {len(csv_paths)} feature subset CSV(s)")
    for p in csv_paths:
        print(f"  {p}")

    # init H2O once
    h2o.init(nthreads=-1, max_mem_size='8G')
    h2o.no_progress()

    # track best metrics per target across all CSVs for summary
    cross_csv = {t: [] for t in targets}

    for csv_path in sorted(csv_paths):
        csv_name    = os.path.splitext(os.path.basename(csv_path))[0]
        csv_summary = run_on_csv(csv_path, targets, max_runtime_secs, max_models, nfolds, seed, split_csv=split_csv)

        for entry in csv_summary:
            cross_csv[entry['target']].append({
                'csv_name':      csv_name,
                'best_r2':       entry['best_r2'],
                'best_spearman': entry['best_spearman'],
                'spearman_pval': entry['spearman_pval'],
            })

    # cross-CSV summary plots and table
    print(f"\nGenerating cross-CSV summary plots...")
    for target, results in cross_csv.items():
        if results:
            plot_cross_csv_summary(results, target, base_dir)

    summary_rows = []
    for target, results in cross_csv.items():
        for r in results:
            summary_rows.append({'target': target, **r})

    if summary_rows:
        summary_df   = pd.DataFrame(summary_rows)
        summary_path = os.path.join(base_dir, 'ablation_summary.csv')
        summary_df.to_csv(summary_path, index=False)
        print(f"Cross-CSV summary saved: {summary_path}")

    h2o.shutdown(prompt=False)
    print("\nDone.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Ablation study: runs H2O AutoML on each feature subset CSV."
    )
    parser.add_argument("--json_file", required=True, help="Path to JSON config.")
    args = parser.parse_args()
    main(args)