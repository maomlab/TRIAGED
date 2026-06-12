##########################
# Plotting
##########################
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ---------------------------------------------------------------------------
# Model classification helpers
# ---------------------------------------------------------------------------

# Model types that support SHAP / variable importance in H2O.
# DeepLearning and StackedEnsemble of DL models are excluded.
_EXPLAINABLE_PREFIXES = (
    'GBM',
    'XGBoost',
    'DRF',
    'RandomForest',
    'XRT',
    'GLM',
)

def is_explainable(model_id: str) -> bool:
    """
    Return True if the model type supports SHAP / variable importance.

    DeepLearning models and StackedEnsembles whose base learners are purely
    DeepLearning are excluded.  All tree-based and linear models are included.
    """
    mid = model_id.upper()
    if 'DEEPLEARNING' in mid:
        return False
    # StackedEnsemble may wrap DL; treat conservatively as explainable
    # because H2O's varimp_plot works for most SE combinations.
    return True
    
def plot_performance(results_df, target, output_dir, csv_name, topN=0.2):
    fig = plt.figure(figsize=(20, 6))
    fig.patch.set_facecolor('#0f0f0f')
    gs = gridspec.GridSpec(1, 3, figure=fig, wspace=0.35)

    # --- leaderboard ---
    ax1 = fig.add_subplot(gs[0])
    ax1.set_facecolor('#1a1a1a')
    models      = results_df['model_id'].str.split('_').str[:3].str.join('_')
    metric_vals = results_df['spearman_rho'].fillna(0)

    # Colour: cyan = best overall, orange = best explainable, blue = rest
    best_overall_val = metric_vals.max()
    colors = []
    for mid, v in zip(results_df['model_id'], metric_vals):
        if v == best_overall_val:
            colors.append('#00e5ff')   # best overall (any type)
        elif is_explainable(mid):
            colors.append('#2979ff')   # explainable model
        else:
            colors.append('#ff6d00')   # non-explainable (e.g. DeepLearning)

    bars = ax1.barh(models, metric_vals, color=colors, edgecolor='none', height=0.6)
    ax1.set_xlabel('Spearman ρ', color='#aaaaaa', fontsize=10)
    ax1.set_title(f'Model Leaderboard\n{target}', color='white', fontsize=11, pad=10)
    ax1.tick_params(colors='#aaaaaa', labelsize=8)
    ax1.spines[:].set_color('#333333')
    ax1.set_xlim(-1, 1)
    ax1.axvline(0, color='#ff4081', linewidth=1.2, linestyle='--', label='random (ρ=0)', alpha=0.8)

    # Legend patches
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#00e5ff', label='best overall'),
        Patch(facecolor='#2979ff', label='explainable'),
        Patch(facecolor='#ff6d00', label='non-explainable (no SHAP)'),
        Patch(facecolor='#ff4081', label='random (ρ=0)'),
    ]
    ax1.legend(handles=legend_elements, fontsize=7, labelcolor='white',
               facecolor='#1a1a1a', edgecolor='#333333', loc='lower right')

    for bar, val in zip(bars, metric_vals):
        ax1.text(val + 0.01, bar.get_y() + bar.get_height() / 2,
                 f'{val:.3f}', va='center', color='white', fontsize=8)

    # --- predicted (x) vs actual (y) — full test set ---
    ax2 = fig.add_subplot(gs[1])
    ax2.set_facecolor('#1a1a1a')
    best     = results_df.iloc[0]
    y_pred   = best.get('y_pred')
    y_actual = best.get('y_actual')

    if y_pred is not None and y_actual is not None:
        y_pred   = np.array(y_pred)
        y_actual = np.array(y_actual)
        ax2.scatter(y_pred, y_actual, color='#00e5ff', alpha=0.4, s=30, edgecolors='none', label='all')
        lim = [min(y_actual.min(), y_pred.min()) - 0.2,
               max(y_actual.max(), y_pred.max()) + 0.2]
        ax2.plot(lim, lim, color='#ff4081', linewidth=1, linestyle='--', label='y=x')
        ax2.set_xlim(lim); ax2.set_ylim(lim)
        ax2.set_xlabel('Predicted', color='#aaaaaa', fontsize=10)
        ax2.set_ylabel('Actual', color='#aaaaaa', fontsize=10)

        rho_ci_low  = best.get('spearman_rho_ci_low', float('nan'))
        rho_ci_high = best.get('spearman_rho_ci_high', float('nan'))

        ax2.set_title(
            f'All Predictions\n'
            f'ρ={best["spearman_rho"]:.3f} [{rho_ci_low:.2f}, {rho_ci_high:.2f}]  '
            f'R²={best["r2"]:.3f}  RMSE={best["rmse"]:.3f}\n'
            f'EF10%={best["ef10"]:.3f}  PR-AUC={best["pr_auc"]:.3f}',
            color='white', fontsize=10, pad=10
        )
        ax2.tick_params(colors='#aaaaaa', labelsize=8)
        ax2.spines[:].set_color('#333333')
        ax2.legend(fontsize=8, labelcolor='white', facecolor='#1a1a1a', edgecolor='#333333')

        # --- predicted (x) vs actual (y) — top-N ---
        ax3 = fig.add_subplot(gs[2])
        ax3.set_facecolor('#1a1a1a')

        n_top          = max(1, int(np.ceil(topN * len(y_actual))))
        pred_order     = np.argsort(y_pred)
        actual_order   = np.argsort(y_actual)
        top_pred_idx   = set(pred_order[:n_top])
        top_actual_idx = set(actual_order[:n_top])
        overlap_idx    = top_pred_idx & top_actual_idx
        pred_only_idx  = top_pred_idx   - top_actual_idx
        truth_only_idx = top_actual_idx - top_pred_idx

        if overlap_idx:
            idx = list(overlap_idx)
            ax3.scatter(y_pred[idx], y_actual[idx], color='#ba68c8', alpha=0.8, s=50,
                        edgecolors='none', label=f'overlap (n={len(idx)})')
        if truth_only_idx:
            idx = list(truth_only_idx)
            ax3.scatter(y_pred[idx], y_actual[idx], color='#ff4081', alpha=0.8, s=50,
                        edgecolors='none', label=f'truth-only (n={len(idx)})')
        if pred_only_idx:
            idx = list(pred_only_idx)
            ax3.scatter(y_pred[idx], y_actual[idx], color='#00e5ff', alpha=0.8, s=50,
                        edgecolors='none', label=f'pred-only (n={len(idx)})')

        ax3.plot(lim, lim, color='#888888', linewidth=1, linestyle='--')
        ax3.set_xlim(lim); ax3.set_ylim(lim)
        ax3.set_xlabel('Predicted', color='#aaaaaa', fontsize=10)
        ax3.set_ylabel('Actual', color='#aaaaaa', fontsize=10)

        rho_top  = best.get('topN_spearman_rho', float('nan'))
        pval_top = best.get('topN_spearman_pval', float('nan'))
        overlap  = best.get('topN_overlap', 0)
        ax3.set_title(
            f'Top {int(topN*100)}% Predictions\n'
            f'ρ(top)={rho_top:.3f} (p={pval_top:.2e})\n'
            f'Overlap: {overlap}/{n_top}',
            color='white', fontsize=10, pad=10
        )
        ax3.tick_params(colors='#aaaaaa', labelsize=8)
        ax3.spines[:].set_color('#333333')
        ax3.legend(fontsize=7, labelcolor='white', facecolor='#1a1a1a', edgecolor='#333333')
    else:
        ax2.text(0.5, 0.5, 'No prediction data', ha='center',
                 va='center', color='#aaaaaa', transform=ax2.transAxes)

    target_slug = target.replace(' ', '_').replace('(', '').replace(')', '').replace('/', '')
    fpath = os.path.join(output_dir, f"{csv_name}_{target_slug}_performance.png")
    plt.savefig(fpath, dpi=150, bbox_inches='tight', facecolor='#0f0f0f')
    plt.close()
    print(f"  Plot saved: {fpath}")

def plot_cross_csv_summary(all_results, target, output_dir, baselines=None):
    """Summary bar chart comparing best metrics across all CSVs with baselines."""
    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    fig.patch.set_facecolor('#0f0f0f')

    labels     = [r['csv_name']      for r in all_results]
    rho_vals   = [r['best_spearman'] for r in all_results]
    r2_vals    = [r['best_r2']       for r in all_results]
    ef_vals    = [r['best_ef10']     for r in all_results]
    prauc_vals = [r['best_pr_auc']   for r in all_results]

    if baselines is None:
        baselines = {}
    rho_baseline    = baselines.get('spearman', 0)
    r2_baseline     = baselines.get('r2', 0)
    ef_baseline     = baselines.get('ef10', 1)
    prauc_baseline  = baselines.get('pr_auc', None)

    # boltz baseline (from per-csv runs)
    boltz_baselines = [r.get('baseline_spearman', float('nan')) for r in all_results]
    valid_boltz     = [b for b in boltz_baselines if not (isinstance(b, float) and np.isnan(b))]
    mean_boltz      = np.mean(valid_boltz) if valid_boltz else None

    panels = [
        (axes[0, 0], rho_vals,   'Best Spearman ρ', (-1, 1),
         rho_baseline,  'random (ρ=0)', mean_boltz, 'Boltz baseline'),
        (axes[0, 1], r2_vals,    'Best R²',         (min(0, min(r2_vals) - 0.1), 1),
         r2_baseline,   'predict mean (R²=0)', None, None),
        (axes[1, 0], ef_vals,    'Best EF10%',      (0, max(max(ef_vals), 1.5) * 1.2),
         ef_baseline,   'random (EF=1)', None, None),
        (axes[1, 1], prauc_vals, 'Best PR-AUC',     (0, 1),
         prauc_baseline,
         f'active fraction ({prauc_baseline:.2f})' if prauc_baseline else None,
         None, None),
    ]

    for ax, vals, ylabel, ylim, baseline, baseline_label, extra_baseline, extra_label in panels:
        ax.set_facecolor('#1a1a1a')
        clean_vals = [v if not (isinstance(v, float) and np.isnan(v)) else 0 for v in vals]
        max_val    = max(clean_vals) if clean_vals else 0
        colors     = ['#00e5ff' if v == max_val and v != 0 else '#2979ff' for v in clean_vals]
        bars       = ax.bar(range(len(labels)), clean_vals, color=colors, edgecolor='none', width=0.6)

        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8, color='#aaaaaa')
        ax.set_ylabel(ylabel, color='#aaaaaa', fontsize=10)
        ax.set_title(f'{ylabel} per Feature Subset\n{target}', color='white', fontsize=11, pad=10)
        ax.tick_params(colors='#aaaaaa')
        ax.spines[:].set_color('#333333')
        ax.set_ylim(*ylim)

        if baseline is not None:
            ax.axhline(baseline, color='#ff4081', linewidth=1.2, linestyle='--',
                       label=baseline_label, alpha=0.8)

        if extra_baseline is not None:
            ax.axhline(extra_baseline, color='#ffc107', linewidth=1.2, linestyle=':',
                       label=f'{extra_label} (ρ={extra_baseline:.3f})', alpha=0.9)

        if baseline is not None or extra_baseline is not None:
            ax.legend(fontsize=7, labelcolor='white',
                      facecolor='#1a1a1a', edgecolor='#333333', loc='upper right')

        ax.axhline(0, color='#333333', linewidth=0.5)
        for bar, val in zip(bars, clean_vals):
            ax.text(bar.get_x() + bar.get_width() / 2, val + (ylim[1] * 0.01),
                    f'{val:.3f}', ha='center', color='white', fontsize=8)

    target_slug = target.replace(' ', '_').replace('(', '').replace(')', '').replace('/', '')
    fpath = os.path.join(output_dir, f"summary_{target_slug}.png")
    plt.tight_layout()
    plt.savefig(fpath, dpi=150, bbox_inches='tight', facecolor='#0f0f0f')
    plt.close()
    print(f"Summary plot saved: {fpath}")
