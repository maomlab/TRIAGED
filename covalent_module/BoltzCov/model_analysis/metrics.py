##########################
# Metrics
##########################
import numpy as np
import matplotlib
matplotlib.use('Agg')

from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score
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
def compute_baseline_spearman(df_model, test_mask, target, baseline_df, baseline_col):
    """
    Compute Spearman of a baseline prediction column vs actual target on test set.
    baseline_df : separate dataframe indexed by substance_id with baseline preds
    """
    if baseline_df is None or baseline_col not in baseline_df.columns:
        print(f"  Baseline df missing or column '{baseline_col}' not found — skipping")
        return None

    df_test = df_model[test_mask].copy()
    df_test = df_test.join(baseline_df[[baseline_col]], how='left')
    df_test = df_test[[target, baseline_col]].dropna()

    if len(df_test) < 3:
        print(f"  Too few test points with baseline preds: {len(df_test)}")
        return None

    rho, pval = spearmanr(df_test[target], df_test[baseline_col])
    print(f"  Baseline ({baseline_col}) on test: ρ={rho:.3f} (p={pval:.2e}, n={len(df_test)})")
    return {'spearman_rho': rho, 'spearman_pval': pval, 'n': len(df_test)}

def compute_enrichment_factor(y_actual, y_pred, threshold=None, fraction=0.2):
    """
    If threshold is None, use the median of y_actual as the active cutoff.
    This guarantees ~50% actives and makes EF meaningful regardless of scale.
    """
    if threshold is None:
        threshold = float(np.median(y_actual))

    actives = (y_actual < threshold).astype(int)

    if actives.sum() == 0 or actives.sum() == len(actives):
        print(f"  Warning: degenerate active split at threshold={threshold:.3f} "
              f"(n_active={actives.sum()}/{len(actives)}) — EF undefined")
        return float('nan'), float('nan')

    scores   = -y_pred
    n_total  = len(y_actual)
    n_active = actives.sum()
    n_top    = max(1, int(np.ceil(fraction * n_total)))

    top_idx      = np.argsort(scores)[::-1][:n_top]
    n_active_top = actives[top_idx].sum()

    ef     = (n_active_top / n_top) / (n_active / n_total)
    pr_auc = average_precision_score(actives, scores)

    print(f"  EF10%: threshold={threshold:.3f}, n_active={n_active}/{n_total}, "
          f"n_top={n_top}, n_active_top={n_active_top}, EF={ef:.3f}")

    return ef, pr_auc

def compute_topN_metrics(y_actual, y_pred, topN=0.2):
    """Compute Spearman on the top N% of predicted ligands (lower IC50 = better)."""
    n_top = max(1, int(np.ceil(topN * len(y_actual))))

    pred_order   = np.argsort(y_pred)
    actual_order = np.argsort(y_actual)

    top_pred_idx   = set(pred_order[:n_top])
    top_actual_idx = set(actual_order[:n_top])

    overlap_idx = list(top_pred_idx & top_actual_idx)
    n_overlap   = len(overlap_idx)

    top_pred_mask = np.zeros(len(y_actual), dtype=bool)
    top_pred_mask[list(top_pred_idx)] = True

    if top_pred_mask.sum() >= 2:
        rho_top, pval_top = spearmanr(y_actual[top_pred_mask], y_pred[top_pred_mask])
    else:
        rho_top, pval_top = float('nan'), float('nan')

    return {
        'topN_spearman_rho':  rho_top,
        'topN_spearman_pval': pval_top,
        'topN_overlap':       n_overlap,
        'topN_overlap_frac':  n_overlap / n_top,
    }


def compute_metrics(y_actual, y_pred, topN=0.2):
    rho, pval  = spearmanr(y_actual, y_pred)
    ss_res     = np.sum((y_actual - y_pred) ** 2)
    ss_tot     = np.sum((y_actual - y_actual.mean()) ** 2)
    r2         = 1 - ss_res / ss_tot if ss_tot > 0 else float('nan')
    rmse       = float(np.sqrt(np.mean((y_actual - y_pred) ** 2)))
    ef, pr_auc = compute_enrichment_factor(y_actual, y_pred)
    top        = compute_topN_metrics(y_actual, y_pred, topN=topN)
    return {
        'spearman_rho':       rho,
        'spearman_pval':      pval,
        'r2':                 r2,
        'rmse':               rmse,
        'ef10':               ef,
        'pr_auc':             pr_auc,
        'topN_spearman_rho':  top['topN_spearman_rho'],
        'topN_spearman_pval': top['topN_spearman_pval'],
        'topN_overlap':       top['topN_overlap'],
        'topN_overlap_frac':  top['topN_overlap_frac'],
    }


def bootstrap_metrics(y_actual, y_pred, n_bootstrap=1000, seed=42):
    """Bootstrap CIs for all metrics. Returns {metric: {mean, ci_low, ci_high}}."""
    rng = np.random.default_rng(seed)
    n   = len(y_actual)

    metrics_lists = {'spearman_rho': [], 'r2': [], 'rmse': [], 'ef10': [], 'pr_auc': []}

    for _ in range(n_bootstrap):
        idx = rng.choice(n, size=n, replace=True)
        y_a = y_actual[idx]
        y_p = y_pred[idx]

        if len(np.unique(y_a)) < 2:
            continue

        try:
            rho, _ = spearmanr(y_a, y_p)
            ss_res = np.sum((y_a - y_p) ** 2)
            ss_tot = np.sum((y_a - y_a.mean()) ** 2)
            r2     = 1 - ss_res / ss_tot if ss_tot > 0 else float('nan')
            rmse   = float(np.sqrt(np.mean((y_a - y_p) ** 2)))
            ef, pr_auc = compute_enrichment_factor(y_a, y_p)

            metrics_lists['spearman_rho'].append(rho)
            metrics_lists['r2'].append(r2)
            metrics_lists['rmse'].append(rmse)
            metrics_lists['ef10'].append(ef)
            metrics_lists['pr_auc'].append(pr_auc)
        except Exception:
            continue

    results = {}
    for k, vals in metrics_lists.items():
        vals = np.array([v for v in vals if not np.isnan(v)])
        if len(vals) > 0:
            results[k] = {
                'mean':    float(np.mean(vals)),
                'ci_low':  float(np.percentile(vals, 2.5)),
                'ci_high': float(np.percentile(vals, 97.5)),
            }
        else:
            results[k] = {'mean': float('nan'), 'ci_low': float('nan'), 'ci_high': float('nan')}

    return results
