import os
import json
import argparse
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score
import plotting
import metrics

import h2o
from h2o.automl import H2OAutoML

_EXPLAINABLE_PREFIXES = (
    'GBM',
    'XGBoost',
    'DRF',
    'RandomForest',
    'XRT',
    'GLM',
)

def is_explainable(model_id: str) -> bool:
    mid = model_id.upper()
    if 'DEEPLEARNING' in mid:
        return False
    return True


##########################
# H2O Training
##########################

def train_h2o(df, target, max_runtime_secs, max_models, nfolds, seed,
              results_dir, csv_name, baseline_df=None, baseline_col='pred_log10ic50'):
    """
    Train H2O AutoML (all model types), evaluate strictly on held-out test set.

    Split safety
    ------------
    The dataframe is split in pandas FIRST into separate train and test frames,
    then each is converted to H2OFrame independently.  This guarantees that
    eval_frame / actuals never contain any training rows, regardless of how H2O
    internally reorders or filters rows during frame conversion.

    Model selection
    ---------------
    * Best overall  — highest Spearman rho on test set, any model type. Saved to disk.
    * Best explainable — highest Spearman rho among non-DeepLearning models.
      SHAP / variable-importance plots generated only for this model.
    """
    all_targets = ['mean_tgcpl_log_ic50 (uM)', 'mean_hscpl_log_ic50 (uM)', 'selectivity']
    drop_cols   = [c for c in all_targets if c != target and c in df.columns]

    id_patterns = ['inchi', 'smiles', 'name', 'synonyms', 'vault_mol_id']
    id_cols     = [c for c in df.columns if any(p in c.lower() for p in id_patterns)]
    if id_cols:
        print(f"  Dropping identifier columns: {id_cols}")
        drop_cols += id_cols

    df_model = df.drop(columns=drop_cols).copy()
    df_model = df_model[df_model[target].notna()]

    if df_model.shape[0] < 10:
        print(f"  Skipping {target} — too few samples ({df_model.shape[0]})")
        return None, None, None

    # -----------------------------------------------------------------------
    # Split pandas-side FIRST, then convert each split to H2OFrame separately.
    # This is the only safe way to guarantee eval_frame contains only test rows.
    # -----------------------------------------------------------------------
    if 'split' in df_model.columns:
        df_model['fold_column'] = df_model['split'].map({'train': 0, 'val': 0, 'test': -1})
        df_model = df_model.drop(columns=['split', 'cluster_id'], errors='ignore')

        train_df = df_model[df_model['fold_column'] == 0].drop(columns=['fold_column']).copy()
        test_df  = df_model[df_model['fold_column'] == -1].drop(columns=['fold_column']).copy()

        if len(test_df) == 0:
            print(f"  Skipping {target} — no test rows after split")
            return None, None, None

        # Ground-truth actuals come from pandas — no H2O row-order dependency
        actuals   = test_df[target].values
        test_mask = df_model['fold_column'] == -1   # boolean Series, indexed by substance_id

        print(f"  train+val={len(train_df)}, test={len(test_df)}")

        # Convert splits independently so H2O row order matches pandas row order
        train_val_hf = h2o.H2OFrame(train_df.reset_index())
        test_hf      = h2o.H2OFrame(test_df.reset_index())

        train_val_hf = train_val_hf.drop('substance_id') \
            if 'substance_id' in train_val_hf.columns else train_val_hf
        test_hf = test_hf.drop('substance_id') \
            if 'substance_id' in test_hf.columns else test_hf

        train_val_hf[target] = train_val_hf[target].asnumeric()
        test_hf[target]      = test_hf[target].asnumeric()

        features   = [c for c in train_val_hf.columns if c != target]
        eval_frame = test_hf   # predictions ONLY ever run on this

        aml = H2OAutoML(
            max_runtime_secs  = max_runtime_secs,
            max_models        = max_models,
            nfolds            = nfolds,
            seed              = seed,
            sort_metric       = 'RMSE',
            verbosity         = None,
        )
        aml.train(
            x                 = features,
            y                 = target,
            training_frame    = train_val_hf,
            leaderboard_frame = test_hf,
        )

    else:
        # No split column — warn loudly, evaluate on full dataset (no held-out set)
        print(f"  WARNING: no 'split' column found — evaluating on FULL dataset (no held-out test).")
        df_model  = df_model.drop(columns=['cluster_id'], errors='ignore')

        full_hf   = h2o.H2OFrame(df_model.reset_index())
        full_hf   = full_hf.drop('substance_id') if 'substance_id' in full_hf.columns else full_hf
        full_hf[target] = full_hf[target].asnumeric()

        actuals    = df_model[target].values
        test_mask  = np.ones(len(df_model), dtype=bool)
        features   = [c for c in full_hf.columns if c != target]
        eval_frame = full_hf

        aml = H2OAutoML(
            max_runtime_secs = max_runtime_secs,
            max_models       = max_models,
            nfolds           = nfolds,
            seed             = seed,
            sort_metric      = 'RMSE',
            verbosity        = None,
        )
        aml.train(x=features, y=target, training_frame=full_hf)

    lb = aml.leaderboard.as_data_frame()

    # -----------------------------------------------------------------------
    # Compute test-set metrics for every model
    # eval_frame is always the held-out test H2OFrame from here on
    # -----------------------------------------------------------------------
    rho_list, pval_list, r2_list, ef_list, prauc_list, rmse_list = [], [], [], [], [], []
    rho_top_list, pval_top_list, overlap_list = [], [], []
    pred_dict       = {}
    explainable_col = []

    print(f"  Computing metrics for {len(lb)} models on test set (n={len(actuals)})...")
    for model_id in lb['model_id']:
        explainable_col.append(metrics.is_explainable(model_id))
        try:
            model = h2o.get_model(model_id)
            preds = model.predict(eval_frame).as_data_frame()['predict'].values

            # Sanity check — lengths must match or metrics are meaningless
            if len(preds) != len(actuals):
                raise ValueError(
                    f"Prediction length {len(preds)} != actuals length {len(actuals)}. "
                    f"eval_frame may contain non-test rows."
                )

            m = metrics.compute_metrics(actuals, preds, topN=0.2)
            rho_list.append(m['spearman_rho'])
            pval_list.append(m['spearman_pval'])
            r2_list.append(m['r2'])
            rmse_list.append(m['rmse'])
            ef_list.append(m['ef10'])
            prauc_list.append(m['pr_auc'])
            rho_top_list.append(m['topN_spearman_rho'])
            pval_top_list.append(m['topN_spearman_pval'])
            overlap_list.append(m['topN_overlap'])
            pred_dict[model_id] = preds

        except Exception as e:
            print(f"  Warning: could not get predictions for {model_id}: {e}")
            for lst in [rho_list, pval_list, r2_list, rmse_list, ef_list, prauc_list,
                        rho_top_list, pval_top_list, overlap_list]:
                lst.append(float('nan'))

    lb['spearman_rho']       = rho_list
    lb['spearman_pval']      = pval_list
    lb['r2']                 = r2_list
    lb['rmse']               = rmse_list
    lb['ef10']               = ef_list
    lb['pr_auc']             = prauc_list
    lb['topN_spearman_rho']  = rho_top_list
    lb['topN_spearman_pval'] = pval_top_list
    lb['topN_overlap']       = overlap_list
    lb['is_explainable']     = explainable_col

    lb = lb.sort_values('spearman_rho', ascending=False).reset_index(drop=True)

    # --- Baseline comparison (test set only via test_mask) ---
    df_for_baseline = df_model.copy()
    baseline = metrics.compute_baseline_spearman(
        df_for_baseline, test_mask, target, baseline_df, baseline_col
    )
    if baseline:
        lb.loc[0, 'baseline_spearman']      = baseline['spearman_rho']
        lb.loc[0, 'baseline_spearman_pval'] = baseline['spearman_pval']

    # -----------------------------------------------------------------------
    # Best OVERALL model (any type)
    # -----------------------------------------------------------------------
    best_overall_idx = 0
    best_overall_id  = lb['model_id'].iloc[best_overall_idx]
    best_overall     = h2o.get_model(best_overall_id)
    print(f"  Best overall model (any type): {best_overall_id}  "
          f"rho={lb['spearman_rho'].iloc[0]:.3f}")

    target_slug = target.replace(' ', '_').replace('(', '').replace(')', '').replace('/', '')
    save_path   = os.path.join(results_dir, f"{csv_name}_{target_slug}_best_model")
    h2o.save_model(model=best_overall, path=save_path, force=True)
    print(f"  Best model saved: {save_path}")

    print(f"  Running bootstrap for best overall model (1000 resamples)...")
    boot = metrics.bootstrap_metrics(actuals, pred_dict[best_overall_id], n_bootstrap=1000)
    for k in ['spearman_rho', 'r2', 'rmse', 'ef10', 'pr_auc']:
        lb.loc[best_overall_idx, f'{k}_ci_low']  = boot[k]['ci_low']
        lb.loc[best_overall_idx, f'{k}_ci_high'] = boot[k]['ci_high']

    print(f"  Bootstrap CIs (95%):")
    print(f"    rho:    [{boot['spearman_rho']['ci_low']:.3f}, {boot['spearman_rho']['ci_high']:.3f}]")
    print(f"    R2:     [{boot['r2']['ci_low']:.3f}, {boot['r2']['ci_high']:.3f}]")
    print(f"    RMSE:   [{boot['rmse']['ci_low']:.3f}, {boot['rmse']['ci_high']:.3f}]")
    print(f"    EF10%:  [{boot['ef10']['ci_low']:.3f}, {boot['ef10']['ci_high']:.3f}]")
    print(f"    PR-AUC: [{boot['pr_auc']['ci_low']:.3f}, {boot['pr_auc']['ci_high']:.3f}]")

    # -----------------------------------------------------------------------
    # Best EXPLAINABLE model — SHAP only here
    # -----------------------------------------------------------------------
    explainable_mask = lb['is_explainable']
    if explainable_mask.any():
        best_explain_idx = lb[explainable_mask].index[0]
        best_explain_id  = lb['model_id'].iloc[best_explain_idx]
        best_explain     = h2o.get_model(best_explain_id)
        print(f"  Best explainable model (SHAP target): {best_explain_id}  "
              f"rho={lb['spearman_rho'].iloc[best_explain_idx]:.3f}")
        lb.loc[best_explain_idx, 'shap_model'] = True
    else:
        print("  Warning: no explainable models found — skipping SHAP.")
        best_explain     = None
        best_explain_id  = None
        best_explain_idx = None

    # -----------------------------------------------------------------------
    # SHAP & explainability — test frame only, best explainable model only
    # -----------------------------------------------------------------------
    explain_dir = os.path.join(results_dir, f"{csv_name}_{target_slug}_explain")
    os.makedirs(explain_dir, exist_ok=True)

    if best_explain is not None:
        print(f"  Generating SHAP/explainability plots for: {best_explain_id}")

        try:
            best_explain.varimp_plot(server=True)
            plt.savefig(os.path.join(explain_dir, 'varimp.png'), dpi=150, bbox_inches='tight')
            plt.close()
        except Exception as e:
            print(f"  Warning: varimp_plot failed: {e}")

        try:
            best_explain.shap_summary_plot(eval_frame)
            plt.savefig(os.path.join(explain_dir, 'shap_summary.png'), dpi=150, bbox_inches='tight')
            plt.close()
        except Exception as e:
            print(f"  Warning: shap_summary_plot failed ({type(e).__name__}): {e}")

        try:
            varimp = best_explain.varimp(use_pandas=True)
            if varimp is not None and len(varimp) > 0:
                top_features = [str(f) for f in varimp['variable'].head(5).tolist()]
                best_explain.partial_plot(eval_frame, cols=top_features, server=True)
                plt.savefig(os.path.join(explain_dir, 'partial_dependence.png'),
                            dpi=150, bbox_inches='tight')
                plt.close()
        except Exception as e:
            print(f"  Warning: partial_plot failed: {e}")

        non_explain_ids = lb.loc[~lb['is_explainable'], 'model_id'].tolist()
        if non_explain_ids:
            note_path = os.path.join(explain_dir, 'shap_skipped_models.txt')
            with open(note_path, 'w') as f:
                f.write(
                    "SHAP was NOT generated for the following models because they "
                    "are not SHAP-compatible (e.g. DeepLearning):\n\n"
                )
                for mid in non_explain_ids:
                    rho_val = lb.loc[lb['model_id'] == mid, 'spearman_rho'].values
                    rho_str = f"{rho_val[0]:.3f}" if len(rho_val) else 'N/A'
                    f.write(f"  {mid}  (rho={rho_str})\n")
                f.write(
                    f"\nSHAP was generated for: {best_explain_id}  "
                    f"(best explainable model, rho="
                    f"{lb['spearman_rho'].iloc[best_explain_idx]:.3f})\n"
                )
            print(f"  SHAP-skipped model list saved: {note_path}")
    else:
        print("  Skipping all SHAP plots — no explainable models available.")

    y_pred_best = pred_dict.get(best_overall_id)
    return lb, y_pred_best, actuals


##########################
# Per-CSV Runner
##########################

def run_on_csv(csv_path, targets, max_runtime_secs, max_models, nfolds, seed,
               split_csv=None, baseline_df=None, baseline_col='pred_log10ic50'):
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
        df     = df.join(splits, how='left')
        print(f"  Split labels merged: {df['split'].value_counts().to_dict()}")

    all_metrics = []
    csv_summary = []

    for target in targets:
        if target not in df.columns:
            print(f"  Target '{target}' not in CSV — skipping")
            continue

        print(f"\n  Training for target: {target}")
        lb, y_pred, y_actual = train_h2o(
            df, target, max_runtime_secs, max_models, nfolds, seed,
            results_dir, csv_name, baseline_df=baseline_df, baseline_col=baseline_col
        )

        if lb is None:
            continue

        lb['target']   = target
        lb['csv_name'] = csv_name

        lb_plot = lb.copy()
        lb_plot['y_pred']   = [y_pred.tolist() if y_pred is not None else None] \
                              + [None] * (len(lb_plot) - 1)
        lb_plot['y_actual'] = [y_actual.tolist()] + [None] * (len(lb_plot) - 1)

        all_metrics.append(lb)

        best = lb.iloc[0]
        csv_summary.append({
            'target':            target,
            'best_model':        best['model_id'],
            'best_r2':           best['r2'],
            'best_spearman':     best['spearman_rho'],
            'spearman_pval':     best['spearman_pval'],
            'best_ef10':         best['ef10'],
            'best_pr_auc':       best['pr_auc'],
            'best_rmse':         best['rmse'],
            'baseline_spearman': best.get('baseline_spearman', float('nan')),
        })

        print(f"  Best model: {best['model_id']}")
        print(f"    rho={best['spearman_rho']:.3f} (p={best['spearman_pval']:.2e})  "
              f"R2={best['r2']:.3f}  RMSE={best['rmse']:.3f}")
        print(f"    EF10%={best['ef10']:.3f}  PR-AUC={best['pr_auc']:.3f}")

        plotting.plot_performance(lb_plot, target, results_dir, csv_name)

        target_slug = target.replace(' ', '_').replace('(', '').replace(')', '').replace('/', '')
        save_cols = ['model_id', 'is_explainable', 'shap_model',
                     'rmse', 'mse', 'mae', 'r2',
                     'spearman_rho', 'spearman_pval', 'ef10', 'pr_auc',
                     'topN_spearman_rho', 'topN_spearman_pval', 'topN_overlap',
                     'spearman_rho_ci_low', 'spearman_rho_ci_high',
                     'r2_ci_low', 'r2_ci_high',
                     'rmse_ci_low', 'rmse_ci_high',
                     'ef10_ci_low', 'ef10_ci_high',
                     'pr_auc_ci_low', 'pr_auc_ci_high',
                     'baseline_spearman', 'baseline_spearman_pval',
                     'target', 'csv_name']
        lb_path = os.path.join(results_dir, f"{csv_name}_{target_slug}_leaderboard.csv")
        lb[[c for c in save_cols if c in lb.columns]].to_csv(lb_path, index=False)
        print(f"  Leaderboard saved: {lb_path}")

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
    split_csv        = cfg.get('SPLIT_CSV', None)
    baseline_csv     = cfg.get('BASELINE_CSV', None)
    baseline_col     = cfg.get('BASELINE_COL', 'pred_log10ic50')

    baseline_df = None
    if baseline_csv:
        baseline_csv = os.path.expandvars(baseline_csv)
        baseline_df  = pd.read_csv(baseline_csv, index_col='substance_id')
        print(f"Loaded baseline: {baseline_csv} ({len(baseline_df)} rows)")
        print(f"  Baseline column: {baseline_col}")

    csv_paths = [
        p for p in glob.glob(os.path.join(base_dir, '**', '*.csv'), recursive=True)
        if 'model_results'     not in p
        and 'manifest'         not in p
        and 'ablation_summary' not in p
    ]

    if not csv_paths:
        raise FileNotFoundError(f"No feature subset CSVs found under: {base_dir}")

    print(f"Found {len(csv_paths)} feature subset CSV(s)")
    for p in csv_paths:
        print(f"  {p}")

    h2o.init(nthreads=-1, max_mem_size='8G')
    h2o.no_progress()

    cross_csv = {t: [] for t in targets}

    for csv_path in sorted(csv_paths):
        csv_name    = os.path.splitext(os.path.basename(csv_path))[0]
        csv_summary = run_on_csv(
            csv_path, targets, max_runtime_secs, max_models, nfolds, seed,
            split_csv=split_csv, baseline_df=baseline_df, baseline_col=baseline_col
        )
        for entry in csv_summary:
            cross_csv[entry['target']].append({
                'csv_name':          csv_name,
                'best_model':        entry['best_model'],
                'best_r2':           entry['best_r2'],
                'best_spearman':     entry['best_spearman'],
                'spearman_pval':     entry['spearman_pval'],
                'best_ef10':         entry['best_ef10'],
                'best_pr_auc':       entry['best_pr_auc'],
                'best_rmse':         entry['best_rmse'],
                'baseline_spearman': entry['baseline_spearman'],
            })

    print(f"\nGenerating cross-CSV summary plots...")
    for target, results in cross_csv.items():
        if results:
            sample_csv = sorted(csv_paths)[0]
            df_sample  = pd.read_csv(sample_csv, index_col='substance_id')
            if target in df_sample.columns:
                active_frac = float((df_sample[target].dropna() < 0).mean())
            else:
                active_frac = None
            baselines = {
                'spearman': 0,
                'r2':       0,
                'ef10':     1,
                'pr_auc':   active_frac,
            }
            plotting.plot_cross_csv_summary(results, target, base_dir, baselines=baselines)

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