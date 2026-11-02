"""
free_wilson.py

Free-Wilson / SAR decomposition analysis for pre-encoded R-group columns.

Assumes R-group columns are already one-hot encoded binary columns with
naming convention:  <rgroup_prefix>_<substituent_id>
e.g.  R1_72903, R3a_10088, R3b_0, ...

The script detects which prefix each column belongs to, fits an additive
linear model, and reports per-substituent contributions with bootstrap CIs.

Usage
-----
    python free_wilson.py \
        --csv        feature_matrix.csv \
        --targets    "mean_tgcpl_log_ic50 (uM)" \
        --rgroup_prefixes R1 R3a R3b \
        --output_dir fw_results/ \
        --split_csv  butina_splits.csv \
        --alpha      0.0 \
        --n_bootstrap 1000

The CSV must have substance_id as index (or a column named substance_id).

JSON Config alternative
-----------------------
    python free_wilson.py --json_file fw_config.json

{
    "CSV":              "/path/to/feature_matrix.csv",
    "TARGETS":          ["mean_tgcpl_log_ic50 (uM)"],
    "RGROUP_PREFIXES":  ["R1", "R3a", "R3b"],
    "OUTPUT_DIR":       "/path/to/fw_results",
    "SPLIT_CSV":        "/path/to/butina_splits.csv",
    "ALPHA":            0.0,
    "N_BOOTSTRAP":      1000,
    "SEED":             42
}
"""

import os
import json
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from scipy.stats import spearmanr
from sklearn.linear_model import Ridge, LinearRegression
from sklearn.metrics import r2_score
from sklearn.model_selection import LeaveOneOut


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def slug(s: str) -> str:
    return (s.replace(' ', '_').replace('(', '').replace(')', '')
             .replace('/', '').replace('%', 'pct'))


def compute_rmse(y_true, y_pred):
    return float(np.sqrt(np.mean((np.asarray(y_true) - np.asarray(y_pred)) ** 2)))


BG       = '#0f0f0f'
PANEL_BG = '#1a1a1a'
CYAN     = '#00e5ff'
BLUE     = '#2979ff'
PINK     = '#ff4081'
AMBER    = '#ffc107'
GREY     = '#aaaaaa'
BORDER   = '#333333'


def _style_ax(ax, title='', xlabel='', ylabel=''):
    ax.set_facecolor(PANEL_BG)
    ax.spines[:].set_color(BORDER)
    ax.tick_params(colors=GREY, labelsize=8)
    if title:  ax.set_title(title,   color='white', fontsize=11, pad=8)
    if xlabel: ax.set_xlabel(xlabel, color=GREY,    fontsize=9)
    if ylabel: ax.set_ylabel(ylabel, color=GREY,    fontsize=9)


# ---------------------------------------------------------------------------
# Design matrix — already one-hot, just select and validate
# ---------------------------------------------------------------------------

def extract_design_matrix(df: pd.DataFrame, rgroup_prefixes: list[str]):
    """
    Select all binary one-hot columns whose name starts with one of the
    given prefixes (e.g. 'R1', 'R3a', 'R3b').

    Returns
    -------
    X            : np.ndarray  (n_compounds, n_substituents)
    feat_names   : list[str]   column names
    prefix_map   : dict        {col_name -> rgroup_prefix}
    prefix_counts: dict        {prefix -> n_columns}
    """
    feat_names = []
    prefix_map = {}

    for col in df.columns:
        for pfx in rgroup_prefixes:
            # match "R1_*", "R3a_*" etc — prefix must be followed by underscore
            if col.startswith(pfx + '_'):
                feat_names.append(col)
                prefix_map[col] = pfx
                break

    if not feat_names:
        raise ValueError(
            f"No columns found matching prefixes {rgroup_prefixes}. "
            f"Example columns in df: {df.columns[:10].tolist()}"
        )

    # NaN means the compound has no assignment for that substituent column → treat as 0 (absent)
    X_df = df[feat_names].fillna(0.0)
    X    = X_df.values.astype(float)

    # Validate: each row should sum to ~1 per prefix (one-hot sanity check)
    prefix_counts = {}
    for pfx in rgroup_prefixes:
        pfx_cols = [c for c in feat_names if prefix_map[c] == pfx]
        prefix_counts[pfx] = len(pfx_cols)
        if len(pfx_cols) > 0:
            row_sums = X_df[pfx_cols].sum(axis=1)
            n_multi  = (row_sums > 1).sum()
            n_zero   = (row_sums == 0).sum()
            if n_multi > 0:
                print(f"  Warning: {n_multi} compounds have >1 active column for prefix '{pfx}' "
                      f"— not strictly one-hot")

    # Final NaN check — should never trigger after fillna, but fail loudly if it does
    n_nan = np.isnan(X).sum()
    if n_nan > 0:
        raise ValueError(
            f"Design matrix still contains {n_nan} NaN values after fillna(0). "
            f"Check for non-numeric content in R-group columns."
        )

    print(f"  Design matrix: {X.shape[0]} compounds × {X.shape[1]} substituent columns")
    for pfx, cnt in prefix_counts.items():
        print(f"    {pfx}: {cnt} columns")

    return X, feat_names, prefix_map, prefix_counts


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------

def fit_model(X: np.ndarray, y: np.ndarray, alpha: float = 0.0):
    if alpha > 0:
        model = Ridge(alpha=alpha, fit_intercept=True)
    else:
        model = LinearRegression(fit_intercept=True)
    model.fit(X, y)
    return model


def loo_cv(X: np.ndarray, y: np.ndarray, alpha: float = 0.0):
    """Leave-one-out CV. Returns q² and per-sample CV predictions."""
    y_cv = np.zeros_like(y, dtype=float)
    for train_idx, test_idx in LeaveOneOut().split(X):
        m = Ridge(alpha=alpha, fit_intercept=True) if alpha > 0 \
            else LinearRegression(fit_intercept=True)
        m.fit(X[train_idx], y[train_idx])
        y_cv[test_idx] = m.predict(X[test_idx])

    ss_res = np.sum((y - y_cv) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    q2     = 1 - ss_res / ss_tot if ss_tot > 0 else float('nan')
    return q2, y_cv


def held_out_metrics(y_actual, y_pred):
    rho, pval = spearmanr(y_actual, y_pred)
    ss_res = np.sum((y_actual - y_pred) ** 2)
    ss_tot = np.sum((y_actual - y_actual.mean()) ** 2)
    r2     = 1 - ss_res / ss_tot if ss_tot > 0 else float('nan')
    rmse   = compute_rmse(y_actual, y_pred)
    return {'r2': r2, 'rmse': rmse, 'rho': rho, 'rho_pval': pval}


def bootstrap_contributions(X, y, feat_names, alpha=0.0, n_bootstrap=1000, seed=42):
    """Bootstrap 95% CIs for each coefficient."""
    rng   = np.random.default_rng(seed)
    n     = len(y)
    coefs = np.full((n_bootstrap, X.shape[1]), np.nan)

    for i in range(n_bootstrap):
        idx = rng.choice(n, size=n, replace=True)
        Xb, yb = X[idx], y[idx]
        col_var = Xb.var(axis=0)
        if (col_var == 0).any():
            # degenerate bootstrap sample — skip
            continue
        m = Ridge(alpha=alpha, fit_intercept=True) if alpha > 0 \
            else LinearRegression(fit_intercept=True)
        m.fit(Xb, yb)
        coefs[i] = m.coef_

    return pd.DataFrame({
        'feature':   feat_names,
        'boot_mean': np.nanmean(coefs,               axis=0),
        'ci_low':    np.nanpercentile(coefs, 2.5,    axis=0),
        'ci_high':   np.nanpercentile(coefs, 97.5,   axis=0),
    })


# ---------------------------------------------------------------------------
# Results table
# ---------------------------------------------------------------------------

def build_contributions_table(model, feat_names, prefix_map, boot_df):
    """
    One row per one-hot column with: rgroup, substituent_id, contribution, CIs.
    """
    rows = []
    for i, feat in enumerate(feat_names):
        pfx   = prefix_map[feat]
        subst = feat[len(pfx) + 1:]   # strip "R1_" → "72903"
        coef  = float(model.coef_[i])

        br       = boot_df[boot_df['feature'] == feat]
        ci_low   = float(br['ci_low'].values[0])  if len(br) else np.nan
        ci_high  = float(br['ci_high'].values[0]) if len(br) else np.nan
        sig      = bool(ci_low > 0 or ci_high < 0) \
                   if not (np.isnan(ci_low) or np.isnan(ci_high)) else False

        rows.append({
            'feature':      feat,
            'rgroup':       pfx,
            'substituent':  subst,
            'contribution': coef,
            'ci_low':       ci_low,
            'ci_high':      ci_high,
            'significant':  sig,
            'abs_contribution': abs(coef),
        })

    df = pd.DataFrame(rows).sort_values(
        ['rgroup', 'contribution'], ascending=[True, False]
    ).reset_index(drop=True)
    return df


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_contributions(contrib_df, output_dir, target, intercept,
                       train_metrics, cv_metrics):
    """Forest plot of substituent contributions, one panel per R-group."""
    rgroups = contrib_df['rgroup'].unique().tolist()
    n_rg    = len(rgroups)

    fig, axes = plt.subplots(1, n_rg, figsize=(7 * n_rg, max(6, len(contrib_df) * 0.18 + 2)),
                             squeeze=False)
    fig.patch.set_facecolor(BG)

    for ax, rg in zip(axes[0], rgroups):
        sub    = contrib_df[contrib_df['rgroup'] == rg].sort_values('contribution', ascending=True)
        n      = len(sub)
        y_pos  = np.arange(n)
        colors = [CYAN if c <= 0 else PINK for c in sub['contribution']]

        xerr_low  = np.clip(sub['contribution'].values - sub['ci_low'].values,  0, None)
        xerr_high = np.clip(sub['ci_high'].values - sub['contribution'].values, 0, None)

        ax.barh(y_pos, sub['contribution'], color=colors, edgecolor='none',
                height=0.6, alpha=0.85)
        ax.errorbar(sub['contribution'].values, y_pos,
                    xerr=[xerr_low, xerr_high],
                    fmt='none', color='white', linewidth=1, capsize=2, capthick=1)

        # mark significant
        for i, (_, row) in enumerate(sub.iterrows()):
            if row['significant']:
                offset = 0.005 * (sub['contribution'].abs().max() or 1)
                ha     = 'left' if row['contribution'] >= 0 else 'right'
                ax.text(row['contribution'] + (offset if ha == 'left' else -offset),
                        i, '*', color=AMBER, fontsize=9, va='center', ha=ha)

        ax.set_yticks(y_pos)
        ax.set_yticklabels(sub['substituent'], fontsize=6, color=GREY)
        ax.axvline(0, color=GREY, linewidth=1, linestyle='--', alpha=0.5)
        _style_ax(ax, title=f'{rg}  ({n} substituents)',
                  xlabel='Contribution to log IC50')

    from matplotlib.patches import Patch
    legend_els = [
    Patch(facecolor=CYAN, label='Beneficial (↓ IC50)'),   # negative contribution
    Patch(facecolor=PINK, label='Detrimental (↑ IC50)'),  # positive contribution
    ]
    axes[0][0].legend(handles=legend_els, fontsize=8, labelcolor='white',
                      facecolor=PANEL_BG, edgecolor=BORDER)

    stats = (
        f"μ (intercept) = {intercept:.3f}\n"
        f"Train  R²={train_metrics['r2']:.3f}  RMSE={train_metrics['rmse']:.3f}  ρ={train_metrics['rho']:.3f}\n"
        f"LOO-CV q²={cv_metrics['r2']:.3f}   RMSE={cv_metrics['rmse']:.3f}  ρ={cv_metrics['rho']:.3f}\n"
        f"* = 95% CI excludes zero"
    )
    fig.text(0.5, -0.02, stats, ha='center', color=GREY, fontsize=9,
             bbox=dict(boxstyle='round,pad=0.5', facecolor=PANEL_BG, edgecolor=BORDER))

    fig.suptitle(f'Free-Wilson Substituent Contributions\n{target}',
                 color='white', fontsize=13, y=1.01)
    plt.tight_layout()
    fpath = os.path.join(output_dir, f"{slug(target)}_contributions.svg")
    plt.savefig(fpath, bbox_inches='tight', facecolor=BG, format='svg')
    plt.close()
    print(f"  Contributions plot: {fpath}")


def plot_predicted_vs_actual(y_all, y_pred_all, y_train, y_cv_train,
                              split_labels, train_labels,
                              output_dir, target, train_metrics, cv_metrics):
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.patch.set_facecolor(BG)

    all_vals = np.concatenate([y_all, y_pred_all, y_cv_train])
    lim = [all_vals.min() - 0.3, all_vals.max() + 0.3]

    # panel 1 — fitted (all compounds)
    ax = axes[0]
    _style_ax(ax, title=f"Fitted  R²={train_metrics['r2']:.3f}  RMSE={train_metrics['rmse']:.3f}\n{target}",
              xlabel='Predicted log IC50', ylabel='Actual log IC50')
    if split_labels is not None:
        for sp, col, mk in [('train', BLUE, 'o'), ('val', AMBER, 's'), ('test', CYAN, '^')]:
            mask = split_labels == sp
            if mask.any():
                ax.scatter(y_pred_all[mask], y_all[mask], color=col, alpha=0.75,
                           s=45, edgecolors='none', label=sp, marker=mk)
        ax.legend(fontsize=8, labelcolor='white', facecolor=PANEL_BG, edgecolor=BORDER)
    else:
        ax.scatter(y_pred_all, y_all, color=CYAN, alpha=0.6, s=40, edgecolors='none')
    ax.plot(lim, lim, color=PINK, linewidth=1, linestyle='--')
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.text(0.03, 0.96, f"ρ = {train_metrics['rho']:.3f}", transform=ax.transAxes,
            color=GREY, fontsize=9, va='top')

    # panel 2 — LOO-CV (train compounds only)
    ax = axes[1]
    _style_ax(ax, title=f"LOO-CV  q²={cv_metrics['r2']:.3f}  RMSE={cv_metrics['rmse']:.3f}\n{target}",
              xlabel='Predicted log IC50', ylabel='Actual log IC50')
    if train_labels is not None:
        for sp, col, mk in [('train', BLUE, 'o'), ('val', AMBER, 's')]:
            mask = train_labels == sp
            if mask.any():
                ax.scatter(y_cv_train[mask], y_train[mask], color=col, alpha=0.75,
                           s=45, edgecolors='none', label=sp, marker=mk)
        ax.legend(fontsize=8, labelcolor='white', facecolor=PANEL_BG, edgecolor=BORDER)
    else:
        ax.scatter(y_cv_train, y_train, color=CYAN, alpha=0.6, s=40, edgecolors='none')
    ax.plot(lim, lim, color=PINK, linewidth=1, linestyle='--')
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.text(0.03, 0.96, f"ρ = {cv_metrics['rho']:.3f}", transform=ax.transAxes,
            color=GREY, fontsize=9, va='top')

    plt.tight_layout()
    fpath = os.path.join(output_dir, f"{slug(target)}_predicted_vs_actual.svg")
    plt.savefig(fpath, bbox_inches='tight', facecolor=BG)
    plt.close()
    print(f"  Predicted vs actual: {fpath}")


def plot_top_contributions(contrib_df, output_dir, target, top_n=20):
    """
    Bar chart of the top_n most impactful substituents (by |contribution|)
    across all R-group positions — useful when there are hundreds of columns.
    """
    top = (contrib_df.nlargest(top_n // 2, 'contribution')
           ._append(contrib_df.nsmallest(top_n // 2, 'contribution'))
           .sort_values('contribution', ascending=True))

    n      = len(top)
    y_pos  = np.arange(n)
    colors = [CYAN if c <= 0 else PINK for c in top['contribution']]
    labels = [f"[{r['rgroup']}] {r['substituent']}" for _, r in top.iterrows()]

    xerr_low  = np.clip(top['contribution'].values - top['ci_low'].values,  0, None)
    xerr_high = np.clip(top['ci_high'].values - top['contribution'].values, 0, None)

    fig, ax = plt.subplots(figsize=(11, max(5, n * 0.38 + 1.5)))
    fig.patch.set_facecolor(BG)
    _style_ax(ax, title=f'Top ±{top_n // 2} Substituents by Contribution\n{target}',
              xlabel='Contribution to log IC50')

    ax.barh(y_pos, top['contribution'], color=colors, edgecolor='none', height=0.65, alpha=0.85)
    ax.errorbar(top['contribution'].values, y_pos,
                xerr=[xerr_low, xerr_high],
                fmt='none', color='white', linewidth=1.1, capsize=3)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=8, color=GREY)
    ax.axvline(0, color=GREY, linewidth=1, linestyle='--', alpha=0.5)

    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor=CYAN, label='Beneficial (↓ IC50)'),
                        Patch(facecolor=PINK, label='Detrimental (↑ IC50)')],
              fontsize=8, labelcolor='white', facecolor=PANEL_BG, edgecolor=BORDER)

    plt.tight_layout()
    fpath = os.path.join(output_dir, f"{slug(target)}_top_contributions.svg")
    plt.savefig(fpath, bbox_inches='tight', facecolor=BG)
    plt.close()
    print(f"  Top contributions plot: {fpath}")


def plot_heatmap(contrib_df, output_dir, target, max_per_rgroup=30):
    """
    Heatmap: rows = top substituents per R-position, cols = R-group.
    Limits to max_per_rgroup per position to keep it readable.
    """
    rgroups = contrib_df['rgroup'].unique().tolist()
    frames  = []
    for rg in rgroups:
        sub = contrib_df[contrib_df['rgroup'] == rg]
        # take top N by |contribution|
        sub = sub.nlargest(max_per_rgroup, 'abs_contribution')
        sub = sub[['substituent', 'rgroup', 'contribution']].copy()
        frames.append(sub)

    all_sub = pd.concat(frames)
    pivot   = all_sub.pivot_table(
        index='substituent', columns='rgroup', values='contribution', aggfunc='first'
    ).reindex(columns=rgroups)

    pivot['_sort'] = pivot.abs().max(axis=1)
    pivot = pivot.sort_values('_sort', ascending=False).drop(columns='_sort')

    n_rows, n_cols = pivot.shape
    fig_h = max(5, 0.28 * n_rows + 1.5)
    fig, ax = plt.subplots(figsize=(max(6, 2.2 * n_cols + 2), fig_h))
    fig.patch.set_facecolor(BG)
    _style_ax(ax, title=f'Contribution Heatmap (top {max_per_rgroup} per position)\n{target}')

    vals = pivot.values.astype(float)
    vmax = np.nanmax(np.abs(vals))
    im   = ax.imshow(vals, cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='auto')

    ax.set_xticks(range(n_cols)); ax.set_xticklabels(pivot.columns, color='white', fontsize=10)
    ax.set_yticks(range(n_rows)); ax.set_yticklabels(pivot.index, color=GREY, fontsize=6)

    for i in range(n_rows):
        for j in range(n_cols):
            v = vals[i, j]
            if not np.isnan(v):
                ax.text(j, i, f'{v:.2f}', ha='center', va='center',
                        color='white' if abs(v) > vmax * 0.5 else '#222222', fontsize=6)

    cbar = plt.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cbar.ax.tick_params(colors=GREY, labelsize=8)
    cbar.set_label('Contribution', color=GREY, fontsize=9)

    plt.tight_layout()
    fpath = os.path.join(output_dir, f"{slug(target)}_heatmap.svg")
    plt.savefig(fpath, bbox_inches='tight', facecolor=BG)
    plt.close()
    print(f"  Heatmap: {fpath}")


def plot_rgroup_variance(contrib_df, output_dir, target):
    """
    Bar chart: std dev of contributions per R-group position.
    Shows which position contributes most SAR diversity.
    """
    stats = (contrib_df.groupby('rgroup')['contribution']
             .agg(['std', 'mean', 'count'])
             .rename(columns={'std': 'contrib_std', 'mean': 'contrib_mean', 'count': 'n_substituents'})
             .reset_index()
             .sort_values('contrib_std', ascending=False))

    fig, ax = plt.subplots(figsize=(8, 4))
    fig.patch.set_facecolor(BG)
    _style_ax(ax, title=f'SAR Diversity per R-Group Position\n{target}',
              xlabel='R-Group Position', ylabel='Std Dev of Contributions')

    bars = ax.bar(stats['rgroup'], stats['contrib_std'],
                  color=BLUE, edgecolor='none', width=0.6)
    for bar, row in zip(bars, stats.itertuples()):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.002,
                f'n={row.n_substituents}', ha='center', color=GREY, fontsize=8)

    plt.tight_layout()
    fpath = os.path.join(output_dir, f"{slug(target)}_rgroup_variance.svg")
    plt.savefig(fpath, bbox_inches='tight', facecolor=BG)
    plt.close()
    print(f"  R-group variance plot: {fpath}")


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def free_wilson_analysis(df: pd.DataFrame, target: str,
                          rgroup_prefixes: list[str], output_dir: str,
                          alpha: float = 0.0, n_bootstrap: int = 1000,
                          seed: int = 42, split_col: str = None,
                          run_loo: bool = True, top_n_plot: int = 20):
    os.makedirs(output_dir, exist_ok=True)

    df_fw = df[df[target].notna()].copy()
    print(f"\nFree-Wilson analysis")
    print(f"  Target:    {target}")
    print(f"  Compounds: {len(df_fw)}")
    print(f"  Prefixes:  {rgroup_prefixes}")

    # --- design matrix (already one-hot) ---
    X, feat_names, prefix_map, prefix_counts = extract_design_matrix(df_fw, rgroup_prefixes)
    y = df_fw[target].values.astype(float)

    # --- drop compounds with no R-group assignment ---
    rgroup_mask = X.sum(axis=1) > 0
    n_dropped   = (~rgroup_mask).sum()
    if n_dropped > 0:
        print(f"  Dropping {n_dropped} compounds with no R-group assignment")
        X     = X[rgroup_mask]
        y     = y[rgroup_mask]
        df_fw = df_fw[rgroup_mask]

    print(f"  Compounds after filtering: {len(y)}")

    # --- split labels (after filtering so indices align) ---
    split_labels = df_fw[split_col].values if split_col and split_col in df_fw.columns else None

    # --- train/test split ---
    if split_labels is not None:
        train_mask = pd.Series(split_labels).str.strip().str.lower().isin(['train', 'val']).values
        test_mask  = pd.Series(split_labels).str.strip().str.lower().eq('test').values
        X_train, y_train = X[train_mask], y[train_mask]
        X_test,  y_test  = X[test_mask],  y[test_mask]
        print(f"  train+val={len(y_train)}, test={len(y_test)}")
    else:
        X_train, y_train = X, y
        X_test,  y_test  = X, y

    # --- fit on train ---
    print(f"  Fitting (alpha={alpha}, {X.shape[1]} features)...")
    model = fit_model(X_train, y_train, alpha=alpha)

    # --- predictions ---
    y_pred_all   = model.predict(X)        # all compounds — for saving and plots
    y_pred_train = model.predict(X_train)  # train only — for train metrics
    y_pred_test  = model.predict(X_test)   # test only  — for test metrics

    train_m = held_out_metrics(y_train, y_pred_train)
    test_m  = held_out_metrics(y_test,  y_pred_test)
    print(f"  Train: R²={train_m['r2']:.3f}  RMSE={train_m['rmse']:.3f}  ρ={train_m['rho']:.3f}")
    print(f"  Test:  R²={test_m['r2']:.3f}   RMSE={test_m['rmse']:.3f}  ρ={test_m['rho']:.3f}")

    # --- LOO-CV on train only ---
    if run_loo and len(y_train) <= 500:
        print(f"  LOO-CV ({len(y_train)} folds)...")
        q2, y_cv_train = loo_cv(X_train, y_train, alpha=alpha)
        cv_m = held_out_metrics(y_train, y_cv_train)
        print(f"  LOO-CV: q²={cv_m['r2']:.3f}  RMSE={cv_m['rmse']:.3f}  ρ={cv_m['rho']:.3f}")
    else:
        if len(y_train) > 500:
            print(f"  Skipping LOO-CV (n={len(y_train)} > 500)")
        y_cv_train = y_pred_train
        cv_m       = train_m

    # --- bootstrap CIs ---
    print(f"  Bootstrapping ({n_bootstrap} resamples)...")
    boot_df = bootstrap_contributions(X_train, y_train, feat_names, alpha=alpha,
                                       n_bootstrap=n_bootstrap, seed=seed)

    # --- contributions table ---
    contrib_df = build_contributions_table(model, feat_names, prefix_map, boot_df)

    print(f"\n  Top 5 beneficial substituents (most negative = most potent):")
    print(contrib_df.nsmallest(5, 'contribution')  
        [['rgroup', 'substituent', 'contribution', 'ci_low', 'ci_high', 'significant']]
        .to_string(index=False))

    print(f"\n  Top 5 detrimental substituents (most positive = least potent):")
    print(contrib_df.nlargest(5, 'contribution')  
        [['rgroup', 'substituent', 'contribution', 'ci_low', 'ci_high', 'significant']]
        .to_string(index=False))

    print(f"\n  Significant substituents (95% CI excludes 0): "
          f"{contrib_df['significant'].sum()} / {len(contrib_df)}")

    # --- per-compound predictions ---
    pred_df = pd.DataFrame({
        'substance_id': df_fw.index,
        'actual':       y,
        'fitted':       y_pred_all,
        'residual':     y - y_pred_all,
    })
    if split_labels is not None:
        pred_df['split'] = split_labels

    # --- save ---
    contrib_path = os.path.join(output_dir, f"{slug(target)}_contributions.csv")
    contrib_df.to_csv(contrib_path, index=False)

    pred_path = os.path.join(output_dir, f"{slug(target)}_predictions.csv")
    pred_df.to_csv(pred_path, index=False)

    summary = {
        'target':          target,
        'n_compounds':     int(len(y)),
        'n_features':      int(X.shape[1]),
        'rgroup_prefixes': rgroup_prefixes,
        'prefix_counts':   prefix_counts,
        'alpha':           alpha,
        'intercept':       float(model.intercept_),
        'r2_train':        train_m['r2'],
        'rmse_train':      train_m['rmse'],
        'rho_train':       train_m['rho'],
        'r2_test':         test_m['r2'],
        'rmse_test':       test_m['rmse'],
        'rho_test':        test_m['rho'],
        'q2_loo':          cv_m['r2'],
        'rmse_loo':        cv_m['rmse'],
        'rho_loo':         cv_m['rho'],
        'n_significant':   int(contrib_df['significant'].sum()),
    }
    summary_path = os.path.join(output_dir, f"{slug(target)}_summary.json")
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)

    print(f"\n  Saved: contributions CSV, predictions CSV, summary JSON")
    # in free_wilson_analysis, before plotting
    if split_labels is not None:
        train_split_labels = split_labels[train_mask]
    else:
        train_split_labels = None

    # --- plots ---
    print(f"  Generating plots...")
    plot_predicted_vs_actual(
        y_all          = y,
        y_pred_all     = y_pred_all,
        y_train        = y_train,
        y_cv_train     = y_cv_train,
        split_labels   = split_labels,       # all 322 — for fitted panel
        train_labels   = train_split_labels, # 273 — for LOO panel
        output_dir     = output_dir,
        target         = target,
        train_metrics  = train_m,
        cv_metrics     = cv_m,
    )
    
    plot_contributions(contrib_df, output_dir, target,
                       model.intercept_, train_m, cv_m)
    plot_top_contributions(contrib_df, output_dir, target, top_n=top_n_plot)
    plot_heatmap(contrib_df, output_dir, target)
    plot_rgroup_variance(contrib_df, output_dir, target)

    return contrib_df, pred_df, summary


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(args):
    if args.json_file:
        with open(args.json_file) as f:
            cfg = json.load(f)
        csv_path         = os.path.expandvars(cfg['CSV'])
        targets          = cfg['TARGETS']
        rgroup_prefixes  = cfg['RGROUP_PREFIXES']
        output_dir       = os.path.expandvars(cfg['OUTPUT_DIR'])
        split_csv        = cfg.get('SPLIT_CSV', None)
        alpha            = cfg.get('ALPHA', 0.0)
        n_bootstrap      = cfg.get('N_BOOTSTRAP', 1000)
        seed             = cfg.get('SEED', 42)
        top_n_plot       = cfg.get('TOP_N_PLOT', 20)
        
    else:
        csv_path         = args.csv
        targets          = args.targets
        rgroup_prefixes  = args.rgroup_prefixes
        output_dir       = args.output_dir
        split_csv        = args.split_csv
        alpha            = args.alpha
        n_bootstrap      = args.n_bootstrap
        seed             = args.seed
        top_n_plot       = args.top_n_plot

    df = pd.read_csv(csv_path, index_col='substance_id')
    print(f"Loaded: {csv_path}  ({df.shape[0]} rows, {df.shape[1]} cols)")

    split_col = None
    if split_csv:
        splits    = pd.read_csv(split_csv, index_col='substance_id')
        df        = df.join(splits[['split']], how='inner')
        split_col = 'split'
        print(f"Compounds after inner join: {len(df)}")
        print(df['split'].value_counts())
        print(f"Split labels merged: {df['split'].value_counts().to_dict()}")

    for target in targets:
        if target not in df.columns:
            print(f"Target '{target}' not in CSV — skipping")
            continue
        target_out = os.path.join(output_dir, slug(target))
        free_wilson_analysis(
            df               = df,
            target           = target,
            rgroup_prefixes  = rgroup_prefixes,
            output_dir       = target_out,
            alpha            = alpha,
            n_bootstrap      = n_bootstrap,
            seed             = seed,
            split_col        = split_col,
            top_n_plot       = top_n_plot,
        )

    print("\nDone.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Free-Wilson SAR analysis (pre-encoded one-hot R-groups).")
    parser.add_argument("--json_file",       default=None, help="JSON config (alternative to flags)")
    parser.add_argument("--csv",             default=None, help="Feature matrix CSV")
    parser.add_argument("--targets",         default=None, nargs='+', help="Activity column(s)")
    parser.add_argument("--rgroup_prefixes", default=None, nargs='+',
                        help="R-group column prefixes e.g. R1 R3a R3b")
    parser.add_argument("--output_dir",      default='fw_results', help="Output directory")
    parser.add_argument("--split_csv",       default=None, help="CSV with split labels")
    parser.add_argument("--alpha",           default=0.0,  type=float, help="Ridge alpha (0=OLS)")
    parser.add_argument("--n_bootstrap",     default=1000, type=int)
    parser.add_argument("--seed",            default=42,   type=int)
    parser.add_argument("--top_n_plot",      default=20,   type=int,
                        help="Number of top/bottom substituents in summary plot")
    args = parser.parse_args()

    if args.json_file is None and (args.csv is None or args.targets is None or args.rgroup_prefixes is None):
        parser.error("Provide either --json_file or all of --csv, --targets, --rgroup_prefixes")

    main(args)