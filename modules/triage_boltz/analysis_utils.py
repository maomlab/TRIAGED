import os, sys, glob
import json
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score, precision_recall_curve, auc as sk_auc
from tqdm import tqdm
import matplotlib  # must import first
import matplotlib.pyplot as plt
matplotlib.use('Agg')
import seaborn as sns # Making graphs 
from plotnine import ggplot, aes, geom_histogram, geom_vline, annotate, theme_minimal, ggtitle

import numpy as np
from scipy import stats
from statsmodels.stats import weightstats as stests
from typing import List, Dict, Optional, Tuple, Union
from logauc_new import compute_logauc_from_dataframe as logAUC_from_dataframe

    
def convert_IC_to_energy(IC: float) -> float:
    """
    From README docs of boltz2:
    Convert IC (IC50) to kcal/mol.
    :param IC: IC50.
    :return: converted kcal/mol estimate.
    """
    return -(6 - IC) * 1.364

def read_boltz_predictions(predictions_dir: str) -> pd.DataFrame:
    """
    Reads prediction JSON files from subdirectories and compiles them into a pandas DataFrame.
    :param predictions_dir: Path to the directory containing subdirectories with JSON files.
    :return: Pandas DataFrame with compiled data.
    """
    data = []
    #print(predictions_dir)
    #print(os.listdir(predictions_dir))
    for compound_name in os.listdir(predictions_dir):
        compound_dir = os.path.join(predictions_dir, compound_name)
        if not os.path.isdir(compound_dir):
            continue

        affinity_file = os.path.join(compound_dir, f"affinity_{compound_name}.json")
        confidence_file = os.path.join(compound_dir, f"confidence_{compound_name}_model_0.json")

        if os.path.exists(affinity_file) and os.path.exists(confidence_file):
            with open(affinity_file, 'r') as af:
                affinity_data = json.load(af)
                affinity_pred_value = affinity_data.get("affinity_pred_value", None)
                affinity_probability_binary = affinity_data.get("affinity_probability_binary", None)
            with open(confidence_file, 'r') as cf:
                confidence_data = json.load(cf)
                confidence_score = confidence_data.get("confidence_score", None)
                ptm = confidence_data.get("ptm", None)
                iptm = confidence_data.get("iptm", None)
                ligand_iptm = confidence_data.get("ligand_iptm", None)
                protein_iptm = confidence_data.get("protein_iptm", None)
                complex_plddt = confidence_data.get("complex_plddt", None)
                complex_iplddt = confidence_data.get("complex_iplddt", None)
                complex_pde = confidence_data.get("complex_pde", None)
                complex_ipde = confidence_data.get("complex_ipde", None)

            energy_value = convert_IC_to_energy(affinity_pred_value) if affinity_pred_value is not None else None

            data.append({
                "compound_ID": compound_name,
                "Affinity Pred Value": affinity_pred_value,
                "Affinity Probability Binary": affinity_probability_binary,
                "Confidence Score": confidence_score,
                "kcal/mol": energy_value,
                "PTM": ptm,
                "IPTM": iptm,
                "Ligand IPTM": ligand_iptm,
                "Protein IPTM": protein_iptm,
                "Complex pLDDT": complex_plddt,
                "Complex iPLDDT": complex_iplddt,
                "Complex PDE": complex_pde,
                "Complex iPDE": complex_ipde
            })

    return pd.DataFrame(data)

#from Bootstrap_TLDR
def enrichment_standard(scores1: Dict[str, List[Union[str, float]]], lig_list: List[str], decoy_list: List[str], metrics: str) -> None:
    mols = scores1.keys()
    ranked_list = [scores1[key] for key in mols]
    ranked_list.sort(key=lambda x: float(x[-1]))

    points = do_roc(ranked_list, lig_list, decoy_list)
    points = interpolate_curve(points)
    auc    = AUC(points)*100
    logauc = logAUC(points)*100

    fig = plt.figure(figsize=(5, 5))
    fig.subplots_adjust(hspace=0.4, wspace=0.4)
    sns.set_style("white")
    sns.set_context("notebook", font_scale=1.0)
    ax1 = fig.add_subplot(1, 1, 1)

    x = np.arange(0,100,0.1)
    if metrics == "AUC":
            ax1.plot(x, x, 'k--')
            ax1.axis([-0.25, 100, 0, 100])
            ax1.set_xlabel(" Decoys Found %")
            ax1.set_ylabel(" Ligands Found %")
            x, y = zip(*points)
            ax1.plot(x, y, linewidth=1.5, label='AUC: %.2f' % auc)
            # ax1.plot(x, y, linewidth=1)
                            
    elif metrics == "logAUC":
            ax1.semilogx(x, x, 'k--')
            ax1.axis([0.1, 100, 0, 100])
            ax1.set_xlabel(" Decoys Found %")
            ax1.set_ylabel(" Ligands Found %")
            x, y = zip(*points)
            ax1.semilogx(x, y, linewidth=1.5, label='logAUC: %.2f' % logauc)
            # ax1.semilogx(x, y, linewidth=1)
    ax1.legend(loc="best")

    fig.suptitle('ROC Plot')
    fig.tight_layout(pad=2.0)
    fig.savefig(f"plot_{metrics}.png")


def compute_vscreen_metrics(
    df_pos: pd.DataFrame, 
    df_neg: pd.DataFrame, 
    df_pred: pd.DataFrame, 
    score_col: str, 
    output_dir: str, 
    compound_col: str = 'compound_ID', 
    return_curves: bool = False, 
    bootstrap: bool = False
) -> Union[Dict[str, Union[float, Dict[str, float]]], Tuple[Dict[str, Union[float, Dict[str, float]]], Dict[str, Tuple[np.ndarray, np.ndarray]]]]:
    # Replace '/' in score_col with 'kcal' for output file names
    sanitized_score_col = score_col.replace('/', 'kcal') if score_col == 'kcal/mol' else score_col

    df_pos = df_pos.copy()
    df_neg = df_neg.copy()
    df_pos['label'] = 1
    df_neg['label'] = 0
    df_truth = pd.concat([df_pos, df_neg], ignore_index=True)
    # Ensure compound_col values are strings and stripped of whitespace for both DataFrames
    df_truth[compound_col] = df_truth[compound_col].astype(str).str.strip()
    df_pred[compound_col] = df_pred[compound_col].astype(str).str.strip()
    df_merged = pd.merge(df_truth, df_pred[[compound_col, score_col]], on=compound_col)

    if f"{score_col}_x" in df_merged.columns and f"{score_col}_y" in df_merged.columns:
        df_merged[score_col] = df_merged[f"{score_col}_y"]
        df_merged.drop(columns=[f"{score_col}_x", f"{score_col}_y"], inplace=True)
    
    print(f"df_pred:")
    print(df_pred)
    print(f"df_truth:")
    print(df_truth)
    print(f"df_merged:")
    print(df_merged)

    # Define which scores are 'lower is better'
    LOWER_IS_BETTER = {
        "DOCK score", "score", "Affinity Pred Value", "kcal/mol"
    }
    # All others are assumed 'higher is better'
    ascending = True if score_col in LOWER_IS_BETTER else False
    df_sorted = df_merged.sort_values(by=score_col, ascending=ascending).reset_index(drop=True)

    n_total = len(df_sorted)
    n_actives = df_sorted['label'].sum()
    n_top1pct = max(1, int(np.ceil(0.01 * n_total)))
    top1pct_hits = df_sorted.iloc[:n_top1pct]['label'].sum()
    ef1 = (top1pct_hits / n_top1pct) / (n_actives / n_total) if n_actives > 0 else 0.0

    auc_score = roc_auc_score(df_sorted['label'], df_sorted[score_col])
    precision, recall, _ = precision_recall_curve(df_sorted['label'], df_sorted[score_col])
    prc_auc = sk_auc(recall, precision)
    avg_prec = average_precision_score(df_sorted['label'], df_sorted[score_col])
    log_auc = logAUC_from_dataframe(df_sorted, score_col, 'label')  # Compute logAUC

    # Compute BEDROC
    try:
        from scipy.special import expit
        alpha = 20.0
        y_true = df_sorted['label'].values
        n = len(y_true)
        ra = np.where(y_true == 1)[0] + 1
        n_actives = len(ra)
        bedroc = 0.0
        if n_actives > 0:
            sum_exp = np.sum(np.exp(-alpha * (ra / n)))
            bedroc = (sum_exp * alpha / (1 - np.exp(-alpha))) / n_actives
    except ImportError:
        bedroc = None

    # Prepare curves if requested
    curves = {}
    if return_curves:
        fpr, tpr, _ = precision_recall_curve(df_sorted['label'], df_sorted[score_col])
        curves['roc_curve'] = (fpr, tpr)  # Ensure correct format for ROC curve
        curves['pr_curve'] = (recall, precision)  # Ensure correct format for PR curve
        curves['bedroc_curve'] = (np.arange(len(y_true)), np.cumsum(y_true))

    metrics = {
        'EF1%': ef1,
        'AUC': auc_score,
        'AUPRC': prc_auc,
        'Average_Precision': avg_prec,
        'LogAUC': log_auc,  # Include logAUC in metrics
        'BEDROC_alpha_20': bedroc,
        'Score_Used': score_col  # Add score descriptor to metrics dictionary
    }

    if bootstrap:
        bootstrap_results = []
        for _ in tqdm(range(1000), desc=f"Bootstrapping metrics using {score_col}"):
            sampled_df = df_sorted.sample(n=len(df_sorted), replace=True)
            sampled_auc = roc_auc_score(sampled_df['label'], sampled_df[score_col])
            sampled_precision, sampled_recall, _ = precision_recall_curve(sampled_df['label'], sampled_df[score_col])
            sampled_prc_auc = sk_auc(sampled_recall, sampled_precision)
            sampled_avg_prec = average_precision_score(sampled_df['label'], sampled_df[score_col])
            sampled_log_auc = logAUC_from_dataframe(sampled_df, score_col, 'label')  # Compute logAUC for bootstrap
            sampled_top1pct_hits = sampled_df.iloc[:n_top1pct]['label'].sum()
            sampled_ef1 = (sampled_top1pct_hits / n_top1pct) / (n_actives / n_total) if n_actives > 0 else 0.0

            # Compute BEDROC for bootstrap sample
            try:
                sampled_y_true = sampled_df['label'].values
                sampled_ra = np.where(sampled_y_true == 1)[0] + 1
                sampled_n_actives = len(sampled_ra)
                sampled_bedroc = 0.0
                if sampled_n_actives > 0:
                    sampled_sum_exp = np.sum(np.exp(-alpha * (sampled_ra / n)))
                    sampled_bedroc = (sampled_sum_exp * alpha / (1 - np.exp(-alpha))) / sampled_n_actives
            except ImportError:
                sampled_bedroc = None

            bootstrap_results.append({
                'EF1%': sampled_ef1,
                'AUC': sampled_auc,
                'AUPRC': sampled_prc_auc,
                'Average_Precision': sampled_avg_prec,
                'LogAUC': sampled_log_auc,  # Include logAUC in bootstrap results
                'BEDROC_alpha_20': sampled_bedroc
            })

        # Calculate mean and std for each metric
        bootstrap_df = pd.DataFrame(bootstrap_results)
        bootstrap_means = bootstrap_df.mean().to_dict()
        bootstrap_stds = bootstrap_df.std().to_dict()

        # Export bootstrap results to CSV
        bootstrap_csv_path = os.path.join(output_dir, f'bootstrap_results_{sanitized_score_col}.csv')
        bootstrap_df.to_csv(bootstrap_csv_path, index=False)

        # Generate histograms for each metric
        for metric in bootstrap_df.columns:
            mean_val = bootstrap_means[metric]
            std_val = bootstrap_stds[metric]
            plot = (
                ggplot(bootstrap_df, aes(x=metric)) +
                geom_histogram(binwidth=0.1, fill="blue", alpha=0.7) +
                geom_vline(xintercept=mean_val, color="red", linetype="dashed") +
                geom_vline(xintercept=mean_val - std_val, color="black", linetype="dashed") +
                geom_vline(xintercept=mean_val + std_val, color="black", linetype="dashed") +
                annotate("text", x=mean_val, y=0, label=f"Mean: {mean_val:.2f}", color="red", ha="center") +
                theme_minimal() +
                ggtitle(f"Bootstrap Histogram for {metric} using {sanitized_score_col}")
            )
            plot_path = os.path.join(output_dir, f"{metric}_bootstrap_histogram_{sanitized_score_col}.png")
            plot.save(plot_path)

        metrics['Bootstrap_Means'] = bootstrap_means
        metrics['Bootstrap_Stds'] = bootstrap_stds
        metrics['Score_Used'] = score_col  # Add score descriptor to bootstrap metrics

    return (metrics, curves) if return_curves else metrics
