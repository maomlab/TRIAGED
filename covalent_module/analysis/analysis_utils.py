# analysis main script 
# conda env: boltz_analysis_env 
# Adapted from Miguel Limcaoco 
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math 
import os
import numpy as np
import json
from scipy.stats import spearmanr
from scipy.special import expit
from sklearn.linear_model import LinearRegression
from sklearn.metrics import roc_auc_score, roc_curve, precision_recall_curve, auc as sk_auc, average_precision_score

def convert_IC_to_energy(IC):
    """
    Convert log10(IC50) measured in uM to kcal/mol.
    :param IC: IC50.
    :return: converted kcal/mol estimate.
    """
    return (6 - IC) * 1.364

def read_boltz_predictions(predictions_dir, reps=False):
    """
    Reads prediction JSON files from subdirectories and compiles them into a pandas DataFrame.
    :param predictions_dir: Path to the directory containing subdirectories with JSON files.
    :return: Pandas DataFrame with compiled data.
    """
    data = []
    for compound_name in os.listdir(predictions_dir):
        compound_id = compound_name.split('_')[-1]
        compound_dir = os.path.join(predictions_dir, compound_name)
        results = [compound_dir,  f"boltz_results_{compound_name}", "predictions", f"{compound_name}"]
        compound_result =  "/".join(results)
        if not os.path.isdir(compound_result):
            continue

        affinity_file = os.path.join(compound_result, f"affinity_{compound_name}.json")
        confidence_file = os.path.join(compound_result, f"confidence_{compound_name}_model_0.json")

        if os.path.exists(affinity_file) and os.path.exists(confidence_file):
            with open(affinity_file, 'r') as af:
                affinity_data = json.load(af)
                affinity_pred_value = affinity_data.get("affinity_pred_value", None)
                ic50_nm = (10 ** affinity_pred_value) * 1000 
                pred_pic50 = -math.log10((10 ** affinity_pred_value) * 1e-6)
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
                "compound_id": compound_id,
                "Pred log10(IC50)": affinity_pred_value,
                "Pred pIC50": pred_pic50,
                "Pred Label (IC50-like)": True if ic50_nm < 1000 else False,
                "Binding Probability": affinity_probability_binary,
                "Pred Label":  True if affinity_probability_binary > 0.5 else False,
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
            
    if reps:
        all_reps = [os.path.join(predictions_dir, f) for f in os.listdir(predictions_dir)] # only replicate dirs should be in here
        data = []
        for rep in all_reps:
            compounds_in_reps = [d for d in os.listdir(rep) if os.path.isdir(os.path.join(rep, d))]
            rep_name = rep.split('/')[-1]
            for compound_name in compounds_in_reps:
                compound_dir = os.path.join(predictions_dir, rep_name, compound_name)
                results = [compound_dir, f"boltz_results_{compound_name}", "predictions", f"{compound_name}"]
                compound_result =  "/".join(results)
        
                if not os.path.isdir(compound_result):
                    continue

                affinity_file = os.path.join(compound_result, f"affinity_{compound_name}.json")
                confidence_file = os.path.join(compound_result, f"confidence_{compound_name}_model_0.json")

                if os.path.exists(affinity_file) and os.path.exists(confidence_file):
                    with open(affinity_file, 'r') as af:
                        affinity_data = json.load(af)
                        affinity_pred_value = affinity_data.get("affinity_pred_value", None)
                        ic50_nm = (10 ** affinity_pred_value) * 1000
                        pred_pic50 = -math.log10((10 ** affinity_pred_value) * 1e-6)
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

                    # energy_value = convert_IC_to_energy(affinity_pred_value) if affinity_pred_value is not None else None
                    compound_id = compound_name.split('_')[-1]

                    data.append({
                        "rep_id": rep_name,
                        "compound_id": compound_id,
                        "Pred log10(IC50)": affinity_pred_value,
                        "Pred pIC50": pred_pic50,
                        "Pred Label (IC50-like)": True if ic50_nm < 1000 else False,
                        "Binding Probability": affinity_probability_binary,
                        "Pred Label":  True if affinity_probability_binary > 0.5 else False,
                        "Confidence Score": confidence_score,
                        "kcal/mol": None,
                        "PTM": ptm,
                        "IPTM": iptm,
                        "Ligand IPTM": ligand_iptm,
                        "Protein IPTM": protein_iptm,
                        "Complex pLDDT": complex_plddt,
                        "Complex iPLDDT": complex_iplddt,
                        "Complex PDE": complex_pde,
                        "Complex iPDE": complex_ipde})
    return pd.DataFrame(data)

def process_invitro(invitro_df, score_col, threshold=1000):
    '''
    Cleans up in-vitro csv data, converts affinity values from nM to uM, and returns DataFrames with labeled actives and inactives.
    :param invitro_df: Pandas DataFrame
        Experimental in vitro affinity measured for a set of ligands. 
    :param score_col (str): Name of the column with in vitro affinity or IC50 values. 
    :threshold (int): Threshold in nM to label actives vs inactives. Default = 1000nM 
    :return: Pandas DataFrame
        Active and inactive ligands labelled based on the threshold.
    '''
    num_nans = invitro_df[score_col].isna().sum()
    print(f"Number of NaN values in '{score_col}': {num_nans}")
    
    invitro_df["is_binder"] = invitro_df[score_col].apply(lambda x: True if x < int(threshold) else False)

    invitro_df["log10(IC50)"] = invitro_df[score_col].apply(lambda x: math.log10(round(x/1000, 3))) # nM -> uM and log10(uM)
    invitro_df["pIC50"] = invitro_df[score_col].apply(lambda x: -math.log10(x * 1e-9)) # nM -> -log10(M)
    invitro_df.drop(score_col, axis=1, inplace=True)
    # remove nans 
    invitro_df.replace(["nan", "NaN"], np.nan, inplace=True) 

    return pd.DataFrame(invitro_df.dropna())

def enrichment_factor(df_truth_pred, score_col, x=0.10):
    '''
    Computes Enrichment Factor given by 
    [n_hits_x/n_x]/[n_hits_T/n_T]
    [(Total Actives in Top x%)/(Total Compounds in Top x%)]/[(Total Actives)/(Total Compounds)]
    :param score_col: Name of the column to use in ranking ligands.
    :param df_truth_pred: Pandas DataFrame with Actual and Predicted score (eg. IC50, Affinity, etc).
    :param x: Threshold for subsetting Top x%.

    :return: List of Enrichment Factors for all Boltz prediction replicates given. 
    '''
    df_sorted = df_truth_pred.sort_values(by=score_col, ascending=False).reset_index(drop=True) # sort by the pred column, worst on top

    n_T = len(df_sorted) # total compound in dataset
    n_x = max(1, int(np.ceil(x * n_T))) # how many in x% of total ligands; garauntees at least 1 compound
    n_hits_x = df_sorted.iloc[:n_x]['label'].sum() # count actives in subset of x% ligands
    n_hits_T = df_sorted['label'].sum()# total hits in dataset

    ef = (n_hits_x/n_x)/(n_hits_T/n_T)

    return ef

def plot_enrichment(ef_dict):
    '''
    Plots enrichment vs Dataset.
    :param ef_dict (dict): Dictionary with {'Dataset Name':[efs]} as key:value pairs.
    :return: Enrichment vs Dataset plot. 
    '''
    categories = list(ef_dict.keys())
    means = [np.mean(ef_dict[cat]) for cat in categories]
    stds  = [np.std(ef_dict[cat]) for cat in categories]
    x_pos = np.arange(len(categories))

    # Create figure
    fig, ax = plt.subplots(figsize=(6,5))
    ax.errorbar(
        x_pos, means, yerr=stds,
        fmt='o', capsize=5, markersize=8, color='blue', label='Mean EF ± SD'
    )

    ax.set_xticks(x_pos)
    ax.set_xticklabels(categories)
    ax.set_ylabel("Average EF")
    ax.set_title("Average Enrichment Factor with Error Bars")
    ax.grid(axis='y', linestyle='--', alpha=0.7)
    ax.legend()

    return fig, ax
    
def affinity_scatter(df_truth_pred, score_col, exp_col, run_name=None):
    '''
    Plots predicted vs experimental affinity.

    :param df_truth_pred: DataFrame containing predicted and experimental values.
    :param score_col: Column name for predicted affinity scores.
    :param exp_col: Column name for experimental affinity scores.
    :param run_name: Optional name for the run, shown in subtitle.

    :return: dict with {fig: ax} of the scatter plot.
    '''

    y_label = exp_col  # experimental
    x_label = score_col  # predicted

    x = df_truth_pred[x_label].values.reshape(-1, 1)
    y = df_truth_pred[y_label].values

    # R squared
    model = LinearRegression()
    model.fit(x, y)
    y_pred = model.predict(x)
    r_squared = model.score(x, y)
    print(f"R²: {r_squared:.2f}")

    # Spearman correlation
    spearman_corr, spearman_p = spearmanr(df_truth_pred[x_label], df_truth_pred[y_label])
    print(f"Spearman r: {spearman_corr:.2f}, p-value: {spearman_p:.3g}")

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(x, y, alpha=0.7, color='grey')
    ax.plot(x, y_pred, color='grey', linestyle=":", label=f"Linear fit (R² = {r_squared:.2f})")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    # Titles
    plt.suptitle(f"Experimental vs \n{score_col} Correlation", fontsize=12, fontweight='bold', y=1.00)
    plt.title(f"for {run_name}", fontsize=10, color='black', y=1.00)

    ax.legend()
    ax.text(
    0.95, 0.90, f"Spearman r = {spearman_corr:.2f}", transform=ax.transAxes, ha='right', va='top')
    # ax.axis('equal')
    ax.grid(True)

    return {fig:ax}

def bedroc_calc(df_truth_pred, score_col, alpha=20):
    
    df_sorted = df_truth_pred.sort_values(by=score_col, ascending=False).reset_index(drop=True)

    y_true = df_sorted['label'].values
    y_scores = df_sorted[score_col].values

    n = len(y_true)
    ra = np.where(y_true == 1)[0] + 1  # ranks of actives
    n_actives = len(ra)
    bedroc = np.nan
    if n_actives > 0:
        sum_exp = np.sum(np.exp(-alpha * ra / n))
        rand_exp = n_actives / n * (1 - np.exp(-alpha))
        bedroc = (sum_exp * alpha / n_actives) / (1 - np.exp(-alpha)) - rand_exp / (1 - np.exp(-alpha))

    bedroc_curve = (np.arange(len(y_true)), np.cumsum(y_true))
    
    return bedroc, bedroc_curve

def calculate_logAUC(df_truth_pred, score_col):
    '''
    Computes the adjusted logAUC for a given score column.
    
    :param df_truth_pred: DataFrame, true labels and predicted scores.
    :param score_col: str, column name of predicted scores.

    :return: float, logAUC value adjusted for random expectation.
    '''
    
    LOGAUC_MAX = 1.0   # this should not change
    LOGAUC_MIN = 0.001 # this can be adjusted for large datasets with strong early enrichment
    RANDOM_LOGAUC = (LOGAUC_MAX - LOGAUC_MIN) / np.log(10) / np.log10(LOGAUC_MAX / LOGAUC_MIN)

    # Sort the DataFrame by scores in descending order
    df_sorted = df_truth_pred.sort_values(by=score_col, ascending=False)

    # Compute cumulative true positive rate (TPR) and false positive rate (FPR)
    total_positives = df_truth_pred['label'].sum()
    total_negatives = len(df_truth_pred) - total_positives

    df_sorted['TPR'] = df_sorted['label'].cumsum() / total_positives
    df_sorted['FPR'] = (~df_sorted['label'].astype(bool)).cumsum() / total_negatives

    # Generate points array (FPR, TPR)
    points = df_sorted[['FPR', 'TPR']].values

    # Filter and normalize points
    npoints = []
    for x in points:
        if (x[0] >= LOGAUC_MIN) and (x[0] <= LOGAUC_MAX):
            npoints.append([x[0], x[1]])

    # Compute log AUC
    area = 0.0
    for point2, point1 in zip(npoints[1:], npoints[:-1]):
        if point2[0] - point1[0] < 0.000001:
            continue

        dx = point2[0] - point1[0]
        dy = point2[1] - point1[1]
        intercept = point2[1] - (dy) / (dx) * point2[0]
        area += dy / np.log(10) + intercept * (np.log10(point2[0]) - np.log10(point1[0]))

    # Return adjusted log AUC
    return area / np.log10(LOGAUC_MAX / LOGAUC_MIN) - RANDOM_LOGAUC

def calculate_metrics(df_truth_pred, score_col, alpha=20):
    '''
    Computes performance metrics and curve data for a given score column.
    
    :param df_truth_pred: DataFrame, true labels and predicted scores.
    :param score_col: str, column name of predicted scores.
    :param alpha: float, BEDROC alpha parameter (default=20).

    :return: tuple(dict, dict), metrics values and curve data.
    '''

    df_sorted = df_truth_pred.sort_values(by=score_col, ascending=False).reset_index(drop=True)

    y_true = df_sorted['label'].values
    y_scores = df_sorted[score_col].values
    # EF
    ef = enrichment_factor(df_truth_pred, score_col, x=0.10)

    # ROC
    if score_col == 'Pred log10(IC50)':
        fpr, tpr, _ = roc_curve(y_true, -y_scores)
        roc_auc = roc_auc_score(y_true, -y_scores)
    else:
        fpr, tpr, _ = roc_curve(y_true, y_scores)
        roc_auc = roc_auc_score(y_true, y_scores)

    # PR curve
    precision, recall, _ = precision_recall_curve(y_true, y_scores)
    pr_auc = sk_auc(recall, precision)
    avg_prec = average_precision_score(y_true, y_scores)

    #BEROC
    bedroc, bedroc_curve = bedroc_calc(df_sorted, score_col, alpha)

    # LogAUC
    log_auc = calculate_logAUC(df_truth_pred, score_col)

    curves = {
        'auc_roc': (fpr, tpr),
        'pr_auc': (recall, precision)
    }
    curves['bedroc'] = bedroc_curve

    metrics = {
        'EF%': round(ef,3),
        'ROC AUC': round(roc_auc,3),
        'PR-AUC': round(pr_auc,3),
        'Average Precision': round(avg_prec,3),
        'LogAUC': round(log_auc,3),
        f'BEDROC': round(bedroc,3),
        'Score Used': score_col
    }

    return metrics, curves

def plot_curves(curves={}, metrics={}, run_name=None):
    '''
    Plots ROC, Precision-Recall, and BEDROC (cumulative hits) curves.

    :param curves: dict with keys ['auc_roc', 'pr_auc', 'bedroc'].
    :param metrics: dict with keys ['ROC AUC', 'PR-AUC', 'Average Precision', 'BEDROC', 'Score Used'].
    
    :return: dict, figure and axis objects for each curve.
    '''
    figs = {}
    for metric in metrics: 
      
        if metric == 'ROC AUC':
            # --- ROC Curve ---
            fig_roc, ax_roc = plt.subplots(figsize=(6, 5))
            fpr, tpr = curves['auc_roc']
            ax_roc.plot(fpr, tpr, color='darkorange', lw=2, label=f"ROC AUC = {metrics['ROC AUC']:.3f}")
            ax_roc.plot([0, 1], [0, 1], color='gray', lw=1, linestyle='--')
            ax_roc.set_xlabel("False Positive Rate")
            ax_roc.set_ylabel("True Positive Rate")
            plt.suptitle(f"ROC Curve using {metrics['Score Used']}", fontsize=12, fontweight='bold', y=1.00)
            plt.title(f"for {run_name}", fontsize=12, color='black', y=1.02)
            ax_roc.legend()
            ax_roc.grid(alpha=0.3)
            figs['roc_curve'] = [fig_roc, ax_roc]
        
        elif metric == 'PR-AUC':
            # --- Precision-Recall Curve ---
            fig_pr, ax_pr = plt.subplots(figsize=(6, 5))
            recall, precision = curves['pr_auc']
            ax_pr.plot(recall, precision, color='blue', lw=2,
                    label=f"PR AUC = {metrics['PR-AUC']:.3f}\nAP = {metrics['Average Precision']:.3f}")
            ax_pr.set_xlabel("Recall")
            ax_pr.set_ylabel("Precision")
            plt.suptitle(f"Precision-Recall Curve\nusing {metrics['Score Used']}", fontsize=12, fontweight='bold', y=1.00)
            plt.title(f"for {run_name}", fontsize=10, color='black', y=0.99)
            ax_pr.legend()
            ax_pr.grid(alpha=0.3)
            figs['pr_curve'] = [fig_pr, ax_pr]

        elif metric == 'BEDROC':
            # --- BEDROC (Cumulative Hits) Curve ---
            fig_bedroc, ax_bedroc = plt.subplots(figsize=(6, 5))
            x, cum_hits = curves['bedroc']
            ax_bedroc.plot(x, cum_hits, color='green', lw=2, label=f"BEDROC = {metrics['BEDROC']:.3f}")
            ax_bedroc.set_xlabel("Compound Rank")
            ax_bedroc.set_ylabel("Cumulative Hits")
            plt.suptitle(f"BEDROC-like Cumulative Hit Curve\nusing {metrics['Score Used']}", fontsize=12, fontweight='bold', y=1.00)
            plt.title(f"for {run_name}", fontsize=10, color='black', y=0.99)
            ax_bedroc.legend()
            ax_bedroc.grid(alpha=0.3)
            figs['bedroc_curve'] = [fig_bedroc, ax_bedroc]

    return figs