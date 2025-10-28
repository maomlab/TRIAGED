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

def convert_IC_to_energy(IC):
    """
    Convert log10(IC50) measured in uM to kcal/mol.
    :param IC: IC50.
    :return: converted kcal/mol estimate.
    """
    return (6 - IC) * 1.364

def read_boltz_predictions(predictions_dir):
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
                "Affinity Pred log10(IC50)": affinity_pred_value,
                "IC50-like (nM)": ic50_nm,
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

    return pd.DataFrame(data)

def process_invitro(invitro_df, affinity_col, threshold=1000):
    '''
    Cleans up in vitro csv data, converts affinity values from nM to uM, and returns DataFrames with actives and inactives.
    :param invitro_df: Pandas DataFrame
        Experimental in vitro affinity measured for a set of ligands. 
    :param affinity_col (str): Name of the column with in vitro affinity or IC50 values. 
    :threshold (int): Threshold in nM to label actives vs inactives. Default = 1000nM 
    :return: Pandas DataFrame
        Active and inactive ligands labelled based on the threshold.
    '''
    num_nans = invitro_df[affinity_col].isna().sum()
    print(f"Number of NaN values in '{affinity_col}': {num_nans}")
    
    invitro_df["is_binder"] = invitro_df[affinity_col].apply(lambda x: True if x < int(threshold) else False)
    
    positive_df = invitro_df[invitro_df['is_binder'] == True]
    negative_df = invitro_df[invitro_df['is_binder'] == False]
    print(f"True Positive compounds: {len(positive_df)}, True Negative compounds: {len(negative_df)}")

    invitro_df["Affinity log10(IC50)"] = invitro_df[affinity_col].apply(lambda x: math.log10(round(x/1000, 3))) # nM to uM and log10(uM)
    invitro_df.drop(affinity_col, axis=1, inplace=True)
    # remove nans 
    invitro_df.replace(["nan", "NaN"], np.nan, inplace=True) 

    return pd.DataFrame(invitro_df.dropna())

def enrichment_factor(invitro_df, score_col, compound_col="compound_id", pred_dfs=[], x=0.10):
    '''
    Computes Enrichment Factor given by 
    [n_hits_x/n_x]/[n_hits_T/n_T]
    [(Total Actives in Top x%)/(Total Compounds in Top x%)]/[(Total Actives)/(Total Compounds)]
    :param invitro_df: Pandas DataFrame.
    :param score_col: Name of the column to use in ranking ligands.
    :param compound_col: Name of the column with compound ids.
    :param pred_dfs: List of Pandas DataFrames.
    :param x: Threshold for subsetting Top x%.

    :return: List of Enrichment Factors for all Boltz prediction replicates given. 
    '''
    positive_df = invitro_df[invitro_df['is_binder'] == True]
    negative_df = invitro_df[invitro_df['is_binder'] == False]

    df_pos = positive_df.copy()
    df_neg = negative_df.copy()
    df_pos['label'] = 1
    df_neg['label'] = 0 

    df_truth = pd.concat([df_pos, df_neg], ignore_index=True)
    efs = []
    for df_pred in pred_dfs:
        df_truth_pred = pd.merge(df_truth, df_pred[[compound_col, score_col]], on=compound_col)
        df_sorted = df_truth_pred.sort_values(by=score_col, ascending=False).reset_index(drop=True) # sort by the pred column, worst on top
   
        n_T = len(df_sorted) # total compound in dataset
        n_x = max(1, int(np.ceil(x * n_T))) # how many in x% of total ligands; garauntees at least 1 compound
        n_hits_x = df_sorted.iloc[:n_x]['label'].sum() # count actives in subset of x% ligands
        n_hits_T = df_sorted['label'].sum()# total hits in dataset

        efs.append((n_hits_x/n_x)/(n_hits_T/n_T))

    return efs

def plot_enrichment(ef_dict, save=False):
    '''
    Plots enrichment vs Dataset.
    :param ef_dict (dict): Dictionary with {'Dataset Name':[efs]} as key:value pairs.
    :return: Enrichment vs Dataset plot. 
    '''
    categories = list(ef_dict.keys())
    means = [np.mean(ef_dict[cat]) for cat in categories]
    stds  = [np.std(ef_dict[cat]) for cat in categories]
    x_pos = np.arange(len(categories))

    # Plot points with error bars
    plt.figure(figsize=(6,5))
    plt.errorbar(x_pos, means, yerr=stds, fmt='o', capsize=5, markersize=8, color='blue')
    plt.xticks(x_pos, categories)
    plt.ylabel("Average EF")
    plt.title("Average Enrichment Factor with Error Bars")
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    if save:
        # Save plot to file
        plt.savefig("average_ef_plot.png", dpi=300, bbox_inches='tight')  
        plt.close()
        
    else:
        plt.show()

