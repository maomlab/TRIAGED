import os
import pandas as pd
import argparse
import json
from analysis_utils import *

def analyze_boltz_preds(invitro_file, boltz_outdir, score_col, plot=False):
    '''
    Compares Boltz predictions to in vitro IC50 data and computes performance metrics.

    :param invitro_file: Path to in vitro IC50 CSV file.
    :param boltz_outdir: Directory containing Boltz prediction outputs.
    :param score_col: Column name of the prediction score by Boltz.

    :return: analysis (dict) with metrics, curves and plots. 
    '''

    df_invitro = pd.read_csv(invitro_file)
    # convert IC50 measured and label actives/inactives
    df_invitro_labeled = process_invitro(invitro_df=df_invitro, score_col="IC50 (nM)") 

    positive_df = df_invitro_labeled[df_invitro_labeled['is_binder'] == True]
    negative_df = df_invitro_labeled[df_invitro_labeled['is_binder'] == False]
    
    print(f"True Positive compounds: {len(positive_df)}, True Negative compounds: {len(negative_df)}")

    df_pos = positive_df.copy()
    df_neg = negative_df.copy()
    df_pos['label'] = 1
    df_neg['label'] = 0 
    
    # final ground truth experimental df
    df_truth = pd.concat([df_pos, df_neg], ignore_index=True)

    # consolidate predictions from boltz as df
    all_predictions = read_boltz_predictions(boltz_outdir)

    # merge the ground truth and predictions on compound_id col 
    # only take score_col from prediction df
    df_truth_pred = pd.merge(df_truth, all_predictions[["compound_id", score_col]], on="compound_id")

    # analysis 
    analysis = {}
    run_name = ': '.join(boltz_outdir.split('/')[-2:])
   
    metrics, curves = calculate_metrics(df_truth_pred=df_truth_pred, score_col=score_col, alpha=20)
    analysis['metrics_curves'] = [metrics, curves]

    if plot: 
        plot = affinity_scatter(df_truth_pred=df_truth_pred, run_name=run_name, score_col=score_col)
        more_figs = plot_curves(run_name=run_name, curves=curves, metrics=metrics)
        analysis['plots'] = [plot, more_figs] # dicts within lists 

    return analysis

def mean_metrics(boltz_outdir, score_col):
    #  fucntion that calcaultes means and deviations in replicates 
    # just returns a dataframe with all the means and deviations for a given score-col?
    all_preds = read_boltz_predictions(boltz_outdir, reps=True) 

    return x

# plot_enrichment(ef_dict) later when you have all reps ; include run name in title of plot, violin plot??? to show spread of all reps?


def main():
    parser = argparse.ArgumentParser(description="Process ground truth and Boltz prediction data. Outputs metrics and plots.")
    parser.add_argument("-i","--input_dir", type=str, required=True, help="Directory containing boltz predictions")
    parser.add_argument("-v","--invitro_data", type=str, required=False, default=None, help="Path to the in vitro data file.")
    parser.add_argument("-p", "--plots", type=bool, action="store_true", default=False, help="Plots curves from analysis and returns figure objects.") 
    parser.add_argument("-m","--mean_metrics", type=bool, action="store_true", default=False, help="Set to True if you have multiple replicate runs and want the means and spread of " )
    args = parser.parse_args()

    