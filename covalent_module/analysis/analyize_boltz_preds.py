import os
import pandas as pd
import argparse
import json
from analysis_utils import *
import matplotlib.pyplot as plt
from IPython.display import display # for jupyter notebook 

def analyze_boltz_preds(invitro_file, boltz_outdir, score_col, exp_col, plot=False):
    '''
    Compares Boltz predictions to in vitro IC50 data and computes performance metrics.

    :param invitro_file: Path to in vitro IC50 CSV file.
    :param boltz_outdir: Directory containing Boltz prediction outputs.
    :param score_col: Column name of the prediction score by Boltz.

    :return all_predictions: Pandas DataFrame of Boltz2 predictions compiled.
    :return analysis: (dict) with metrics, curves and plots. 
    '''

    df_invitro = pd.read_csv(invitro_file)
    # convert IC50 measured and label actives/inactives
    df_invitro_labeled = process_invitro(invitro_df=df_invitro, score_col=exp_col) 

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
    boltz_preds_df = read_boltz_predictions(boltz_outdir)

    # merge the ground truth and predictions on compound_id col 
    # only take score_col from prediction df
    df_truth_pred = pd.merge(df_truth, boltz_preds_df[["compound_id", score_col]], on="compound_id")

    # analysis 
    analysis_dict = {}
    run_name = ': '.join(boltz_outdir.split('/')[-2:])
   
    metrics, curves = calculate_metrics(df_truth_pred=df_truth_pred, score_col=score_col, alpha=20)
    analysis_dict['metrics_curves'] = [metrics, curves]

    if plot: 
        plot = affinity_scatter(df_truth_pred=df_truth_pred, run_name=run_name, score_col=score_col)
        more_figs = plot_curves(run_name=run_name, curves=curves, metrics=metrics)
        analysis_dict['plots'] = [plot, more_figs] # dicts within lists 

    return df_truth_pred, analysis_dict

def mean_metrics(boltz_reps_outdir, score_col):
    '''
    Computes the mean and standard deviation of a score across replicate predictions.

    :param boltz_reps_outdir: Directory containing Boltz prediction replicates.
    :param score_col: Column name of the score to summarize.
    
    :return: Tuple of (all predictions DataFrame, statistics DataFrame with mean and std per compound).
    '''
    all_pred_reps = read_boltz_predictions(boltz_reps_outdir, reps=True) 
    stats = all_pred_reps.groupby('compound_id')[score_col].agg(['mean', 'std']).reset_index()
    
    return all_pred_reps, stats 


def analyze_mean_preds(invitro_file, stats_df, score_col, exp_col, run_name=None, plot=False):
    
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

    stats_df.rename(columns={'mean': score_col}, inplace=True)
    df_truth_pred = pd.merge(df_truth, stats_df[["compound_id", score_col]], on="compound_id")
   
    metrics, curves = calculate_metrics(df_truth_pred=df_truth_pred, score_col=score_col, alpha=20)
    analysis_dict = {
        'metrics': metrics,
        'curves': curves
    }

    if plot: 
        scatter = affinity_scatter(df_truth_pred=df_truth_pred, run_name=run_name, score_col=score_col, exp_col=exp_col)
        curve_figs = plot_curves(run_name=run_name, curves=curves, metrics=metrics)
        analysis_dict['plots'] = {
            'scatter': scatter,
            'curves': curve_figs
        }

    return df_truth_pred, analysis_dict


def view_plot(fig, save_path=None, show=True, close=False):
    '''
    Show or save a matplotlib figure.
    :param fig: matplotlib.figure.Figure
    When using analyize_boltz_preds(), analysis_dict is an output that has plots. 
    
    Usage eg. 
    boltz_preds_df, analysis_dict = analyize_boltz_preds(...)
    plots = analysis_dict['plots'][0]
    fig, ax = list(plots.items())[0] 
    show_or_save_plot(fig)

    another eg.
    analysis_dict = analyze_mean_preds(...)
    fig, ax = analysis_dict['plots']['curves']['roc_curve']
    show_or_save_plot(fig)
    '''
    if save_path:
        fig.savefig(save_path, bbox_inches='tight', dpi=300)
        print(f"✅ Saved plot to {save_path}")

    if show:
        display(fig) 
    if close:
        plt.close(fig)