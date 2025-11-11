import os
import pandas as pd
import argparse
import json
import sys
from analysis_utils import *
import matplotlib.pyplot as plt
from IPython.display import display # for jupyter notebook 

def analyze_boltz_preds(invitro_file, boltz_outdir, score_col, exp_col="IC50 (nM)", plot=False):
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
   
    metrics, curves = calculate_metrics(df_truth_pred=df_truth_pred, score_col=score_col)
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
   
    metrics, curves = calculate_metrics(df_truth_pred=df_truth_pred, score_col=score_col)
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


def view_plot(fig, save_path=None, show=True, close=False, run_name=None):
    '''
    Show or save a matplotlib figure.
    :param fig: matplotlib.figure.Figure
    When using analyize_boltz_preds(), analysis_dict is an output that has plots. 
    
    Usage eg. 
    boltz_preds_df, analysis_dict = analyize_boltz_preds(...)
    plots = analysis_dict['plots'][0]
    fig, ax = list(plots.items())[0] 
    view_plot(fig)

    another eg.
    analysis_dict = analyze_mean_preds(...)
    fig, ax = analysis_dict['plots']['curves']['roc_curve']
    view_plot(fig)
    '''
    if run_name is None:
        run_name = 'plot.png'

    if save_path:
        fig.savefig(os.path.join(save_path, run_name), bbox_inches='tight', dpi=300)
        print(f"✅ Saved plot to {save_path}")

    if show:
        display(fig) 
    if close:
        plt.close(fig)

def topN_affinity_scatter(truth_pred_df, analysis_dict, score_col, topN=0.1, write_output=None, run_name=None):
    '''
    Plots topN ligands predicted by boltz vs topN experimentally ranked. 
    Uses truth_pred_df, analysis_dict output from analyze_mean_preds only!
    :param topN (int): Must be less than 1 (topN=topN%/100)
    :param write_output (str): If path given, the merged dataframe showing which ligands of topN come 
    from predicted, experimental, or both will be written to given path. Figure too will output here.
    '''
    if topN > 1:
        print('topN must be less than 1 (topN=topN%/100)')
        sys.exit(1)

    if run_name is None:
        run_name = 'merged_ligs'

    if score_col == "Pred log10(IC50)":
        # most negative val needs to be top/best for log(ic50)
        df_sorted_pred = truth_pred_df.sort_values(by=score_col, ascending=True) 
        df_sorted_truth = truth_pred_df.sort_values(by="log10(IC50)", ascending=True)
    else:
        df_sorted_pred = truth_pred_df.sort_values(by=score_col, ascending=False) 
        df_sorted_truth = truth_pred_df.sort_values(by="log10(IC50)", ascending=True) # always using log10(IC50) from experiments

    # predicted topN by boltz
    topN_pred = df_sorted_pred.head(int(topN * len(df_sorted_pred)))
    topN_x_pred = topN_pred[['compound_id', score_col]]
    topN_y_pred =  topN_pred[['compound_id', "log10(IC50)"]]

    # experimental topN
    topN_truth = df_sorted_truth.head(int(topN * len(df_sorted_truth)))
    topN_x_truth= topN_truth[['compound_id', score_col]]
    topN_y_truth=  topN_truth[['compound_id', "log10(IC50)"]]

    # find common ligands in topN of both predicted and true values
    df_truth_pred = pd.merge(topN_truth[['compound_id', "log10(IC50)"]], topN_pred[['compound_id', score_col]], on="compound_id")

    if write_output: 
        df_truth_pred2 = pd.merge(
        topN_truth[['compound_id', "log10(IC50)"]],
        topN_pred[['compound_id', score_col]],
        on="compound_id",
        how="outer",       # use outer to include all from both
        indicator=True     # adds a column "_merge"
    )

        os.makedirs(write_output, exist_ok=True)
        df_truth_pred2.to_csv(os.path.join(write_output, f'{run_name}_source.csv'), index=False)

    fig = list(analysis_dict['plots']['scatter'].keys())[0] 
    ax  = list(analysis_dict['plots']['scatter'].values())[0]   

    # overlap compounds in topN
    ax.scatter(
    df_truth_pred[score_col],
    df_truth_pred['log10(IC50)'],
    color='purple',
    label=f'top{topN*100}% overlap',
    alpha=0.7
    )

    ax.legend()
    ax.axis('equal')

    # only in ground truth ranked not predicted
    commons = df_truth_pred['compound_id'].tolist()
    topN_x_truth_filtered = topN_x_truth[~topN_x_truth['compound_id'].isin(commons)]
    topN_y_truth_filtered = topN_y_truth[~topN_y_truth['compound_id'].isin(commons)] 

    # plot on existing fig,ax
    ax.scatter(
    topN_x_truth_filtered[score_col],
    topN_y_truth_filtered['log10(IC50)'],
    color='red',
    label=f'top{topN*100}% truth',
    alpha=0.7
    )

    # add legend
    ax.legend()

    # only in ground predicted ranked not ground truth
    topN_x_pred_filtered = topN_x_pred[~topN_x_pred['compound_id'].isin(commons)]
    topN_y_pred_filtered = topN_y_pred[~topN_y_pred['compound_id'].isin(commons)]

    ax.scatter(
    topN_x_pred_filtered[score_col],
    topN_y_pred_filtered['log10(IC50)'],
    color='blue',
    label=f'top{topN*100}% pred',
    alpha=0.7
    )

    # add legend
    ax.legend()
    if write_output:
        fig.savefig(os.path.join(write_output, f'{run_name}_top{topN*100}_scatter.png'), dpi=300, bbox_inches="tight")

    if write_output is None:
        display(fig)


def plot_combined_curves(systems, curve_type='roc', write_output=None, run_name=None):
    """
    Plot combined ROC or PR curves from multiple systems on one figure.

    :param systems: dict of system_name -> {'curves': {...}, 'metrics': {...}}
    :param curve_type: 'roc' or 'pr'
    """
    if run_name is None:
        run_name = 'all_runs'

    fig, ax = plt.subplots(figsize=(6, 5))

    if curve_type == 'roc':
        for sys_name, data in systems.items():
            fpr, tpr = data['curves']['auc_roc']
            auc = data['metrics']['ROC AUC']
            ax.plot(fpr, tpr, lw=2, label=f"{sys_name} (AUC={auc:.3f})")

        ax.plot([0, 1], [0, 1], color='gray', lw=1, linestyle='--')
        ax.set_xlabel("False Positive Rate")
        ax.set_ylabel("True Positive Rate")
        plt.title(f"{run_name} ROC Curves")
    
    elif curve_type == 'pr':
        for sys_name, data in systems.items():
            recall, precision = data['curves']['pr_auc']
            auc = data['metrics']['PR-AUC']
            ax.plot(recall, precision, lw=2, label=f"{sys_name} (PR-AUC={auc:.3f})")

        ax.set_xlabel("Recall")
        ax.set_ylabel("Precision")
        plt.title(f"{run_name} Precision-Recall Curves")

    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()
  
    if write_output:
        file_name = run_name.replace(" ", "_").lower()
        fig.savefig(os.path.join(write_output, f'{file_name}_roc.png'), dpi=300, bbox_inches="tight")

    if write_output is None:
        display(fig)