import os
import pandas as pd
import sys
from BoltzCov.analysis.analysis_utils import *
import matplotlib.pyplot as plt
from IPython.display import display # for jupyter notebook 

def analyze_boltz_preds(invitro_file, boltz_outdir, score_col, topN, 
                        protein_name=None, exp_col=None, plot=False):
    """
    Compare Boltz predictions to in vitro IC50 data and compute performance metrics.

    :param invitro_file: Path to experimental CSV file
    :param boltz_outdir: Directory containing Boltz prediction outputs
    :param score_col: Column name of the prediction score
    :param topN: Fraction for top-N enrichment (e.g., 0.1 = 10%)
    :param protein_name: Protein name to construct column 'mean_{protein}_log_ic50 (uM)'
    :param exp_col: Explicit experimental column name (overrides protein_name)
    :param plot: Generate plots if True
    
    :return df_truth_pred: DataFrame with predictions and experimental labels
    :return analysis_dict: Dictionary with metrics, curves, and plots
    """
    
    # Determine experimental column
    if exp_col is None:
        if protein_name is None:
            raise ValueError("Must provide either 'protein_name' or 'exp_col'")
        exp_col = f"mean_{protein_name.lower()}_log_ic50 (uM)"
    
    # Load experimental data
    df_invitro = pd.read_csv(invitro_file)
    
    if exp_col not in df_invitro.columns:
        raise KeyError(f"Column '{exp_col}' not found. Available: {list(df_invitro.columns)}")
    
    # Label actives/inactives
    df_invitro_labeled = process_invitro(invitro_df=df_invitro, exp_col=exp_col)
    
    # Separate actives and inactives
    df_actives = df_invitro_labeled[df_invitro_labeled['is_binder']].copy()
    df_inactives = df_invitro_labeled[~df_invitro_labeled['is_binder']].copy()
    
    print(f"Active compounds: {len(df_actives)}, Inactive compounds: {len(df_inactives)}")
    
    # Assign binary labels
    df_actives['label'] = 1
    df_inactives['label'] = 0
    
    # Combine into ground truth dataset
    df_truth = pd.concat([df_actives, df_inactives], ignore_index=True)
    
    #  dont need this anymore tehcnially 
    boltz_preds_df = read_boltz_predictions(boltz_outdir)
    
    # Merge predictions with ground truth
    df_truth_pred = pd.merge(
        df_truth, 
        boltz_preds_df[["substance_id", score_col]], 
        on="substance_id"
    )
    
    # Generate run name from directory path
    run_name = ': '.join(boltz_outdir.split('/')[-2:])
    
    # Compute performance metrics
    metrics, curves = calculate_metrics(
        topN=topN, 
        df_truth_pred=df_truth_pred, 
        score_col=score_col
    )
    
    analysis_dict = {'metrics_curves': [metrics, curves]}
    
    # Generate plots if requested
    if plot:
        scatter_plot = affinity_scatter(
            df_truth_pred=df_truth_pred, 
            run_name=run_name, 
            score_col=score_col, 
            exp_col=exp_col
        )
        curve_plots = plot_curves(
            run_name=run_name, 
            curves=curves, 
            metrics=metrics
        )
        analysis_dict['plots'] = [scatter_plot, curve_plots]
    
    return df_truth_pred, analysis_dict

def mean_metrics(boltz_reps_outdir, score_col):
    '''
    Computes the mean and standard deviation of a score across replicate predictions.

    :param boltz_reps_outdir: Directory containing Boltz prediction replicates.
    :param score_col: Column name of the score to summarize.
    
    :return: Tuple of (all predictions DataFrame, statistics DataFrame with mean and std per compound).
    '''
    all_pred_reps = read_boltz_predictions(boltz_reps_outdir, reps=True) 
    stats = all_pred_reps.groupby('substance_id')[score_col].agg(['mean', 'std']).reset_index()
    
    return all_pred_reps, stats 


def analyze_mean_preds(invitro_file, stats_df, score_col, exp_col, topN, run_name=None, plot=False):
    df_invitro = pd.read_csv(invitro_file)
    # convert IC50 measured and label actives/inactives
    df_invitro_labeled = process_invitro(invitro_df=df_invitro, exp_col=exp_col) 

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
    df_truth_pred = pd.merge(df_truth, stats_df[["substance_id", score_col]], on="substance_id")
   
    metrics, curves = calculate_metrics(topN=topN, df_truth_pred=df_truth_pred, score_col=score_col)
    analysis_dict = {
        'metrics': metrics,
        'curves': curves
    }

    if plot: 
        exp_col = f'log_{exp_col}'
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

def topN_affinity_scatter(truth_pred_df, analysis_dict, score_col, topN, exp_col, write_output=None, run_name=None):
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

    if score_col == "pred_log10ic50":
        # most negative val needs to be top/best for log(ic50)
        df_sorted_pred = truth_pred_df.sort_values(by=score_col, ascending=True) 
        df_sorted_truth = truth_pred_df.sort_values(by=exp_col, ascending=True)
    else:
        df_sorted_pred = truth_pred_df.sort_values(by=score_col, ascending=False) 
        df_sorted_truth = truth_pred_df.sort_values(by=exp_col, ascending=True) 

    # predicted topN by boltz
    topN_pred = df_sorted_pred.head(int(topN * len(df_sorted_pred)))
    topN_x_pred = topN_pred[['substance_id', score_col]]
    topN_y_pred =  topN_pred[['substance_id', exp_col]]

    # experimental topN
    topN_truth = df_sorted_truth.head(int(topN * len(df_sorted_truth)))
    topN_x_truth= topN_truth[['substance_id', score_col]]
    topN_y_truth=  topN_truth[['substance_id', exp_col]]

    # find common ligands in topN of both predicted and true values
    df_truth_pred = pd.merge(topN_truth[['substance_id', exp_col]], topN_pred[['substance_id', score_col]], on="substance_id")

    if write_output: 
        df_truth_pred2 = pd.merge(
        topN_truth[['substance_id', exp_col]],
        topN_pred[['substance_id', score_col]],
        on="substance_id",
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
    df_truth_pred[exp_col],
    color='purple',
    label=f'top{topN*100}% overlap',
    alpha=0.7
    )

    ax.legend()
    ax.axis('equal')

    # only in ground truth ranked not predicted
    commons = df_truth_pred['substance_id'].tolist()
    topN_x_truth_filtered = topN_x_truth[~topN_x_truth['substance_id'].isin(commons)]
    topN_y_truth_filtered = topN_y_truth[~topN_y_truth['substance_id'].isin(commons)] 

    # plot on existing fig,ax
    ax.scatter(
    topN_x_truth_filtered[score_col],
    topN_y_truth_filtered[exp_col],
    color='red',
    label=f'top{topN*100}% truth',
    alpha=0.7
    )

    # add legend
    ax.legend()

    # only in ground predicted ranked not ground truth
    topN_x_pred_filtered = topN_x_pred[~topN_x_pred['substance_id'].isin(commons)]
    topN_y_pred_filtered = topN_y_pred[~topN_y_pred['substance_id'].isin(commons)]

    ax.scatter(
    topN_x_pred_filtered[score_col],
    topN_y_pred_filtered[exp_col],
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


def plot_combined_curves(systems, curve_type='roc', write_output=None, run_name=None, colors=None):
    """
    Plot combined ROC or PR curves from multiple systems on one figure.

    :param systems: dict of system_name -> {'curves': {...}, 'metrics': {...}}
    :param curve_type: 'roc' or 'pr'
    :param colors: dict mapping system names to colors, or list of colors, or None for default
    """
    if run_name is None:
        run_name = 'all_runs'

    fig, ax = plt.subplots(figsize=(6, 5))

    # Handle colors
    if colors is None:
        # Use default matplotlib color cycle
        color_list = None
    elif isinstance(colors, dict):
        # Colors provided as dictionary mapping system names to colors
        color_list = [colors.get(sys_name, None) for sys_name in systems.keys()]
    elif isinstance(colors, list):
        # Colors provided as list
        color_list = colors
    else:
        color_list = None

    if curve_type == 'roc':
        for idx, (sys_name, data) in enumerate(systems.items()):
            fpr, tpr = data['curves']['auc_roc']
            auc = data['metrics']['ROC AUC']
            color = color_list[idx] if color_list else None
            ax.plot(fpr, tpr, lw=2, color=color, label=f"{sys_name} (AUC={auc:.3f})")

        ax.plot([0, 1], [0, 1], color='gray', lw=1, linestyle='--')
        ax.set_xlabel("False Positive Rate")
        ax.set_ylabel("True Positive Rate")
        plt.title(f"{run_name} ROC Curves")
    
    elif curve_type == 'pr':
        for idx, (sys_name, data) in enumerate(systems.items()):
            recall, precision = data['curves']['pr_auc']
            auc = data['metrics']['PR-AUC']
            no_skill = data['metrics']['No-Skill']
            color = color_list[idx] if color_list else None
            ax.plot([0, 1], [no_skill, no_skill], linestyle='--', color=color)
            ax.plot(recall, precision, lw=2, color=color, label=f"{sys_name} (PR-AUC={auc:.3f})")

        ax.set_xlabel("Recall")
        ax.set_ylabel("Precision")
        plt.title(f"{run_name} Precision-Recall Curves")
    
    elif curve_type == 'logauc':
        # Get LOGAUC_MIN from first system (assuming all use same value)
        first_system = next(iter(systems.values()))
        LOGAUC_MIN = first_system['metrics'].get('LOGAUC_MIN', 0.001)
        LOGAUC_MAX = 1.0
        
        for idx, (sys_name, data) in enumerate(systems.items()):
            df_sorted = data['curves']['logauc'].copy()
            score_col = data['metrics']['Score Used']
            
            # Sort by scores
            if score_col == 'Pred log10(IC50)':
                df_sorted = df_sorted.sort_values(by=score_col, ascending=True)
            else:
                df_sorted = df_sorted.sort_values(by=score_col, ascending=False)
            
            # Compute cumulative TPR and FPR
            total_positives = df_sorted['label'].sum()
            total_negatives = len(df_sorted) - total_positives
            df_sorted['TPR'] = df_sorted['label'].cumsum() / total_positives
            df_sorted['FPR'] = (~df_sorted['label'].astype(bool)).cumsum() / total_negatives
            
            # Filter to range
            mask = (df_sorted['FPR'] >= LOGAUC_MIN) & (df_sorted['FPR'] <= LOGAUC_MAX)
            fpr = df_sorted.loc[mask, 'FPR'].values
            tpr = df_sorted.loc[mask, 'TPR'].values
            
            # Get logAUC value
            logauc = data['metrics']['LogAUC']
            
            # Plot
            color = color_list[idx] if color_list else None
            ax.plot(fpr, tpr, lw=2, color=color, label=f"{sys_name} (logAUC={logauc:.3f})")
        
        # Add random expectation line (only once)
        random_fpr = np.linspace(LOGAUC_MIN, LOGAUC_MAX, 100)
        random_tpr = random_fpr
        ax.plot(random_fpr, random_tpr, linestyle='--', color='gray', lw=1, label='Random')
        
        ax.set_xscale('log')
        ax.set_xlim(LOGAUC_MIN, LOGAUC_MAX)
        ax.set_ylim(0, 1)
        ax.set_xlabel("False Positive Rate (log scale)")
        ax.set_ylabel("True Positive Rate")
        plt.title(f"{run_name} LogAUC Curves")

    ax.legend(bbox_to_anchor=(0.5, 1.05), loc='lower center', ncol=2)
    ax.grid(alpha=0.3)
    # plt.tight_layout()
  
    if write_output:
        file_name = run_name.replace(" ", "_").lower()
        fig.savefig(os.path.join(write_output, f'{file_name}_{curve_type}.png'), dpi=300, bbox_inches="tight")

    return fig, ax