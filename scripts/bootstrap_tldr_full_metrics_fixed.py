#!/usr/bin/python
from __future__ import print_function
import collections
from tqdm import tqdm
""" 
Bootstrap ranked list from docking
To obtain statistical significance between AUC and/or LogAUC
"""

#==== 11-24-2020 YY
#==== add flag for bootstrapping vs. standard logAUC plot

#==== 3-25-2020 YY
#==== change plot & Z-test


import os, sys, glob
import argparse 

import matplotlib  # must import first
import matplotlib.pyplot as plt
matplotlib.use('Agg')
import seaborn as sns # Making graphs 

import numpy   as np
from scipy import stats
from statsmodels.stats import weightstats as stests
import pandas  as pd  # For data frame


from sklearn.metrics import roc_auc_score, average_precision_score, precision_recall_curve, auc as sk_auc

def compute_vscreen_metrics(df_pos, df_neg, df_pred, score_col='score', compound_col='compound_id', return_curves=False):
    df_pos = df_pos.copy()
    df_neg = df_neg.copy()
    df_pos['label'] = 1
    df_neg['label'] = 0
    df_truth = pd.concat([df_pos, df_neg], ignore_index=True)
    df_merged = pd.merge(df_truth, df_pred[[compound_col, score_col]], on=compound_col)
    df_sorted = df_merged.sort_values(by=score_col, ascending=False).reset_index(drop=True)

    n_total = len(df_sorted)
    n_actives = df_sorted['label'].sum()
    n_top1pct = max(1, int(np.ceil(0.01 * n_total)))
    top1pct_hits = df_sorted.iloc[:n_top1pct]['label'].sum()
    ef1 = (top1pct_hits / n_top1pct) / (n_actives / n_total) if n_actives > 0 else 0.0

    auc_score = roc_auc_score(df_sorted['label'], df_sorted[score_col])
    precision, recall, _ = precision_recall_curve(df_sorted['label'], df_sorted[score_col])
    prc_auc = sk_auc(recall, precision)
    avg_prec = average_precision_score(df_sorted['label'], df_sorted[score_col])

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
        'BEDROC_alpha_20': bedroc
    }

    return (metrics, curves) if return_curves else metrics





my_parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# my_parser.add_argument('one_sys', 
        # help='absolute path to the folder.')
my_parser.add_argument('-l', '-lig', dest='lig_file', required=True,
        help='File contain ligand names.', default='ligands.name')
my_parser.add_argument('-d', '-dec', dest='dec_file', required=True,
        help='File contain decoys names.', default='decoys.name')

my_parser.add_argument('-s1', dest='ref', required=True, 
        # help='Folder with method 1 (ref) score file: extract_all.sort.uniq.txt')
        help='File with ranked list of molecules in the format of (MOLECULE_NAME DOCKING_SCORE)')
my_parser.add_argument('-s2', dest='new', 
        # help='Folder(s) with new method(s) score. extract_all.sort.uniq.txt')
        help='File with ranked list for new method to compare to. In the same format')
my_parser.add_argument('-m', '-metrics', dest='metrics', choices=['AUC', 'logAUC', 'EF1%', 'BEDROC', 'AUPRC', 'avg_precision', 'all'],
        help='choose to use AUC or logAUC as the metrics.', default='both')

# bootstrapping section
my_parser.add_argument('-b', '-bootstrap',  dest='bootstrap',  
        action = 'store_true',
        help='Flag for bootstrap. ')
my_parser.add_argument('-n', '-num',  dest='num', type=int, 
        choices=range(10,1001), metavar="{10..1000}",
        help='Number of bootstrap replicate. ', default=50)
my_parser.add_argument('-t', '-test', dest='test', choices=['ttest_ind', 'ttest_rel', 'ztest'],
        help='choose stats test to report p value. Default=', default='ztest')
my_parser.add_argument('-debug', dest='debug', action='store_true',
        help='Enable debug verbose mode.')

def debug_print(message):
    if debug_mode:
        print(f"[DEBUG] {message}")

def parse_input():
        #=================#
        # Parse arguments #
        #=================#
        opt = my_parser.parse_args()

        global logging_start, logging_end, debug_mode
        logging_start = "\n" + "#"*50 + "\n" 
        logging_end   = "\n" + "#"*50 + "\n"
        debug_mode = opt.debug

        print(logging_start)
        debug_print("Parsing input arguments...")

        print(f'Ligands file:                     {opt.lig_file}')
        print(f'Decoys file:                      {opt.dec_file}')
        debug_print(f"Ligands file: {opt.lig_file}, Decoys file: {opt.dec_file}")

        if not os.path.isfile(opt.lig_file) or not os.path.isfile(opt.dec_file):
                print("Ligand/decoys files not exits... Check...")
                print(f'Ligands file:                     {opt.lig_file}')
                print(f'Decoys file:                      {opt.dec_file}')
                sys.exit()

        # compare two methods/scoring functions
        if opt.new:
                print(f'Reference:                        {opt.ref}')
                print(f'New method:                       {opt.new}')
                print(f'Number of bootstrap replicate:    {opt.num}')
                print(f'Test for p-value:                 {opt.test}')
                debug_print(f"Reference file: {opt.ref}, New method file: {opt.new}")
                debug_print(f"Bootstrap replicates: {opt.num}, Test: {opt.test}")
                if not os.path.isfile(opt.new) or not os.path.isfile(opt.ref):
                        print("File with ranked list not exits... Check...")
                        print(f'Ranked list from reference method: {opt.ref}')
                        print(f'Ranked list from new method:       {opt.new}')
                        sys.exit()

        else:
                print(f'Rank list:                        {opt.ref}')
                if opt.bootstrap:
                        opt.num = 1000
                        print(f'Number of bootstrap replicate:    {opt.num}')
                        debug_print(f"Bootstrap replicates set to default: {opt.num}")

        print(f'Metrics:                          {opt.metrics}')
        debug_print(f"Metrics selected: {opt.metrics}")

        global num_bootstrap_replicate
        num_bootstrap_replicate = opt.num

        print(logging_end)

        return opt

def export_curves(df_pos, df_neg, df_scores, output_prefix="whole_dataset"):
    """Export ROC, PR, and BEDROC curves for the whole dataset."""
    print(f"Exporting curves for the whole dataset with prefix '{output_prefix}'...")
    metrics, curves = compute_vscreen_metrics(df_pos, df_neg, df_scores, return_curves=True)
    print(curves)
    if 'roc_curve' in curves:
        fpr, tpr = curves['roc_curve']  # Unpack correctly formatted ROC curve
        plt.figure()
        plt.plot(fpr, tpr, label="ROC Curve")
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate")
        plt.title("ROC Curve")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"{output_prefix}_roc_curve.png")
        plt.close()
        print(f"Saved ROC curve as '{output_prefix}_roc_curve.png'.")

    if 'pr_curve' in curves:
        recall, precision = curves['pr_curve']  # Unpack correctly formatted PR curve
        plt.figure()
        plt.plot(recall, precision, label="PR Curve")
        plt.xlabel("Recall")
        plt.ylabel("Precision")
        plt.title("Precision-Recall Curve")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"{output_prefix}_pr_curve.png")
        plt.close()
        print(f"Saved PR curve as '{output_prefix}_pr_curve.png'.")

    if 'bedroc_curve' in curves:
        x, y = curves['bedroc_curve']
        plt.figure()
        plt.plot(x, y, label="BEDROC Curve")
        plt.xlabel("X-axis")
        plt.ylabel("Y-axis")
        plt.title("BEDROC Curve")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"{output_prefix}_bedroc_curve.png")
        plt.close()
        print(f"Saved BEDROC curve as '{output_prefix}_bedroc_curve.png'.")

def main():

        opt = parse_input()
        debug_print("Starting main function...")

        #======================#
        # get ligands and decoys names from file  #
        #======================#
        debug_print("Reading ligands and decoys...")
        lig_list, decoy_list = read_ligands_decoys(opt.lig_file, opt.dec_file)

        # read reference method scores
        debug_print("Reading reference scores...")
        ref_scores = read_score(opt.ref)

        # Prepare data for metrics computation
        debug_print("Preparing data for metrics computation")

        df_pos = pd.DataFrame({'compound_id': lig_list})
        df_neg = pd.DataFrame({'compound_id': decoy_list})
        df_scores = pd.DataFrame([[k, float(v[-1])] for k, v in ref_scores.items()],
                                 columns=['compound_id', 'score'])
        print(df_scores.head())
        print(df_pos.head())
        print(df_neg.head())
        # Export curves for the whole dataset
        debug_print("Export curves for the whole dataset")

        export_curves(df_pos, df_neg, df_scores, output_prefix="whole_dataset")

        # ----------------------------- #
        # just logAUC
        # ----------------------------- #
        if not opt.bootstrap:
                debug_print("Processing without bootstrapping...")
                if opt.metrics == "AUC":
                        debug_print("Calling enrichment_standard with AUC...")
                        enrichment_standard(ref_scores, lig_list, decoy_list, "AUC")
                elif opt.metrics == "logAUC":
                        debug_print("Calling enrichment_standard with logAUC...")
                        enrichment_standard(ref_scores, lig_list, decoy_list, "logAUC")
                elif opt.metrics == "both":
                        debug_print("Calculating both AUC and logAUC...")
                        mols = ref_scores.keys()
                        ranked_list = [ref_scores[key] for key in mols]
                        ranked_list.sort(key=lambda x: float(x[-1]))
                        points = do_roc(ranked_list, lig_list, decoy_list)
                        points = interpolate_curve(points)
                        auc    = AUC(points)*100
                        logauc = logAUC(points)*100
                        plot_ROC(points, auc, logauc)
                        sys.exit(0)
                elif opt.metrics in ["EF1%", "BEDROC", "AUPRC", "avg_precision", "all"]:
                    debug_print("Computing additional metrics...")
                    df_pos = pd.DataFrame({'compound_id': lig_list})
                    df_neg = pd.DataFrame({'compound_id': decoy_list})
                    df_scores = pd.DataFrame([[k, float(v[-1])] for k, v in ref_scores.items()],
                                             columns=['compound_id', 'score'])
                    result = compute_vscreen_metrics(df_pos, df_neg, df_scores)
                    for k, v in result.items():
                        if opt.metrics == "all" or opt.metrics.lower() in k.lower():
                            print(f"{k}: {v:.4f}" if v is not None else f"{k}: Not computed")
                    sys.exit(0)

        # ----------------------------- #
        #  bootstrap on single rank list
        # ----------------------------- #
        if not opt.new: 
                debug_print("Bootstrapping on single rank list...")
                if opt.metrics == "AUC":
                        debug_print("Calling bootstrap_single_keep_ratio with AUC...")
                        bootstrap_single_keep_ratio(ref_scores, lig_list, decoy_list, "AUC")
                elif opt.metrics == "logAUC":
                        debug_print("Calling bootstrap_single_keep_ratio with logAUC...")
                        bootstrap_single_keep_ratio(ref_scores, lig_list, decoy_list, "logAUC")
                elif opt.metrics == "both":
                        debug_print("Calling bootstrap_single_keep_ratio with both AUC and logAUC...")
                        bootstrap_single_keep_ratio(ref_scores, lig_list, decoy_list, "AUC")
                        bootstrap_single_keep_ratio(ref_scores, lig_list, decoy_list, "logAUC")
                elif opt.metrics == "all":
                        debug_print("Calling bootstrap_single_keep_ratio with both AUC and logAUC...")
                        bootstrap_single_keep_ratio(ref_scores, lig_list, decoy_list, "AUC")
                        bootstrap_single_keep_ratio(ref_scores, lig_list, decoy_list, "logAUC")
                        debug_print("Computing additional metrics...")
                        df_pos = pd.DataFrame({'compound_id': lig_list})
                        df_neg = pd.DataFrame({'compound_id': decoy_list})
                        df_scores = pd.DataFrame([[k, float(v[-1])] for k, v in ref_scores.items()],
                                             columns=['compound_id', 'score'])
                        result = compute_vscreen_metrics(df_pos, df_neg, df_scores)
                        for k, v in result.items():
                                if opt.metrics == "all" or opt.metrics.lower() in k.lower():
                                        print(f"{k}: {v:.4f}" if v is not None else f"{k}: Not computed")
                sys.exit(0)

        # -------------------------------- #
        #  bootstrap comparing two methods
        # -------------------------------- #
        debug_print("Bootstrapping to compare two methods...")
        name_list,  ref_list      = [], []
        auc_list,   logauc_list   = [], []
        metrics,    metrics_label = [], []
        mean_list_auc,    medi_list_auc,    rmse_list_auc    = [], [], []
        mean_list_logauc, medi_list_logauc, rmse_list_logauc = [], [], []
        p_ind_list_auc,    p_rel_list_auc,    p_zte_list_auc    = [], [], []
        p_ind_list_logauc, p_rel_list_logauc, p_zte_list_logauc = [], [], []

        new_scores = read_score(opt.new)

        med_auc,    mean_auc,    rmse_auc,    p2_auc,    p4_auc,    p5_auc,    ref_auc,    new_auc, \
        med_logauc, mean_logauc, rmse_logauc, p2_logauc, p4_logauc, p5_logauc, ref_logauc, new_logauc = \
        bootstrap_2methods_keep_ratio(ref_scores, new_scores, lig_list, decoy_list)

        medi_list_auc.append( med_auc)
        mean_list_auc.append(mean_auc)
        rmse_list_auc.append(rmse_auc)

        p_rel_list_auc.append(p2_auc)
        p_ind_list_auc.append(p4_auc)
        p_zte_list_auc.append(p5_auc)

        medi_list_logauc.append( med_logauc)
        mean_list_logauc.append(mean_logauc)
        rmse_list_logauc.append(rmse_logauc)

        p_rel_list_logauc.append(p2_logauc)
        p_ind_list_logauc.append(p4_logauc)
        p_zte_list_logauc.append(p5_logauc)

        names = ["DOCK"]*2*num_bootstrap_replicate
        name_list += names
        # print(len(name_list))

        ref_name = ["REF"]*num_bootstrap_replicate + ["NEW"]*num_bootstrap_replicate
        ref_list += ref_name
        # print(len(ref_list))

        auc_list += ref_auc + new_auc
        logauc_list += ref_logauc + new_logauc

        # metrics += auc_list 
        # metrics += logauc_list
        # print(len(metrics))

        # labels = ["AUC"]*2*num_bootstrap_replicate + ["logAUC"]*2*num_bootstrap_replicate
        # metrics_label += labels
        # print(len(metrics_label))

        data = {'score':   name_list,
                        'method':  ref_list,  
                        'AUC': auc_list,
                        'logAUC':  logauc_list}

        df = pd.DataFrame(data)
        name = 'Bootstrap'

        f = open(name+"_out.csv", 'a')
        f.write(f'#    AUC Mean,             {",".join(str(x) for x in mean_list_auc)}\n')
        f.write(f'#    AUC median,           {",".join(str(x) for x in medi_list_auc)}\n')
        f.write(f'#    AUC RMSE,             {",".join(str(x) for x in rmse_list_auc)}\n')
        f.write(f'#    AUC P-val T-test Rel, {",".join(str(x) for x in p_rel_list_auc)}\n')
        f.write(f'#    AUC P-val T-test Ind, {",".join(str(x) for x in p_ind_list_auc)}\n')
        f.write(f'#    AUC P-val Z-test,     {",".join(str(x) for x in p_zte_list_auc)}\n')

        f.write(f'# logAUC Mean,             {",".join(str(x) for x in mean_list_logauc)}\n')
        f.write(f'# logAUC median,           {",".join(str(x) for x in medi_list_logauc)}\n')
        f.write(f'# logAUC RMSE,             {",".join(str(x) for x in rmse_list_logauc)}\n')
        f.write(f'# logAUC P-val T-test Rel, {",".join(str(x) for x in p_rel_list_logauc)}\n')
        f.write(f'# logAUC P-val T-test Ind, {",".join(str(x) for x in p_ind_list_logauc)}\n')
        f.write(f'# logAUC P-val Z-test,     {",".join(str(x) for x in p_zte_list_logauc)}\n')

        df.to_csv(f, index=False)
        f.close()

        m_auc_list, m_logauc_list = mean_list_auc, mean_list_logauc
        if opt.test == "ttest_rel":
                p_auc_list    = p_rel_list_auc
                p_logauc_list = p_rel_list_logauc
        elif opt.test == "ttest_ind":
                p_auc_list    = p_ind_list_auc
                p_logauc_list = p_ind_list_logauc
        elif opt.test == "ztest":
                p_auc_list    = p_zte_list_auc
                p_logauc_list = p_zte_list_logauc

        violin_plot_2methods(df, name, m_auc_list, p_auc_list, m_logauc_list, p_logauc_list)

def enrichment_standard(scores1, lig_list, decoy_list, metrics):
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



def bootstrap_single_keep_ratio(score_dict, lig_list, decoy_list, metric):
    print(f"Bootstrapping with metric: {metric}")
    all_scores = pd.DataFrame([[k, float(v[-1])] for k, v in score_dict.items()],
                              columns=['compound_id', 'score'])

    df_pos = pd.DataFrame({'compound_id': lig_list})
    df_neg = pd.DataFrame({'compound_id': decoy_list})

    boot_results = {m: [] for m in ['EF1%', 'AUC', 'AUPRC', 'Average_Precision', 'BEDROC_alpha_20']}

    for _ in tqdm(range(num_bootstrap_replicate)):
        sampled_pos = df_pos.sample(frac=1.0, replace=True)
        sampled_neg = df_neg.sample(frac=1.0, replace=True)
        res = compute_vscreen_metrics(sampled_pos, sampled_neg, all_scores)
        for k in boot_results:
            boot_results[k].append(res[k])

    # Convert to DataFrame
    boot_df = pd.DataFrame(boot_results)

    # Plotting
    plt.figure(figsize=(10, 6))
    sns.violinplot(data=boot_df)
    plt.title(f"Bootstrap Distribution of Screening Metrics ({num_bootstrap_replicate} reps)")
    plt.ylabel("Score")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig("bootstrap_metrics.png")
    print("Saved bootstrap plot as 'bootstrap_metrics.png'.")

    # Export CSV
    boot_df.to_csv("bootstrap_metrics.csv", index=False)
    print("Saved bootstrap metrics to 'bootstrap_metrics.csv'.")

    # Print summary stats
    for metric, values in boot_results.items():
        mean_val = np.mean(values)
        std_val = np.std(values)
        print(f"{metric}: Mean = {mean_val:.4f}, Std = {std_val:.4f}")
def bootstrap_2methods_keep_ratio(scores1, scores2, lig_list, decoy_list):

        mols1 = scores1.keys()
        mols2 = scores2.keys()
        common_mols = list( set(mols1) & set(mols2) )

        print(f'\nRef: {len(mols1):d}, New: {len(mols2):d}, Common: {len(common_mols):d}\n')

        # with open(f"not_match_{method_name}", "w") as f:
        #       for ele in list( set(mols1) - set(common_mols) ):
        #               f.write(f"New: {ele}\n")
        #       for ele in list( set(mols2) - set(common_mols) ):
        #               f.write(f"Ref: {ele}\n")


        lig_scores, decoy_scores = [], []
        lig_scores = [scores2[key] for key in list(set(lig_list) & set(common_mols))]
        decoy_scores = [scores2[key] for key in list(set(decoy_list) & set(common_mols))]

        print(f'lig_list:{len(lig_list)}, decoy_list:{len(decoy_list)}')
        print(f'lig_scores:{len(lig_scores)}, decoy_scores:{len(decoy_scores)}')

        # create list for auc/logauc
        list_auc1,    list_auc2,    list_d_auc    = [], [], []
        list_logauc1, list_logauc2, list_d_logauc = [], [], []
        for i in tqdm(range(num_bootstrap_replicate), ncols=100):
                bootstrap_lig_indeces   = np.random.randint(0, len(lig_scores), size=len(lig_scores))
                bootstrap_decoy_indeces = np.random.randint(0, len(decoy_scores), size=len(decoy_scores))
                new_lig_scores   = [ lig_scores[i] for i in bootstrap_lig_indeces]
                new_decoy_scores = [ decoy_scores[i] for i in bootstrap_decoy_indeces]
                new_scores2 = new_lig_scores + new_decoy_scores

                # always bootstrap from score2 and find matched zincID in score1
                # same mols from score2
                mol_names = [ i[0] for i in new_scores2 ] 
                new_scores1 = [ scores1[mol] for mol in mol_names ]

                new_scores1.sort(key=lambda x: float(x[-1]))
                new_scores2.sort(key=lambda x: float(x[-1]))

                points = do_roc(new_scores1, lig_list, decoy_list)
                points = interpolate_curve(points)
                auc1 = AUC(points)*100
                logauc1 = logAUC(points)*100
                # if i == 1:
                #       plot_ROC("plotROC_STD", points, auc1, logauc1)

                points = do_roc(new_scores2, lig_list, decoy_list)
                npoints = interpolate_curve(points)
                auc2 = AUC(npoints)*100
                logauc2 = logAUC(npoints)*100
                # if i == 1:
                #       plot_ROC("plotROC_0.5LD", npoints, auc2, logauc2)

                list_auc1.append(auc1)
                list_auc2.append(auc2)
                list_d_auc.append(auc2-auc1)

                list_logauc1.append(logauc1)
                list_logauc2.append(logauc2)
                list_d_logauc.append(logauc2-logauc1)

        med_auc,  med_logauc  = np.median( list_d_auc ), np.median( list_d_logauc )
        mean_auc, mean_logauc = np.mean( list_d_auc ), np.mean( list_d_logauc )
        rmse_auc, rmse_logauc = rmse(np.array(list_auc2), np.array(list_auc1)), rmse(np.array(list_logauc2), np.array(list_logauc1))

        # print(f'              AUC  Med: {med_auc:.2f}, Mean: {mean_auc:.2f}, RMSE: {rmse_auc:.2f}')
        # print(f'ttest_ind     AUC  P-val: {p_auc:.1E}, T-val: {t_auc:.2f}')
        t1_auc, p1_auc = stats.wilcoxon(list_auc1,    list_auc2)
        t2_auc, p2_auc = stats.ttest_rel(list_auc1,   list_auc2)
        t3_auc, p3_auc = stats.mannwhitneyu(list_auc1,list_auc2)
        t4_auc, p4_auc = stats.ttest_ind(list_auc1,   list_auc2)
        t5_auc, p5_auc = stests.ztest(list_auc1,   x2=list_auc2, value=0, alternative='two-sided')


        print(f'\n\nlogAUC  Med: {med_logauc:.2f}, Mean: {mean_logauc:.2f}, RMSE: {rmse_logauc:.2f}')

        # Paired T-test
        t1_logauc, p1_logauc = stats.wilcoxon(list_d_logauc)
        t2_logauc, p2_logauc = stats.ttest_rel(list_logauc1,list_logauc2)
        # Independent T-test
        t3_logauc, p3_logauc = stats.mannwhitneyu(list_logauc1,list_logauc2)
        t4_logauc, p4_logauc = stats.ttest_ind(list_logauc1,list_logauc2)
        # Independent Z-test 
        t5_logauc, p5_logauc = stests.ztest(list_logauc1,x2=list_logauc2, value=0, alternative='two-sided')

        print(f'Wilcoxon       P-val: {p1_logauc:.1E}, T-val: {t1_logauc:.2f}')
        print(f'ttest_rel      P-val: {p2_logauc:.1E}, T-val: {t2_logauc:.2f}')
        print(f'Mann-WhitneyU  P-val: {p3_logauc:.1E}, T-val: {t3_logauc:.2f}')
        print(f'ttest_ind      P-val: {p4_logauc:.1E}, T-val: {t4_logauc:.2f}')
        print(f'Z-test         P-val: {p5_logauc:.1E}, T-val: {t5_logauc:.2f}')
        return  med_auc, mean_auc, rmse_auc, p2_auc, p4_auc, p5_auc, list_auc1, list_auc2, \
                        med_logauc, mean_logauc, rmse_logauc, p2_logauc, p4_logauc, p5_logauc, list_logauc1, list_logauc2

# median can be mean depends on the returned values from bootstrap function
def violin_plot_2methods(df, name, med_AUC, pVal_AUC, med_logAUC, pVal_logAUC):

        # specific for the website: only comparing two method at once
        num=1
        pos = np.arange(num)

        upperLabels1 = []
        for i in range(num):
                label_text = "$\Delta\mu$= " + f"{med_AUC[i]:.2f}\n"
                if pVal_AUC[i] < 0.05:
                        label_text += "$p$< .05"
                else:
                        label_text += "$p$= " + f"{pVal_AUC[i]:.2f}"
                upperLabels1.append(label_text)
        upperLabels2 = []
        for i in range(num):
                label_text = "$\Delta\mu$= " + f"{med_logAUC[i]:.2f}\n"
                if pVal_logAUC[i] < 0.05:
                        label_text += "$p$< .05"
                else:
                        label_text += "$p$= " + f"{pVal_logAUC[i]:.2f}"
                upperLabels2.append(label_text)


        fig = plt.figure(figsize=(8, 5))
        sns.set_style("whitegrid")
        sns.set_context("notebook", font_scale=1.)

        ax1 = fig.add_subplot(1,2,1)
        ax2 = fig.add_subplot(1,2,2)

        sns.violinplot(x="score", y="AUC", hue="method",
                data=df, split=True, inner='quartile', ax=ax1)
        min_val1 = df["AUC"].min()
        max_val1 = df["AUC"].max()
        add_val1 = (max_val1-min_val1)/5.0
        ax1.set(ylim=(min_val1-add_val1*2.5, max_val1+add_val1*2))
        for tick, label in zip(range(num), ax1.get_xticklabels()):
                ax1.text(pos[tick], min_val1-add_val1*2.2, upperLabels1[tick], horizontalalignment='center')

        sns.violinplot(x="score", y="logAUC", hue="method",
                data=df, split=True, inner='quartile', ax=ax2)
        min_val2 = df["logAUC"].min()
        max_val2 = df["logAUC"].max()
        add_val2 = (max_val2-min_val2)/5.0
        ax2.set(ylim=(min_val2-add_val2*2.5, max_val2+add_val2*2))

        for tick, label in zip(range(num), ax2.get_xticklabels()):
                ax2.text(pos[tick], min_val2-add_val2*2.2, upperLabels2[tick], horizontalalignment='center')

        ax1.set_xlabel("")
        ax2.set_xlabel("")
        fig.tight_layout(pad=2.0)
        fig.suptitle('%s'%name)
        fig.savefig("%s.png"%(name))

def read_ligands_decoys(lig_file, decoy_file):
        ligand_list, decoy_list = [], []
        with open(lig_file, "r") as f:
                for oneline in f:
                        ligand_list.append(oneline.strip())
        with open(decoy_file, "r") as f:
                for oneline in f:
                        decoy_list.append(oneline.strip())
        print("%d ligands and %d decoys read in." % (len(ligand_list), len(decoy_list)) )

        return ligand_list, decoy_list

def read_score(score_file):
        out_dict = collections.OrderedDict()
        with open(score_file, "r") as f:
                for oneline in f:
                        ele = oneline.split()
                        out_dict[ele[0]] = ele # id:line
        return out_dict

def do_roc(scores, lig_list, decoy_list, nbins=10000):

        num_data = len(scores)
        binsize = int(num_data/nbins) + 1
        num_lig = len(lig_list)
        num_dec = len(decoy_list)

        found_ligand = 0
        results = []
        for i in range(num_data):
                if i % binsize == 0:
                        results.append([i-found_ligand, found_ligand])
                # if scores[i][2].split('.')[0] in lig_list: 
                if scores[i][0] in lig_list: 
                        found_ligand += 1
        results.append([num_data - found_ligand, found_ligand])
        results.append([num_dec, num_lig])

        points = []
        for x in results:
                fpr = x[0]*100.0/num_dec
                tpr = x[1]*100.0/num_lig
                points.append([fpr, tpr])
        return points

def interpolate_curve(points):
        i = 0
        while i < len(points) and points[i][0] < 0.1:
                i += 1
        slope = (points[i][1] - points[i-1][1])/(points[i][0] - points[i-1][0])
        intercept = points[i][1] - slope * points[i][0]
        point_one =  [0.100001, (slope * 0.100001 + intercept)]
        npoints = [x for x in points]
        npoints.insert(i, point_one)
        return npoints

def AUC(points):
        """Calulate the area under the curve using trapezoid rule."""
        auc = 0.0
        for point2, point1 in zip(points[1:], points[:-1]):
                #print(point2, point1)
                base = (point2[0] - point1[0]) / 100.0
                height = ( (point2[1] - point1[1])/2.0 + point1[1] ) / 100.0
                auc += (base*height)
        return auc

def logAUC(points):
        # constants
        ## if you modify also change in plots.py        
        LOGAUC_MAX = 1.0   ## this should not change
        LOGAUC_MIN = 0.001 ## this you may want to change if you database is large and you have strong early enrichment. 
        RANDOM_LOGAUC = (LOGAUC_MAX-LOGAUC_MIN)/np.log(10)/np.log10(LOGAUC_MAX/LOGAUC_MIN)


        """Compute semilog x AUC minus the perfectly random semilog AUC."""
        # assumes we have previously interpolated to get y-value at x = 0.1% 
        # generate new points array clamped between 0.1% and 100%

        npoints = []
        for x in points:
                if (x[0] >= LOGAUC_MIN*100) and (x[0] <= LOGAUC_MAX*100):
                        npoints.append( [x[0]/100 , x[1]/100] )

        area = 0.0
        for point2, point1 in zip(npoints[1:], npoints[:-1]):
                if point2[0] - point1[0] < 0.000001:
                        continue

                # segment area computed as integral of log transformed equation
                dx = point2[0]-point1[0]
                dy = point2[1]-point1[1]
                intercept = point2[1] - (dy)/(dx) * point2[0]
                area += dy/np.log(10) + intercept*(np.log10(point2[0])-np.log10(point1[0]))

        # print("logAUC", area, area/np.log10(LOGAUC_MAX/LOGAUC_MIN) - RANDOM_LOGAUC)
        return area/np.log10(LOGAUC_MAX/LOGAUC_MIN) - RANDOM_LOGAUC

def write_points(filename, points, auc):
        """Write curve points to filename."""
        f = open(filename, 'w')
        f.write("#AUC\t%.2f\tLogAUC\t%.2f\n" % (auc[0]*100, auc[1]*100))
        [f.write("%.4f\t%.4f\n" % (x[0], x[1])) for x in points]
        f.close()

def plot_ROC(points, auc, logauc):
        x, y = zip(*points)

        fig = plt.figure(figsize=(10,5), dpi=120, facecolor='w', edgecolor='k')
        fig.subplots_adjust(hspace=0.4, wspace=0.4)
        ax1 = fig.add_subplot(1, 2, 1)
        ax1.plot(x, y, linewidth=1.5, label='AUC: %.2f' % auc)
        ax2 = fig.add_subplot(1, 2, 2)
        ax2.semilogx(x, y, linewidth=1.5, label='logAUC: %.2f' % logauc)
        x = np.arange(0,100,0.1)
        ax1.plot(x, x, 'k--')
        ax1.axis([-0.25, 100, 0, 100])
        ax2.semilogx(x, x, 'k--')
        ax2.axis([0.1, 100, 0, 100])

        ax1.set_xlabel(" Decoys Found %")
        ax1.set_ylabel(" Ligands Found %")
        ax2.set_xlabel("Decoys Found % ")
        ax2.set_ylabel("Ligands Found %")
        ax1.legend(loc="best")
        ax2.legend(loc="best")
        fig.suptitle('ROC Plot')
        fig.savefig("plot_roc.png")
        plt.clf()

def rmse(predictions, targets):
        return np.sqrt(((predictions - targets) ** 2).mean())


if __name__ == '__main__':
        main()

