import pandas as pd
#####################
# Boltz Preds/Confs
#####################

def get_pred_only(avg_preds_path, boltz_confidence_feats, prot, mechanism):
    """
    Loads averaged boltz predictions, drops confidence cols and protein col.
    Prefixes remaining columns with protein + mechanism.
    """
    predictions_avg = pd.read_csv(avg_preds_path)
    drop_cols = ['protein'] + (boltz_confidence_feats or [])
    df = predictions_avg.drop(columns=[c for c in drop_cols if c in predictions_avg.columns])
    prefix = f"{prot}_{mechanism}"
    df = df.rename(columns={c: f"{prefix}_{c}" for c in df.columns if c != 'substance_id'})
    return df


def get_conf_only(avg_preds_path, boltz_pred_feats, prot, mechanism):
    """
    Loads averaged boltz predictions, drops pred cols and protein col.
    Prefixes remaining confidence columns with protein + mechanism.
    """
    predictions_avg = pd.read_csv(avg_preds_path)
    drop_cols = ['protein'] + (boltz_pred_feats or [])
    df = predictions_avg.drop(columns=[c for c in drop_cols if c in predictions_avg.columns])
    prefix = f"{prot}_{mechanism}"
    df = df.rename(columns={c: f"{prefix}_{c}" for c in df.columns if c != 'substance_id'})
    return df


def get_exp_data(exp_readouts):
    """
    Loads experimental IC50 and selectivity data.
    Returns copies to avoid SettingWithCopyWarning downstream.
    """
    exp_reads = pd.read_csv(exp_readouts)
    exp_ic50 = exp_reads[['substance_id', 'mean_tgcpl_log_ic50 (uM)', 'mean_hscpl_log_ic50 (uM)']].copy()
    exp_sele = exp_reads[['substance_id', 'selectivity']].copy()
    return exp_ic50, exp_sele
