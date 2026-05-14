import pandas as pd
from rdkit import Chem
import os
import argparse
import r_decomp
import get_preds
import get_files
import ifp_processing

# This script only collates all possible feature data.
# Does NOT perform ablation or linear/ML studies — that is handled downstream.

#####################
# Build Feature Matrix
#####################

def build_feature_matrix(
    substance_ids,
    exp_ic50,
    exp_sele,
    ifp_res_blocks,
    ifp_count_blocks,
    pred_blocks,
    conf_blocks,
    metadata_df,
    metadata_rdkit,
    one_hot_r,
    metadata_feats,
):
    """
    Merges ALL available feature blocks on substance_id via left join.
    substance_ids from exp_ic50 serves as the anchor — all blocks aligned to it.
    Ablation/column exclusion is handled downstream in a separate script.
    """
    base = pd.DataFrame({'substance_id': substance_ids})
    blocks = []

    # experimental data
    blocks.append(exp_ic50)
    blocks.append(exp_sele)

    # interaction fingerprints
    blocks.extend(ifp_res_blocks)
    blocks.extend(ifp_count_blocks)

    # boltz predictions and confidence
    blocks.extend(pred_blocks)
    blocks.extend(conf_blocks)

    # physicochemical metadata — filter to requested columns only
    valid_meta = [f for f in metadata_feats if f in metadata_df.columns]
    missing_meta = [f for f in metadata_feats if f not in metadata_df.columns]
    if missing_meta:
        print(f"Warning: metadata columns not found and skipped: {missing_meta}")
    blocks.append(metadata_df[['substance_id'] + valid_meta])

    # rdkit-derived ligand metadata
    blocks.append(metadata_rdkit)

    # r group one-hot encoding
    blocks.append(one_hot_r)

    X = base.copy()
    for block in blocks:
        X = X.merge(block, on='substance_id', how='left')

    X = X.set_index('substance_id')
    return X


#####################
# Main
#####################

def main(args):
    import ipdb; ipdb.set_trace()
    cfg = get_files.read_json_args(args.json_file)

    record_path            = cfg['record_path']
    metadata_feats         = cfg['metadata_feats']
    boltz_pred_feats       = cfg['boltz_pred_feats']
    boltz_confidence_feats = cfg['boltz_confidence_feats']
    sdf_file               = cfg['sdf_file']

    mech_files = get_files.get_mechanism_files(cfg)

    all_ifp_res    = []
    all_ifp_counts = []
    all_preds      = []
    all_confs      = []
    metadata_rdkit = None  # ligand-level, only need once

    for mech, files in mech_files.items():
        print(f"\nLoading {mech} data...")

        # per-residue binary fingerprints
        all_ifp_res.append(ifp_processing.one_hot_per_res(files['tgcpl_ifp_res'], mech))
        all_ifp_res.append(ifp_processing.one_hot_per_res(files['hscpl_ifp_res'], mech))

        # interaction type counts
        tgcpl_counts, rdkit_meta = ifp_processing.interaction_counts(files['tgcpl_ifp_type'], mech)
        hscpl_counts, _          = ifp_processing.interaction_counts(files['hscpl_ifp_type'], mech)
        all_ifp_counts.extend([tgcpl_counts, hscpl_counts])

        # rdkit metadata is ligand-level — only take from first mechanism
        if metadata_rdkit is None:
            metadata_rdkit = rdkit_meta

        # boltz predictions and confidence
        all_preds.append(get_preds.get_pred_only(files['tgcpl'], boltz_confidence_feats, 'tgcpl', mech))
        all_preds.append(get_preds.get_pred_only(files['hscpl'], boltz_confidence_feats, 'hscpl', mech))
        all_confs.append(get_preds.get_conf_only(files['tgcpl'], boltz_pred_feats, 'tgcpl', mech))
        all_confs.append(get_preds.get_conf_only(files['hscpl'], boltz_pred_feats, 'hscpl', mech))

    # r group decomposition
    print("\nRunning R group decomposition...")
    suppl     = Chem.SDMolSupplier(sdf_file, removeHs=False)
    r_one_hot = r_decomp.run_decomp(suppl)

    # experimental data — substance_id from here is the anchor
    print("\nLoading experimental data...")
    exp_ic50, exp_sele = get_preds.get_exp_data(os.path.join(record_path, 'experiment_readouts.csv'))

    # physicochemical metadata
    metadata_df = pd.read_csv(os.path.join(record_path, 'metadata.csv'))

    print("\nBuilding feature matrix...")
    X = build_feature_matrix(
        substance_ids    = exp_ic50['substance_id'],
        exp_ic50         = exp_ic50,
        exp_sele         = exp_sele,
        ifp_res_blocks   = all_ifp_res,
        ifp_count_blocks = all_ifp_counts,
        pred_blocks      = all_preds,
        conf_blocks      = all_confs,
        metadata_df      = metadata_df,
        metadata_rdkit   = metadata_rdkit,
        one_hot_r        = r_one_hot,
        metadata_feats   = metadata_feats,
    )

    outpath = os.path.join(record_path, 'feature_matrix.csv')
    X.to_csv(outpath)
    print(f"\nFeature matrix shape: {X.shape}")
    print(f"Saved to: {outpath}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Featurizer: collates all feature data for downstream SHAP/ablation analysis."
    )
    parser.add_argument(
        "--json_file",
        required=True,
        help="Path to JSON configuration file.",
    )
    args = parser.parse_args()
    main(args)