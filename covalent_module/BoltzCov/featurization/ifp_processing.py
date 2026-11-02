import pandas as pd

##########################
# Residue Interactions
##########################

def one_hot_per_res(fing_file, mechanism):
    """
    One hot encodes per residue interaction per protein target.
    Returns a DataFrame with columns [substance_id, (protID_mechanism_ResID_InteractionType)]
    eg. 1 if X type of interaction occurs at residue A in protein P under mechanism M.
    """
    df_fing = pd.read_csv(fing_file)

    if '5MAJ' in fing_file:
        prot = 'hscpl'
    elif '3F75' in fing_file:
        prot = 'tgcpl'
    else:
        raise ValueError(f"Cannot determine protein from filename: {fing_file}")

    prefix = f"{prot}_{mechanism}"
    df_fing['residue_interaction'] = prefix + '_' + df_fing['residue'] + '_' + df_fing['interaction_type']

    result = (df_fing.groupby(['substance_id', 'residue_interaction'])
                     .size()
                     .unstack(fill_value=0)
                     .clip(upper=1)
                     .reset_index())
    result.columns.name = None
    return result


def interaction_counts(fing_file, mechanism):
    """
    Returns per-interaction-type counts per ligand, prefixed by protein + mechanism.
    Also returns rdkit metadata (ligand-level, not prefixed).
    """
    df = pd.read_csv(fing_file)

    if '5MAJ' in fing_file:
        prot = 'hscpl'
    elif '3F75' in fing_file:
        prot = 'tgcpl'
    else:
        raise ValueError(f"Cannot determine protein from filename: {fing_file}")

    prefix = f"{prot}_{mechanism}"
    df_inter_counts = df[['substance_id', 'saltbridges', 'hbonds', 'pication',
                           'pistack', 'halogen', 'waterbridge', 'hydrophobic', 'metal']].copy()
    df_inter_counts = df_inter_counts.rename(columns={
        c: f"{prefix}_{c}" for c in df_inter_counts.columns if c != 'substance_id'
    })

    # rdkit metadata is ligand-level — no protein/mechanism prefix needed
    metadata_rdkit = df[['substance_id', 'molwt', 'numheavy', 'numrotbonds',
                          'numrings', 'hydrophobicatoms', 'hbondacceptors']].copy()

    return df_inter_counts, metadata_rdkit
