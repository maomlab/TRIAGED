import pandas as pd
from rdkit import Chem
import numpy as np
from rdkit.Chem import rdRGroupDecomposition
from rdkit.Chem import PandasTools
import os

#####################
# R group decomp
#####################

def core_matching(suppl, core1=Chem.MolFromSmarts("n1c(N2CCOCC2)nc(C#N)nc1[*:3]")):
    """Filters supplier for molecules matching the triazine core."""
    core_hits = []
    names = []
    fails = []

    for mol in suppl:
        if mol is None:
            print('warning: could not read molecule')
            continue
        if mol.HasSubstructMatch(core1):
            core_hits.append(mol)
            name = mol.GetProp("Molecule Name")
            names.append(name)
        else:
            fails.append(mol)

    print(f"Core matching: {len(core_hits)} hits, {len(fails)} misses")
    return core_hits, fails, names


def r_decomp(core_hits, names, core1=Chem.MolFromSmarts("n1c(N2CCOCC2)nc(C#N)nc1[*:3]")):
    """
    R group decomposition on core hits.
    Filters names by same success indices as decomp to ensure alignment.
    """
    ps = rdRGroupDecomposition.RGroupDecompositionParameters()
    ps.onlyMatchAtRGroups = False
    ps.removeAllHydrogenRGroups = True
    ps.removeAllHydrogenRGroupsAndLabels = True
    ps.removeHydrogensPostMatch = True

    crgd, fails = rdRGroupDecomposition.RGroupDecompose([core1], core_hits, options=ps, asRows=False)
    print(f"R group decomp failed: {len(fails)}")

    # filter both core_hits and names by same non-fail indices
    success_idx = [i for i in range(len(core_hits)) if i not in fails]
    successful_hits  = [core_hits[i] for i in success_idx]
    successful_names = [names[i] for i in success_idx]

    decomp_df1 = PandasTools.RGroupDecompositionToFrame(crgd, successful_hits)

    # verify alignment before inserting substance_id
    assert len(successful_names) == len(decomp_df1), \
        f"Misalignment: {len(successful_names)} names vs {len(decomp_df1)} rows in decomp"

    decomp_df1.insert(0, "substance_id", successful_names)
    decomp_df1 = decomp_df1.drop(columns=['R2', 'R4', 'R5'], errors='ignore')
    decomp_df1 = decomp_df1.dropna(subset=['R3'])

    return decomp_df1


def convert_mol_smi_df(decomp_df1):
    """Converts RDKit Mol objects in R group columns to SMILES strings."""
    decomp_df1_smi = decomp_df1.copy()

    for col in decomp_df1_smi.columns:
        if col.startswith("R"):
            decomp_df1_smi[col] = decomp_df1_smi[col].apply(
                lambda m: Chem.MolToSmiles(m) if (m is not None and isinstance(m, Chem.Mol)) else None
            )

    return decomp_df1_smi


def decomp_diethyl_amine(decomp_df1_smi, map_num=3):
    """Fragments R3 at the diethylamine N to get R3a and R3b sub-groups."""
    def get_r_groups(smi, map_num=3):
        mol = Chem.MolFromSmiles(smi)

        att_atom = next(a for a in mol.GetAtoms() if a.GetAtomMapNum() == map_num)

        n_atom = next((a for a in att_atom.GetNeighbors() if a.GetSymbol() == 'N'), None)
        if n_atom is None:
            return [None]

        bonds_to_break = [
            mol.GetBondBetweenAtoms(n_atom.GetIdx(), nbr.GetIdx()).GetIdx()
            for nbr in n_atom.GetNeighbors()
            if nbr.GetAtomMapNum() != map_num
        ]

        fragmented = Chem.FragmentOnBonds(mol, bonds_to_break, addDummies=True)
        frags = Chem.GetMolFrags(fragmented, asMols=True)
        return [Chem.MolToSmiles(f) for f in frags]

    diethylamine_grps = {}
    for _, row in decomp_df1_smi.iterrows():
        smi = row['R3']
        name = row['substance_id']
        groups = get_r_groups(smi)
        diethylamine_grps[name] = groups

    return diethylamine_grps


def parse_frags(frag_list):
    """Extracts R3a and R3b from fragment list, excluding the attachment point fragment."""
    legs = [f for f in frag_list if '[*:3]' not in f]
    r3a = legs[0] if len(legs) > 0 else None
    r3b = legs[1] if len(legs) > 1 else None
    return r3a, r3b


def insert_diethyl_df(diethylamine_grps, decomp_df1_smi):
    """Inserts R3a and R3b columns into the decomp dataframe."""
    grps = {}
    for k, v in diethylamine_grps.items():
        if None in v:
            print(f"Skipping {k} — None in fragments: {v}")
        else:
            grps[k] = parse_frags(v)

    frag_df = pd.DataFrame.from_dict(
        grps,
        orient='index',
        columns=['R3a', 'R3b']
    )
    frag_df.index.name = 'substance_id'
    frag_df = frag_df.reset_index()

    all_frags = decomp_df1_smi.merge(frag_df, on='substance_id', how='left')
    all_frags.drop(columns=['R3'], inplace=True)
    return all_frags


def strip_rgroup_labels(series):
    """Strips R group attachment labels from SMILES strings in a Series."""
    series = series.str.replace(r"\[\*\:\d+\]", "[*]", regex=True)
    series = series.str.replace(r"\[\d+\*\]", "[*]", regex=True)
    return series


def count_unique_rgroups(all_frags, rgroup_cols=["R1", "R3a", "R3b"]):
    all_r = pd.concat([all_frags[col] for col in rgroup_cols]).dropna().astype(str)
    stripped = strip_rgroup_labels(all_r)
    counts = stripped.value_counts().reset_index()
    counts.columns = ["R_group", "count"]
    return counts


def count_unique_rgroups_per_position(all_frags, rgroup_cols=["R1", "R3a", "R3b"]):
    group_pos_count = {}
    for col in rgroup_cols:
        stripped = strip_rgroup_labels(all_frags[col].dropna().astype(str))
        counts = stripped.value_counts().reset_index()
        counts.columns = ["R_group", "count"]
        group_pos_count[col] = counts
    return group_pos_count


def assign_ids(group_pos_count):
    """Assigns random unique integer IDs to each R group at each position."""
    R1  = group_pos_count['R1'].copy()
    R3a = group_pos_count['R3a'].copy()
    R3b = group_pos_count['R3b'].copy()

    np.random.seed(42)
    R1['R_id']  = np.random.choice(range(10000, 99999), size=len(R1),  replace=False)
    R3a['R_id'] = np.random.choice(range(10000, 99999), size=len(R3a), replace=False)
    R3b['R_id'] = np.random.choice(range(10000, 99999), size=len(R3b), replace=False)
    return R1, R3a, R3b


def merge_r_ids(R1, R3a, R3b, all_frags):
    """Merges R group IDs into the fragment dataframe by position."""
    for col in ["R1", "R3a", "R3b"]:
        all_frags[col] = strip_rgroup_labels(all_frags[col])

    all_frags = all_frags.merge(
        R1[["R_group", "R_id"]].rename(columns={"R_group": "R1", "R_id": "R1_id"}),
        on="R1", how="left"
    )
    all_frags = all_frags.merge(
        R3a[["R_group", "R_id"]].rename(columns={"R_group": "R3a", "R_id": "R3a_id"}),
        on="R3a", how="left"
    )
    all_frags = all_frags.merge(
        R3b[["R_group", "R_id"]].rename(columns={"R_group": "R3b", "R_id": "R3b_id"}),
        on="R3b", how="left"
    )

    for col in ["R1_id", "R3a_id", "R3b_id"]:
        all_frags[col] = all_frags[col].fillna(0).astype(int)

    return all_frags


def write_rgroup_csvs(r_group_ids_df, R1, R3a, R3b, output_dir):
    """
    Writes two CSVs to output_dir:
    1. compound_rgroup_mapping.csv — per-compound R-group SMILES and IDs
    2. rgroup_lookup.csv — per-position R-group ID to SMILES mapping
    """
    os.makedirs(output_dir, exist_ok=True)

    # per-compound mapping
    mapping = r_group_ids_df[[
    'substance_id',
    'R1',  'R1_id',
    'R3a', 'R3a_id',
    'R3b', 'R3b_id'
    ]].copy()
    mapping = mapping.rename(columns={
        'R3a': 'R2', 'R3a_id': 'R2_id',
        'R3b': 'R3', 'R3b_id': 'R3_id'
    })
    mapping_path = os.path.join(output_dir, 'compound_rgroup_mapping1.csv')
    mapping.to_csv(mapping_path, index=False)
    print(f"  Compound R-group mapping written: {mapping_path}  ({len(mapping)} entries)")

    # per-position lookup tables
    for name, df, label in [('R1', R1, 'R1'), ('R3a', R3a, 'R2'), ('R3b', R3b, 'R3')]:
        out = df[['R_id', 'R_group', 'count']].sort_values('R_id').reset_index(drop=True)
        path = os.path.join(output_dir, f"{label}_id_lookup.csv")
        out.to_csv(path, index=False)
        print(f"  R-group lookup written: {path}  ({len(out)} entries)")


def make_one_hot_r(r_group_ids_df):
    """One hot encodes R group IDs per position. Carries substance_id."""
    design_matrix = pd.get_dummies(
        r_group_ids_df[["R1_id", "R3a_id", "R3b_id"]].astype(str),
        prefix=["R1", "R2", "R3"]
    )
    design_matrix = design_matrix.astype(int)
    design_matrix.insert(0, 'substance_id', r_group_ids_df['substance_id'].values)
    return design_matrix


def run_decomp(suppl, output_dir=None):
    """Full R group decomposition pipeline. Returns one-hot encoded R group matrix."""
    core_hits, _, names = core_matching(suppl)
    decomp_df1          = r_decomp(core_hits, names)
    decomp_df1_smi      = convert_mol_smi_df(decomp_df1)
    diethylamine_grps   = decomp_diethyl_amine(decomp_df1_smi, map_num=3)
    all_frags           = insert_diethyl_df(diethylamine_grps, decomp_df1_smi)
    group_pos_count     = count_unique_rgroups_per_position(all_frags)
    R1, R3a, R3b        = assign_ids(group_pos_count)
    r_group_ids_df      = merge_r_ids(R1, R3a, R3b, all_frags)
    if output_dir is not None:
        write_rgroup_csvs(r_group_ids_df, R1, R3a, R3b, output_dir)
    r_one_hot           = make_one_hot_r(r_group_ids_df)
    return r_one_hot


if __name__ == '__main__':
    suppl = list(Chem.SDMolSupplier('/home/ymanasa/turbo/ymanasa/opt/tgcpl-campaign/martin_ligs/records/martin_ligs.sdf'))
    run_decomp(suppl, output_dir='/home/ymanasa/turbo/ymanasa/opt/tgcpl-campaign/martin_ligs/records')