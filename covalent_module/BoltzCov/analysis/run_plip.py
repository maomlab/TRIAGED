"""
PLIP Batch Analysis for Boltz Predictions with Replicate Structure

Handles directory structure:
    base_dir/
    ├── tgcpl_cov_1/
    │   └── cifs/
    │       └── 20260204_164012_IRUP6_VMCC-0000752.cif
    └── tgcpl_cov_2/
        └── cifs/
            └── ...

- Uses predictions_csv to determine which substance_ids have successful Boltz predictions
- Filters compounds already present in the fingerprints CSV (by substance_id)
- Auto-detects protein type (tgcpl/hscpl) from directory name
- Appends new results to existing fingerprints CSV
- All paths configured via JSON file

Usage:
------
    python plip_batch_analysis.py --json_file config.json

JSON Configuration Format:
--------------------------
{
    "directory":        "/path/to/replicate/root",
    "protein_name":     "3F75",
    "receptor_type":    "protein",
    "outdir":           "/path/to/output",
    "predictions_csv":  "/path/to/records/predictions.csv",
    "selection_method": "first",
    "verbose":          true
}

Authors: Manasa Yadavalli
"""

import os
import re
import csv
import sys
import json
import shutil
import argparse
import tempfile
from collections import defaultdict

import pandas as pd
import gemmi
from plip.basic import config
from plip.structure.preparation import PDBComplex


# ── JSON config ────────────────────────────────────────────────────────────────

def read_json_args(json_file):
    """
    Reads JSON config and returns an argparse.Namespace.

    Required keys: directory, protein_name, receptor_type, outdir, predictions_csv
    Optional keys: selection_method (default: 'first'), verbose (default: false)
    """
    with open(json_file, 'r') as jf:
        cfg = json.load(jf)

    required = ['directory', 'protein_name', 'receptor_type', 'outdir', 'predictions_csv']
    missing = [k for k in required if not cfg.get(k)]
    if missing:
        raise ValueError(f"Missing required JSON keys: {', '.join(missing)}")

    valid_receptor_types = ['dna', 'rna', 'protein']
    if cfg['receptor_type'] not in valid_receptor_types:
        raise ValueError(f"receptor_type must be one of {valid_receptor_types}, got: {cfg['receptor_type']}")

    valid_selection_methods = ['first', 'random', 'best_confidence']
    selection_method = cfg.get('selection_method', 'first')
    if selection_method not in valid_selection_methods:
        raise ValueError(f"selection_method must be one of {valid_selection_methods}, got: {selection_method}")

    return argparse.Namespace(
        directory=os.path.expandvars(cfg['directory']),
        protein_name=cfg['protein_name'],
        receptor_type=cfg['receptor_type'],
        outdir=os.path.expandvars(cfg['outdir']),
        predictions_csv=os.path.expandvars(cfg['predictions_csv']),
        selection_method=selection_method,
        verbose=cfg.get('verbose', False),
    )


# ── Helpers ────────────────────────────────────────────────────────────────────

# Mapping between short name (subdir filter) and PDB ID (predictions.csv protein col)
PROTEIN_ALIASES = {
    'tgcpl': '3F75',
    '3f75':  'tgcpl',
    'hscpl': '5MAJ',
    '5maj':  'hscpl',
}

def resolve_protein_names(protein_name):
    """
    Accepts either form (tgcpl/hscpl or 3F75/5MAJ) and returns
    (subdir_filter, pdb_id) for use in subdir matching and predictions filtering.

    Examples:
        'tgcpl' -> ('tgcpl', '3F75')
        '3F75'  -> ('tgcpl', '3F75')
        'hscpl' -> ('hscpl', '5MAJ')
        '5MAJ'  -> ('hscpl', '5MAJ')
    """
    lookup = protein_name.lower()
    alias = PROTEIN_ALIASES.get(lookup)
    if alias is None:
        raise ValueError(
            f"Unknown protein_name: '{protein_name}'. "
            f"Use one of: tgcpl, 3F75, hscpl, 5MAJ"
        )
    # Determine which is the subdir filter and which is the PDB ID
    if lookup in ('tgcpl', 'hscpl'):
        subdir_filter = lookup
        pdb_id = alias
    else:
        subdir_filter = alias
        pdb_id = protein_name.upper()
    return subdir_filter, pdb_id


def get_predicted_substance_ids(records, pdb_id, COVALENT):
    """
    Returns set of substance_ids from predictions_csv that have a successful
    Boltz prediction for the given protein (matched by PDB ID e.g. 3F75).
    """
    if COVALENT: 
        predictions_csv = os.path.join(records, 'predictions.csv')
    else:
        predictions_csv = os.path.join(records, 'noncov_predictions.csv')

    if not os.path.exists(predictions_csv):
        raise FileNotFoundError(f"predictions_csv not found: {predictions_csv}")

    df = pd.read_csv(predictions_csv)

    if 'substance_id' not in df.columns:
        raise ValueError(f"'substance_id' column not found in {predictions_csv}")

    if 'protein' in df.columns:
        df = df[df['protein'] == pdb_id]
        print(f"Predictions for {pdb_id}: {len(df)}")
    else:
        print(f"[WARNING] No 'protein' column in predictions_csv — using all {len(df)} rows.")

    return set(df['substance_id'].astype(str))


def get_already_processed(fingerprint_csv):
    """
    Returns set of substance_ids already in fingerprints CSV.
    Returns empty set if file does not exist.
    """
    if not os.path.exists(fingerprint_csv):
        print(f"[INFO] Fingerprints CSV not found, will create: {fingerprint_csv}")
        return set()
    df = pd.read_csv(fingerprint_csv)
    if 'substance_id' not in df.columns:
        print(f"[WARNING] 'substance_id' column not found in {fingerprint_csv}, processing all compounds.")
        return set()
    already = set(df['substance_id'].astype(str))
    print(f"Already processed: {len(already)} compounds (will be skipped)")
    return already


def extract_substance_id_from_filename(filename):
    """
    Extracts substance_id from CIF filename.
    Expected format: 20260204_164012_IRUP6_VMCC-0000752.cif
    Returns the last underscore-separated token before .cif (e.g. VMCC-0000752).
    """
    basename = os.path.splitext(os.path.basename(filename))[0]
    return basename.split('_')[-1]


def extract_ligand_code(filepath):
    """
    Extracts first 3 chars of ligand code from filename for PLIP interaction set lookup.
    From: PROT_LIG_model_0.cif -> first 3 chars of LIG. Falls back to 'UNK'.
    """
    basename = os.path.basename(filepath)
    match = re.search(r'_([A-Za-z0-9]+)_model_0\.cif$', basename)
    if match:
        return match.group(1)[:3]
    return "UNK"


def select_representative_structure(cif_files, selection_method='first'):
    """
    Selects one representative CIF from multiple replicates.
    Methods: 'first' (alphabetical), 'random', 'best_confidence' (falls back to first).
    """
    if len(cif_files) == 1:
        return cif_files[0]
    if selection_method == 'random':
        import random
        return random.choice(cif_files)
    return sorted(cif_files)[0]


def append_to_csv(csv_path, headers, rows):
    """
    Appends rows to CSV. Creates file with headers if it does not exist.
    """
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    file_exists = os.path.exists(csv_path)
    with open(csv_path, 'a', newline='') as csvfile:
        writer = csv.writer(csvfile)
        if not file_exists:
            writer.writerow(headers)
        writer.writerows(rows)


# ── CIF/PDB conversion ─────────────────────────────────────────────────────────

def convert_cif_to_pdb(cif_path, pdb_path, chain_map=None):
    """Converts mmCIF to PDB using gemmi, fixing long chain names."""
    doc = gemmi.cif.read_file(cif_path)
    block = doc.sole_block()
    structure = gemmi.make_structure_from_block(block)
    for model in structure:
        for chain in model:
            if len(chain.name) != 1:
                new_chain = chain_map.get(chain.name, 'X') if chain_map else 'X'
                chain.name = new_chain
    structure.write_pdb(pdb_path)


def convert_pdb_to_pdb(source_pdb_path, target_pdb_path):
    """Standardizes/reformats a PDB file via gemmi."""
    structure = gemmi.read_structure(source_pdb_path)
    structure.write_pdb(target_pdb_path)


# ── PLIP interaction extraction ────────────────────────────────────────────────

def get_interactions(interactions):
    """Returns counts of each interaction type as a list."""
    return [
        len(interactions.saltbridge_lneg + interactions.saltbridge_pneg),
        len(interactions.hbonds_ldon + interactions.hbonds_pdon),
        len(interactions.pication_laro + interactions.pication_paro),
        len(interactions.pistacking),
        len(interactions.halogen_bonds),
        len(interactions.water_bridges),
    ]


def get_interacting_residues(interactions):
    """Returns dict of {restype+resnr: interaction_type} for all interactions."""
    residues = {}
    for sb in interactions.saltbridge_lneg + interactions.saltbridge_pneg:
        residues[f'{sb.restype}{sb.resnr}'] = 'saltbridge'
    for hb in interactions.hbonds_ldon + interactions.hbonds_pdon:
        residues[f'{hb.restype}{hb.resnr}'] = 'hbond'
    for pc in interactions.pication_laro + interactions.pication_paro:
        residues[f'{pc.restype}{pc.resnr}'] = 'pication'
    for ps in interactions.pistacking:
        residues[f'{ps.restype}{ps.resnr}'] = 'pistack'
    for ha in interactions.halogen_bonds:
        residues[f'{ha.restype}{ha.resnr}'] = 'halogen'
    for wb in interactions.water_bridges:
        residues[f'{wb.restype}{wb.resnr}'] = 'waterbridge'
    return residues


# ── CIF file discovery ─────────────────────────────────────────────────────────

def find_all_cif_files_in_replicates(root_dir, protein_name, predicted_ids, already_processed, VERBOSE):
    """
    Finds CIF files across replicate subdirectories.
    Only includes substance_ids that:
      1. Are in predicted_ids (have a successful Boltz prediction)
      2. Are NOT in already_processed (not yet in fingerprints CSV)

    Expected structure:
        root_dir/
        ├── <protein>_rep_1/
        │   └── cifs/
        │       └── date_hash_runID_SUBSTANCEID.cif
        └── <protein>_rep_2/
            └── cifs/
                └── ...

    :return: Dict mapping substance_id -> list of CIF file paths
    """
    print(f"Searching for CIF files in {root_dir}...")

    compound_files = defaultdict(list)

    if not os.path.isdir(root_dir):
        print(f"Directory not found: {root_dir}")
        return compound_files

    prot_root = os.path.join(root_dir, protein_name)
    subdirs = [d for d in os.listdir(prot_root) if os.path.isdir(os.path.join(prot_root, d))]
    if not subdirs:
        print(f"No subdirectories found in {prot_root}")
        return compound_files

    print(f"Found {len(subdirs)} subdirectories: {', '.join(subdirs[:3])}{'...' if len(subdirs) > 3 else ''}")

    for subdir in subdirs: # rep1, rep2.. 
        if subdir == 'rep1':
            cifs_dir = os.path.join(prot_root, subdir, 'cifs')
            if not os.path.isdir(cifs_dir):
                if VERBOSE:
                    print(f"   {subdir}: no cifs/ subdirectory, skipping.")
                continue
            cif_files = [f for f in os.listdir(cifs_dir) if f.endswith('.cif')]
            if VERBOSE:
                print(f"   {subdir}: {len(cif_files)} CIF files")

            for cif_file in cif_files:
                substance_id = extract_substance_id_from_filename(cif_file)

                if substance_id not in predicted_ids:
                    if VERBOSE:
                        print(f"      Skipping {substance_id} (not in predictions.csv)")
                    continue

                if substance_id in already_processed:
                    if VERBOSE:
                        print(f"      Skipping {substance_id} (already in fingerprints CSV)")
                    continue

                compound_files[substance_id].append(os.path.join(cifs_dir, cif_file))

    total = sum(len(v) for v in compound_files.values())
    print(f"Found {total} CIF files for {len(compound_files)} new compounds to process")
    return compound_files


# ── Main ───────────────────────────────────────────────────────────────────────

def main(args):
    fp_headers = [
        "substance_id", "name", "smiles", "inchi", "molwt",
        "numheavy", "numrotbonds", "numrings",
        "hydrophobicatoms", "hbondacceptors",
        "saltbridges", "hbonds", "pication",
        "pistack", "halogen", "waterbridge"
    ]
    res_headers = ["substance_id", "name", "residue", "interaction_type"]

    VERBOSE = args.verbose
    COVALENT = args.COVALENT
    # Receptor config
    if args.receptor_type == "protein":
        config.DNARECEPTOR = False
    elif args.receptor_type in ["rna", "dna"]:
        config.DNARECEPTOR = True
    else:
        print(f"Invalid receptor_type: {args.receptor_type}. Use: protein, rna, dna")
        sys.exit(1)

    # Resolve protein_name to subdir_filter (tgcpl/hscpl) and pdb_id (3F75/5MAJ)
    subdir_filter, pdb_id = resolve_protein_names(args.protein_name)
    if COVALENT:
        fingerprint_csv = os.path.join(args.records, f'{pdb_id}_ifps.csv')
        residues_csv = os.path.join(args.records, f'{pdb_id}_ifps_residues.csv')
    else:
        fingerprint_csv = os.path.join(args.records, f'{pdb_id}_noncov_ifps.csv')
        residues_csv = os.path.join(args.records, f'{pdb_id}_noncov_ifps_residues.csv')

    os.makedirs(args.outdir, exist_ok=True)

    print(f"Protein: {subdir_filter.upper()} ({pdb_id})")
    print(f"Fingerprints CSV: {fingerprint_csv}")

    # Load substance_ids with successful predictions (filter by PDB ID)
    predicted_ids = get_predicted_substance_ids(args.records, pdb_id, args.COVALENT)
    print(f"Substance IDs with predictions: {len(predicted_ids)}")

    # Load substance_ids already in fingerprints CSV
    already_processed = get_already_processed(fingerprint_csv)

    # Find CIF files to process (filter subdirs with reps)
    compound_files = find_all_cif_files_in_replicates(
        args.directory, subdir_filter, predicted_ids, already_processed, VERBOSE
    )

    if not compound_files:
        print("No new compounds to process.")
        return []

    print(f"\nAnalyzing {len(compound_files)} compounds...")

    temp_dir = tempfile.mkdtemp()
    collected_data = []
    residue_data = []
    errors = []

    for substance_id, cif_files in sorted(compound_files.items()):
        target_file = select_representative_structure(cif_files, args.selection_method)

        if VERBOSE:
            print(f"\nProcessing {substance_id} ({len(cif_files)} replicates)")
            print(f"  Using: {target_file}")

        name = os.path.splitext(os.path.basename(target_file))[0]
        lig_code = extract_ligand_code(target_file)
        filetype = target_file.rsplit('.', 1)[-1].lower()
        target_pdb = os.path.join(temp_dir, f"{substance_id}.pdb")

        try:
            if filetype == 'cif':
                convert_cif_to_pdb(target_file, target_pdb)
            elif filetype == 'pdb':
                convert_pdb_to_pdb(target_file, target_pdb)
            else:
                print(f"  Unsupported file type: {filetype}")
                errors.append(target_file)
                continue

            my_mol = PDBComplex()
            my_mol.load_pdb(target_pdb)
            my_mol.analyze()

            interaction_key = f"{lig_code}:X:1"
            if interaction_key not in my_mol.interaction_sets:
                if my_mol.interaction_sets:
                    interaction_key = list(my_mol.interaction_sets.keys())[0]
                    if VERBOSE:
                        print(f"  Using interaction set: {interaction_key}")
                else:
                    print(f"  No interactions found for {substance_id}")
                    errors.append(target_file)
                    continue

            interactions = my_mol.interaction_sets[interaction_key]

            collected_data.append([
                substance_id,
                name,
                interactions.ligand.smiles.strip() if interactions.ligand.smiles else "",
                interactions.ligand.inchikey.strip() if interactions.ligand.inchikey else "",
                interactions.ligand.molweight,
                interactions.ligand.heavy_atoms,
                interactions.ligand.num_rot_bonds,
                interactions.ligand.num_rings,
                len(interactions.ligand.hydroph_atoms),
                len(interactions.ligand.hbond_acc_atoms),
                *get_interactions(interactions),
            ])

            for resid, int_type in get_interacting_residues(interactions).items():
                residue_data.append([substance_id, name, resid, int_type])

            if VERBOSE:
                print(f"  Done: {substance_id}")

        except Exception as e:
            print(f"  Error processing {substance_id}: {e}")
            errors.append(target_file)
            continue

    shutil.rmtree(temp_dir)

    # Append results
    if collected_data:
        print(f"\nAppending {len(collected_data)} new results to {fingerprint_csv}...")
        append_to_csv(fingerprint_csv, fp_headers, collected_data)
        append_to_csv(residues_csv, res_headers, residue_data)
    else:
        print("\nNo new results to write.")

    print(f"\nDone.")
    print(f"   Processed: {len(collected_data)} compounds")
    print(f"   Failed:    {len(errors)} structures")
    print(f"   Output:    {fingerprint_csv}")

    return errors


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run PLIP on Boltz predictions. Skips already-processed compounds."
    )
    parser.add_argument("--json_file", required=True, help="Path to JSON config file.")
    cli_args = parser.parse_args()

    args = read_json_args(cli_args.json_file)
    main(args)