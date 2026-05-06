'''
Boltz Covalent Docking Pipeline
Author: Manasa Yadavalli | Env: cdd_pkl

Orchestrates covalent/non-covalent docking with Boltz-2 via SLURM.
Ligands can be sourced from CDD Vault or provided locally.

Ligand input (choose one):
  CDD mode   — CDD_API_KEY + VAULT_ID + READOUT_QUERY + MOL_QUERY
  Local mode — SMILES_INPUT (single SMILES) or SMILES_CSV (CSV with [substance_id, smiles])

Usage:
    python main.py --json_file config.json
'''
import os
import pandas as pd
import json
import shutil
import time
import argparse
from BoltzCov.update_ligands import update_predictions
from BoltzCov.run_boltz import pull_boltz2_weights

# ── helpers ───────────────────────────────────────────────────────────────────

_REQUIRED_CSV_COLS = {'substance_id', 'smiles'}


def _load_local_ligands(smiles_input: str | None, smiles_csv: str | None) -> pd.DataFrame:
    """
    Build a ligand DataFrame from local inputs (no CDD required).

    Accepts either:
      - smiles_input : a bare SMILES string  → single-row DataFrame
      - smiles_csv   : path to a CSV with at minimum [substance_id, smiles]
                       optional columns: inchi_key, name

    Returns a DataFrame with columns [substance_id, smiles] and optionally
    [inchi_key] – the same columns consumed downstream by the docking step.
    """
    if smiles_csv:
        smiles_csv = os.path.expandvars(smiles_csv)
        if not os.path.exists(smiles_csv):
            raise FileNotFoundError(f"SMILES_CSV not found: {smiles_csv}")
        df = pd.read_csv(smiles_csv)
        missing = _REQUIRED_CSV_COLS - set(df.columns)
        if missing:
            raise ValueError(
                f"SMILES_CSV is missing required columns: {missing}. "
                f"Found: {list(df.columns)}"
            )
        # keep only the columns the pipeline cares about
        keep = ['substance_id', 'smiles'] + [
            c for c in ('inchi_key', 'name') if c in df.columns
        ]
        return df[keep].reset_index(drop=True)

    if smiles_input:
        from rdkit import Chem
        from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey
        mol = Chem.MolFromSmiles(smiles_input)
        if mol is None:
            raise ValueError(f"Could not parse SMILES string: {smiles_input}")
        inchi = MolToInchi(mol)
        inchi_key = InchiToInchiKey(inchi)
        return pd.DataFrame([{'substance_id': 'LIG_0001', 'smiles': smiles_input, 'inchi': inchi, 'inchi_key': inchi_key}])
    
    raise ValueError(
        "Local mode requires either SMILES_INPUT (a SMILES string) or "
        "SMILES_CSV (path to a CSV with columns [substance_id, smiles])."
    )

# ── config reader ─────────────────────────────────────────────────────────────

def read_json_args(json_file):
    """
    Reads the JSON input file with arguments required to dock with Boltz-2.

    Returns a dict so callers are not position-sensitive.
    """
    with open(json_file, 'r') as jf:
        a = json.load(jf)

    def ep(key, default=None):
        """expandvars helper that also handles None."""
        val = a.get(key, default)
        return os.path.expandvars(val) if isinstance(val, str) else val

    return dict(
        record_path     = ep("RECORD_PATH"),
        # CDD (all optional)
        cdd_api_key     = a.get("CDD_API_KEY"),
        vault_id        = a.get("VAULT_ID"),
        readout_query   = a.get("READOUT_QUERY"),
        mol_query       = a.get("MOL_QUERY"),
        syn_include     = a.get("SYN_INCLUDE"),
        syn_exclude     = a.get("SYN_EXCLUDE"),
        # Local SMILES (alternative to CDD)
        smiles_input    = a.get("SMILES_INPUT"),
        smiles_csv      = ep("SMILES_CSV"),
        # Structure / docking
        pdb             = ep("PDB"),
        res_idx         = a.get("RES_IDX"),
        ligand_chain    = a.get("LIGAND_CHAIN"),
        msa_path        = ep("MSA_PATH"),
        run_cache       = ep("RUN_CACHE"),
        boltz_cache     = ep("BOLTZ_CACHE"),
        slurm_template  = ep("SLURM_TEMPLATE"),
        min_replicates  = int(a.get("min_replicates") or 3),
        VERBOSE         = a.get("VERBOSE", False),
        output_dir      = ep("OUTPUT"),
        COVALENT        = a.get("COVALENT", False),
    )


# ── validation ────────────────────────────────────────────────────────────────

def _validate_args(cfg: dict) -> None:
    """
    Raise ValueError for any combination of missing required arguments.
    Handles CDD mode, local mode, covalent vs non-covalent.
    """
    cdd_mode = bool(cfg['cdd_api_key'] and cfg['vault_id'])
    local_mode = bool(cfg['smiles_input'] or cfg['smiles_csv'])

    if not cdd_mode and not local_mode:
        raise ValueError(
            "No ligand source provided. Supply either:\n"
            "  • CDD_API_KEY + VAULT_ID  (CDD Vault mode)\n"
            "  • SMILES_INPUT            (single SMILES string)\n"
            "  • SMILES_CSV              (path to CSV with [substance_id, smiles])"
        )

    if cdd_mode and local_mode:
        raise ValueError(
            "Ambiguous ligand source: both CDD credentials and local SMILES input "
            "were supplied. Please provide only one."
        )

    # params required regardless of mode
    shared = {k: cfg[k] for k in ('pdb', 'boltz_cache', 'run_cache', 'slurm_template', 'output_dir')}

    # extra params required only for covalent docking
    if cfg['COVALENT']:
        shared.update({k: cfg[k] for k in ('res_idx', 'ligand_chain')})

    # CDD mode needs query dicts too
    if cdd_mode:
        shared.update({k: cfg[k] for k in ('readout_query', 'mol_query')})

    missing = [k for k, v in shared.items() if v is None]
    if missing:
        raise ValueError(f"Missing required arguments: {', '.join(missing)}")


# ── main ──────────────────────────────────────────────────────────────────────

def main(args):
    """
    Top priority:
    1) Pull ligands – from CDD Vault *or* from a local SMILES / CSV input
    2) Update compound records
    3) Dock ligands not previously attempted
    4) Update relevant results in compound records
    """
    cfg = read_json_args(args.json_file)
    _validate_args(cfg)

    # unpack for readability
    record_path    = cfg['record_path']
    cdd_api_key    = cfg['cdd_api_key']
    vault_id       = cfg['vault_id']
    readout_query  = cfg['readout_query']
    mol_query      = cfg['mol_query']
    syn_include    = cfg['syn_include']
    syn_exclude    = cfg['syn_exclude']
    smiles_input   = cfg['smiles_input']
    smiles_csv     = cfg['smiles_csv']
    pdb            = cfg['pdb']
    res_idx        = cfg['res_idx']
    ligand_chain   = cfg['ligand_chain']
    msa_path       = cfg['msa_path']
    run_cache      = cfg['run_cache']
    boltz_cache    = cfg['boltz_cache']
    slurm_template = cfg['slurm_template']
    min_replicates = cfg['min_replicates']
    VERBOSE        = cfg['VERBOSE']
    output_dir     = cfg['output_dir']
    COVALENT       = cfg['COVALENT']

    cdd_mode = bool(cdd_api_key and vault_id)

    # ── cache / weight setup ──────────────────────────────────────────────────
    if os.path.exists(run_cache):
        print("[WARNING] Run cache exists and will be deleted.")
        print("Cancel in 5 seconds to prevent deletion.")
        time.sleep(5)
        shutil.rmtree(run_cache)
    os.makedirs(run_cache)

    if not os.path.exists(boltz_cache):
        os.makedirs(boltz_cache)
        from pathlib import Path
        print("[WARNING] Boltz weights not found. Downloading weights and PDB ligand files…")
        pull_boltz2_weights.download_boltz2(Path(boltz_cache))

    # ── step 1: acquire ligands ───────────────────────────────────────────────
    if cdd_mode:
        from BoltzCov.update_ligands import pull_cdd_ligs
        print(
            "1. Pulling ligands from CDD Vault using the following queries:\n"
            f"   Readout query : {readout_query}\n"
            f"   Molecule query: {mol_query}"
        )
        readouts, molecules = pull_cdd_ligs.cdd_query(
            API_KEY=cdd_api_key,
            VAULT_ID=vault_id,
            readout_query=readout_query,
            mol_query=mol_query,
        )
        new_readouts_df, new_metadata_df = pull_cdd_ligs.get_ic50s(
            readouts, molecules,
            syn_include=syn_include,
            syn_exclude=syn_exclude,
        )
    else:
        source = cfg['smiles_csv'] or "SMILES_INPUT"
        print(f"1. Loading ligands from local source: {source}")
        new_metadata_df = _load_local_ligands(smiles_input, smiles_csv)
        new_readouts_df = None  # not available in local mode
        if VERBOSE:
            print(f"   Loaded {len(new_metadata_df)} compound(s).")

    mode_label = "Covalent" if COVALENT else "Non-Covalent"

    if cdd_mode:
        # ── step 2 (CDD only): update compound records ────────────────────────
        print("2. Updating compound records.")

        if record_path is None:
            record_path = os.path.join(output_dir, 'records')
            if VERBOSE:
                print(f"   No RECORD_PATH supplied – writing records to {record_path}")
        os.makedirs(record_path, exist_ok=True)

        old_meta = os.path.join(record_path, 'metadata.csv')
        old_exp  = os.path.join(record_path, 'experiment_readouts.csv')

        meta_empty = not os.path.exists(old_meta) or os.path.getsize(old_meta) <= 1
        exp_empty  = not os.path.exists(old_exp)  or os.path.getsize(old_exp)  <= 1

        if meta_empty or exp_empty:
            if VERBOSE:
                print(
                    f"   metadata.csv or experiment_readouts.csv not found in {record_path}. "
                    "Writing fresh records."
                )
            print("   You have 5 seconds to cancel to avoid overwriting existing records.")
            time.sleep(5)
            pd.DataFrame(new_readouts_df).to_csv(old_exp,  index=False)
            pd.DataFrame(new_metadata_df).to_csv(old_meta, index=False)
            if VERBOSE:
                print(f"   Fresh records written to {record_path}.")
        else:
            print("   You have 5 seconds to cancel to avoid overwriting existing records.")
            time.sleep(5)
            if VERBOSE:
                print(f"   Merging new CDD data into existing records in {record_path}.")
            old_meta_df = pd.read_csv(old_meta)
            old_exp_df  = pd.read_csv(old_exp)
            updated_metadata, updated_readouts = pull_cdd_ligs.update_local_data(
                old_meta_df, old_exp_df, new_metadata_df, new_readouts_df
            )
            if VERBOSE:
                print(
                    f"   Update will add ~{len(updated_metadata) - len(old_meta_df)} metadata rows "
                    f"and ~{len(updated_readouts) - len(old_exp_df)} readout rows."
                )
                confirm = input("   Proceed? (y/n): ")
                if confirm.lower() != 'y':
                    raise ValueError("Execution cancelled by user.")
            pd.DataFrame(updated_metadata).to_csv(old_meta, index=False)
            pd.DataFrame(updated_readouts).to_csv(old_exp,  index=False)

        if VERBOSE:
            print("-SUCCESS- Record files updated.")

        # ── step 3 (CDD only): filter already-docked compounds ───────────────
        protein_name  = os.path.splitext(os.path.basename(pdb))[0]
        pred_filename = 'predictions.csv' if COVALENT else 'noncov_predictions.csv'
        pred_rec      = os.path.join(record_path, pred_filename)
        error_csv     = os.path.join(record_path, 'errored.csv')

        print(f"3. Checking existing predictions and performing {mode_label} docking with Boltz-2.")

        metadata_df = pd.read_csv(old_meta)
        if os.path.exists(pred_rec):
            pred_df = pd.read_csv(pred_rec)
            dock_compounds = update_predictions.fetch_new(
                pred_df, metadata_df, protein_name, min_replicates=min_replicates
            )
            if os.path.exists(error_csv):
                error_df = pd.read_csv(error_csv)
                dock_compounds = update_predictions.check_attempted(error_df, dock_compounds)
            elif VERBOSE:
                print("   No errored compounds found.")
        else:
            if VERBOSE:
                print(f"   {pred_filename} not found – attempting to dock all compounds.")
            dock_compounds = metadata_df[['substance_id', 'smiles'] + (
                ['inchi_key'] if 'inchi_key' in metadata_df.columns else []
            )]

    else:
        # ── local mode: dock the input directly, touch no record files ────────
        print(f"2. Local mode – skipping record update.")
        print(f"3. Performing {mode_label} docking with Boltz-2.")
        dock_compounds = new_metadata_df

    if VERBOSE:
        print(f"   Docking {len(dock_compounds)} compound(s). Cancel in 5 seconds to abort.")
        time.sleep(5)

    # ── step 4: submit docking jobs ───────────────────────────────────────────
    from BoltzCov.run_boltz import submit_job

    if COVALENT:
        final_status = submit_job.run_boltz_cov(
            prot_file=pdb,
            ligand_df=dock_compounds,
            boltz_cache=boltz_cache,
            run_cache=run_cache,
            res_idx=res_idx,
            ligand_chain=ligand_chain,
            VERBOSE=VERBOSE,
            slurm_template=slurm_template,
            msa_path=msa_path,
        )
    else:
        final_status = submit_job.run_boltz_noncov(
            prot_file=pdb,
            ligand_df=dock_compounds,
            boltz_cache=boltz_cache,
            run_cache=run_cache,
            VERBOSE=VERBOSE,
            slurm_template=slurm_template,
            msa_path=msa_path,
        )

    if final_status == "COMPLETED":
        print("Jobs completed!")
    else:
        print(f"Job failed with status: {final_status}")


# ── entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Boltz-2 covalent/non-covalent docking pipeline."
    )
    parser.add_argument(
        "--json_file",
        required=True,
        help="Path to JSON configuration file.",
    )
    args = parser.parse_args()
    main(args)