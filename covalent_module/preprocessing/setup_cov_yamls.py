import os
import sys
import yaml 
import pandas as pd
import subprocess
import argparse
from rdkit import Chem

# use ccd_pkl env

def check_smiles(smiles: str):
    '''
    Attempts to load and sanitize a SMILES string using RDKit.
    Returns a canonicalized SMILES string if successful, otherwise None.
    :param smiles: str 
        The input SMILES string.

    :return: str or None
        A valid, canonical SMILES or None if the molecule is invalid.
    '''
    try:
        # Attempt to parse without sanitizing
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            print(f"[ERROR] MolFromSmiles failed for: {smiles}")
            return None
        
        # Attempt sanitization (includes valence check, aromaticity, Hs)
        Chem.SanitizeMol(mol)

        # Return the canonical SMILES
        return Chem.MolToSmiles(mol, canonical=True)

    except Exception as e:
        print(f"[ERROR] Sanitization failed for SMILES: {smiles}\n{e}")
        return None

class LiteralList(list):
        pass
    
def literal_list_representer(dumper, data):
    return dumper.represent_sequence("tag:yaml.org,2002:seq", data, flow_style=True)
yaml.add_representer(LiteralList, literal_list_representer)

def create_boltz_yamls(csv_file, output_dir, msa_path):
    '''
    Creates YAML files from a CSV of ligands and proteins.

    :param csv_file: Path to input CSV file
    :param output_dir: Directory to write YAML files
    :param msa_path: Optional path to MSA file

    :return: List of paths to created YAML files
    '''
    VERBOSE = os.environ.get("VERBOSE", "FALSE").upper() == "TRUE"
    CCD_DB = os.environ.get("CCD_DB", output_dir)

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    if not csv_file or not os.path.isfile(csv_file):
        print(f"[ERROR] CSV file '{csv_file}' does not exist.")
        sys.exit(1)
    
    # load csv and check columns
    csvfile = pd.read_csv(csv_file)
    required_columns = {"SMILES", "compound_id", "vault_id" ,"WH_Type", "Lig_Atom", "Prot_ID", "Prot_Seq", "Res_Idx", "Res_Name", "Res_Atom"}
    missing = required_columns - set(csvfile.columns)
    if missing:
        print(f"[ERROR] CSV file is missing these columns: {missing}. \
              Use make_input_csv.py to generate the correct format.")
        sys.exit(1)

    invalid_compounds = []
    yaml_files = []
    for _, row in csvfile.iterrows(): # per ligand yaml is made 
        # ligand info 
        ccd = row["compound_id"]
        smiles = row["SMILES"]
        smiles = check_smiles(smiles) # returns conancial smiles or None
        if smiles is None:
            print(f"[ERROR] Invalid SMILES for compound ID {ccd}: {row['SMILES']}")
            invalid_compounds.append(ccd)
            continue
        # check if ccd pkl file exists
        ccd_file = os.path.join(CCD_DB, f"{ccd}.pkl")
        if not os.path.isfile(ccd_file):
            print(f"[ERROR] '{ccd_file}' does not exist for {ccd}. \
                  Please use preprocessing script (/preprocessing/make_input_csv.py) to generate it.")
            invalid_compounds.append(ccd)
            continue
        lig_atom = row["Lig_Atom"] 

        # protein info 
        pdb_name = row["Prot_ID"]
        sequence = row["Prot_Seq"]
        res_idx = row["Res_Idx"]
        res_name = row["Res_Name"]
        prot_atom = row["Res_Atom"]

        protein_data = {
            "id": "A",
            "sequence": sequence,
            "modification": [{"position": int(res_idx), "ccd": res_name}],
        }

        if msa_path is not None:
            protein_data["msa"] = msa_path

        data = {
            "sequences": [
                {"protein": protein_data},
                {"ligand": {"id": "LIG", "ccd": str(ccd)}},
            ],
            "constraints": [
                {
                    "bond": {
                        "atom1": ["A", int(res_idx), prot_atom],
                        "atom2": ["LIG", int(1), lig_atom],
                    }
                }
            ],
            "properties": [
                {"affinity": {"binder": "LIG"}}
            ],
        }
        yaml_file = os.path.join(output_dir, f"{pdb_name}_{ccd}.yaml") # should be unique for each ligand and less than 5 char
        with open(yaml_file, "w") as f:
            yaml.safe_dump(
                data, 
                f,
                sort_keys=False,
                indent=4,
                width=4096,  # prevents wrapping long strings
                default_flow_style=False
            )

        yaml_files.append(yaml_file)
    
    if invalid_compounds:
        print(f"[WARNING] The following compound IDs were skipped: {invalid_compounds}")

    if VERBOSE: print("[DONE] yamls were written in", output_dir)
    return yaml_files # list of all yamls created 