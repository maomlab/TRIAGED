import os
import sys
import yaml 
import pandas as pd
import subprocess
import argparse
from rdkit import Chem

# tested last: 10/03/25 in test_yamls/
# use ccd_pkl env

def ensure_environment_variables():
    '''
    Ensures necessary environment variables are set. If not, runs setup_enviorment.sh.
    '''
    setup_script = os.path.join(os.path.dirname(__file__), "/home/ymanasa/turbo/ymanasa/opt/maom_boltz/covalent_scripts/setup_enviorment.sh")
    
    command = f"bash -c 'source {setup_script} && env'"
    proc = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True, executable="/bin/bash")
    output, _ = proc.communicate()
    
    for line in output.decode().splitlines():
        key, _, value = line.partition("=")
        os.environ[key] = value

    print("Environment variables set successfully.")

def check_smiles(smiles: str, verbose: bool = False):
    '''
    Attempts to load and sanitize a SMILES string using RDKit.
    Returns a canonicalized SMILES string if successful, otherwise None.
    
    :param smiles: str 
        The input SMILES string.
    :param verbose: bool
        If True, print debug messages on failure.

    :return: str or None
        A valid, canonical SMILES or None if the molecule is invalid.
    '''
    try:
        # Attempt to parse without sanitizing
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            if verbose:
                print(f"[ERROR] MolFromSmiles failed for: {smiles}")
            return None
        
        # Attempt sanitization (includes valence check, aromaticity, Hs)
        Chem.SanitizeMol(mol)

        # Return the canonical SMILES
        return Chem.MolToSmiles(mol, canonical=True)

    except Exception as e:
        if verbose:
            print(f"[ERROR] Sanitization failed for SMILES: {smiles}\n{e}")
        return None

class LiteralList(list):
        pass
def literal_list_representer(dumper, data):
    return dumper.represent_sequence("tag:yaml.org,2002:seq", data, flow_style=True)
yaml.add_representer(LiteralList, literal_list_representer)

def create_boltz_yamls(csv_file, output_dir):
    '''
    Creates .yaml files based on the input CSV.
    :param csv_file: str
        Path to the input CSV file. 
    :param output_dir: str
        Path to the output directory.

    :return: path to last yaml file created
    '''
    # Ensure environment variables are set
    ensure_environment_variables()
    ccd_db = os.getenv("CCD_DB", "/home/$USER/.bolts/mols") 

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    if not csv_file or not os.path.isfile(csv_file):
        print(f"[ERROR] CSV file '{csv_file}' does not exist.")
        sys.exit(1)
    
    # load csv and check columns
    csvfile = pd.read_csv(csv_file)
    required_columns = {"Compound_ID", "SMILES", "CCD", "WH_Type", "Lig_Atom", "Prot_ID", "Prot_Seq", "Res_Idx", "Res_Name", "Res_Atom"}
    missing = required_columns - set(csvfile.columns)
    if missing:
        print(f"[ERROR] CSV file is missing these columns: {missing}")
        sys.exit(1)

    invalid_compounds = []
    yaml_files = []
    for _, row in csvfile.iterrows(): # per ligand yaml is made 
        # ligand info 
        ccd = row["CCD"]
        smiles = row["SMILES"]
        smiles = check_smiles(smiles, verbose=True) # returns conancial smiles or None
        if smiles is None:
            print(f"Invalid SMILES for compound ID {ccd}: {row['SMILES']}")
            invalid_compounds.append(ccd)
            continue
        # check if ccd mol file exists
        ccd_file = os.path.join(ccd_db, f"{ccd}.pkl")
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

        data = {
            "sequences": [
                {
                    "protein": {
                        "id": "A",
                        "sequence": sequence,
                        "modification": [
                            {"position": int(res_idx), "ccd": res_name}
                        ],
                    }
                },
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
        yaml_file = os.path.join(output_dir, f"{pdb_name}_{ccd}.yaml") # should be unique for each ligand
        with open(yaml_file, "w") as f:
            yaml.safe_dump(
                data, f,
                sort_keys=False,
                indent=4,
                width=4096,  # prevents wrapping long strings
                default_flow_style=False
            )

        yaml_files.append(yaml_file)
    
    if invalid_compounds:
        print(f"[WARNING] The following compound IDs were skipped: {invalid_compounds}")

    print("[DONE] yamls were written in", output_dir)

    return yaml_files # list of all yamls created 

def main():
    parser = argparse.ArgumentParser(description="Creates YAML files for covalent docking with Boltz2 in output_dir/. Refer to README.md for csv format.")
    parser.add_argument("-i","--input_csv_file", type=str, required=True,help="Path to the input CSV file.")
    parser.add_argument("-o","--output_directory", type=str, required=True,help="Path to the output directory for generated yamls.")

    args = parser.parse_args()

    create_boltz_yamls(csv_file=args.input_csv_file, output_dir=args.output_directory)

if __name__ == "__main__":
    main()
    