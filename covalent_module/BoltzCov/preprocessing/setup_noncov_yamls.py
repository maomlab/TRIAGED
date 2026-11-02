import os
import sys
import yaml 
from rdkit import Chem
from BoltzCov.preprocessing import make_csv_for_yaml

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

def create_boltz_yamls(prot_file, ligand_df, output_dir, msa_path=None):
    '''
    Creates YAML files from a CSV of ligands and proteins.

    :param csv_file: Path to input CSV file
    :param output_dir: Directory to write YAML files
    :param msa_path: Optional path to MSA file

    :return: List of paths to created YAML files
    '''
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    #protein seq
    sequence, _, _ = make_csv_for_yaml.process_protein(prot_file, idx=None, lig_chain='A')
    if len(sequence) < 1:
        raise ValueError('sequence did not build from given pdb')
    
    invalid_compounds = []
    yaml_files = []
    for _, row in ligand_df.iterrows():  # per ligand yaml is made 
        smiles = row['smiles']

        smiles = check_smiles(smiles) # returns conancial smiles or None
        if smiles is None:
            print(f"[ERROR] Invalid SMILES for compound {row['substance_id']}: {row['SMILES']}")
            invalid_compounds.append({row['substance_id']})
            continue

        protein_data = {
            "id": "A",
            "sequence": sequence,
        }

        if msa_path is not None:
            protein_data["msa"] = msa_path

        data = {
            "sequences": [
                {"protein": protein_data},
                {"ligand": {"id": "LIG", "smiles": smiles}},
            ],
            "properties": [
                {"affinity": {"binder": "LIG"}}
            ],
        }
        yaml_file = os.path.join(output_dir, f"{row['substance_id']}.yaml") # should be unique for each ligand 
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
        print(f"[WARNING] The following compounds were skipped: {invalid_compounds}")

    return yaml_files # list of all yamls created 