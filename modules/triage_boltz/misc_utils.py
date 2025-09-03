import os
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
import csv
from typing import Union
from mmseqs2 import run_mmseqs2
import subprocess
def check_smiles(smiles: str, verbose: bool = True) -> Union[str, None]:
    """
    Attempts to load and sanitize a SMILES string using RDKit.
    Returns a canonicalized SMILES string if successful, otherwise None.
    
    Parameters:
    - smiles (str): The input SMILES string.
    - verbose (bool): If True, print debug messages on failure.

    Returns:
    - str or None: A valid, canonical SMILES or None if the molecule is invalid.
    """
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
    
def ensure_environment_variables():
    """
    Ensures necessary environment variables are set. If not, runs setup_enviorment.sh.
    """
    if not os.getenv("PROJECT_DIR") or os.getenv(""):
        print("Environment variable PROJECT_DIR is not set. Running setup_enviorment.sh...")

        setup_script = os.path.join(os.path.dirname(__file__), "setup_enviorment.sh")
        subprocess.run(f"source {setup_script}", shell=True, executable="/bin/bash", check=True)
        print("Environment variables set successfully.")
    else:
        print("Environment variables are already set. Continuing...")

def sanitize_compound_id(compound_ID: str) -> str:
    """
    Adjust the compound_ID to make it file-system friendly.
    Parameters:
    - compound_ID (str): The compound ID to sanitize.
    Returns:
    - str: The sanitized compound ID.
    """
    if any(char in compound_ID for char in [' ', '"', "'", ",", "(", ")", "/", "\\"]):
        sanitized_compound_id = (
            compound_ID.replace(" ", "_")  # Replace spaces with underscores
                       .replace('"', '')    # Remove double quotes
                       .replace("'", '')    # Remove single quotes
                       .replace(",", '')    # Remove commas
                       .replace("(", '')    # Remove parentheses
                       .replace(")", '')    # Remove parentheses
                       .replace("/", '_')   # Replace slashes with underscores
                       .replace("\\", '_')  # Replace backslashes with underscores
        )
        print(f"Found invalid characters in compound ID '{compound_ID}'. Sanitized to '{sanitized_compound_id}'.")
    else:
        sanitized_compound_id = compound_ID
    return sanitized_compound_id


def parse_input_csv(csv_file: str) -> list[dict]:
    """
    Parse a CSV file with varying numbers of compound_ID, SMILES, and optional num columns.
    Also supports optional cofactor_ID_X, cofactor_sequence_X, and cofactor_num_X columns.
    If num_X or cofactor_num_X is missing or empty, it will be set to 1.

    Parameters:
    - csv_file: Path to the CSV file.

    Returns:
    - A list of dictionaries, where each dictionary contains compound_ID, SMILES, num, and optional cofactor data.
    """
    parsed_data = []


    with open(csv_file, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            compound_data = []
            #print(f"Processing row: {row}")
            #print(f"Row keys: {list(row.keys())}")
            for key in row.keys():
                if key.startswith("compound_ID") or key.startswith("protein_ID") or key.startswith("dna_ID") or key.startswith("rna_ID") and row[key]:
                    if (len(key) > len("compound_ID")) or (len(key) > len("protein_ID")) or (len(key) > len("dna_ID")) or (len(key) > len("rna_ID")):
                        if key.startswith("compound_ID"):
                            suffix = key[len("compound_ID"):]
                        elif key.startswith("dna_ID"):
                            suffix = key[len("dna_ID"):]
                        elif key.startswith("rna_ID"):
                            suffix = key[len("rna_ID"):]
                        elif key.startswith("protein_ID"):
                            suffix = key[len("protein_ID"):]
                        smiles_key = f"SMILES{suffix}"
                        num_key = f"compound_num{suffix}"
                        inchi_key = f"InChI{suffix}"
                        protein_id_key = f"protein_ID{suffix}"
                        protein_sequence_key = f"protein_sequence{suffix}"
                        protein_num_key = f"protein_num{suffix}"
                        pdb_path_key = f"pdb_path{suffix}"
                        pdb_id_key = f"pdb_id{suffix}"
                        dna_sequence_key = f"dna_sequence{suffix}"
                        dna_num_key = f"dna_num{suffix}"
                        rna_sequence_key = f"rna_sequence{suffix}"
                        rna_num_key = f"rna_num{suffix}"
                        constraints_key = f"pocket_constraints{suffix}"
                    else:
                        smiles_key = "SMILES"
                        num_key = "compound_num"
                        inchi_key = "InChI"
                        protein_id_key = "protein_ID"
                        protein_sequence_key = "protein_sequence"
                        protein_num_key = "protein_num"
                        pdb_path_key = "pdb_path"
                        pdb_id_key = "pdb_id"
                        dna_sequence_key = "dna_sequence"
                        dna_num_key = "dna_num"
                        rna_sequence_key = "rna_sequence"
                        rna_num_key = "rna_num"
                        constraints_key = "pocket_constraints"

                    compound_entry = None
                    if key.startswith("compound_ID") and row[key]:
                        if smiles_key in row and row[smiles_key]:
                            compound_entry = {
                                "compound_ID": row[key].strip(),
                                "SMILES": row[smiles_key].strip(),
                                "compound_num": row[num_key].strip() if num_key in row and row[num_key] else "1",
                                "inchi": row[inchi_key].strip() if inchi_key in row and row[inchi_key] else None
                            }
                    if key.startswith("protein_ID") and row[key]:
                        if protein_id_key in row and row[protein_id_key]:
                            compound_entry = {
                                "protein_ID": row[protein_id_key].strip(),
                                "protein_sequence": row[protein_sequence_key].strip() if protein_sequence_key in row and row[protein_sequence_key] else None,
                                "protein_num": row[protein_num_key].strip() if protein_num_key in row and row[protein_num_key] else "1",
                                "pdb_path": row[pdb_path_key].strip() if pdb_path_key in row and row[pdb_path_key] else None,
                                "pdb_id": row[pdb_id_key].strip() if pdb_id_key in row and row[pdb_id_key] else None
                            }
                    if key.startswith("dna_ID") and row[key]:
                        if dna_sequence_key in row and row[dna_sequence_key]:
                            compound_entry = {
                                "dna_ID": row[key].strip(),
                                "dna_sequence": row[dna_sequence_key].strip(),
                                "dna_num": row[dna_num_key].strip() if dna_num_key in row and row[dna_num_key] else "1"
                            }
                    if key.startswith("rna_ID") and row[key]:
                        if rna_sequence_key in row and row[rna_sequence_key]:
                            compound_entry = {
                                "rna_ID": row[key].strip(),
                                "rna_sequence": row[rna_sequence_key].strip(),
                                "rna_num": row[rna_num_key].strip() if rna_num_key in row and row[rna_num_key] else "1"
                            }
                    # Parse constraints if present
                    if compound_entry is not None and constraints_key in row and row[constraints_key]:
                        # Expecting constraints as a string, e.g. "A:100/CA;B:50/N"
                        constraints_str = row[constraints_key].strip()
                        constraints = {}
                        for item in constraints_str.split(';'):
                            item = item.strip()
                            if not item:
                                continue
                            # Format: chain:res_idx/atom_name, e.g. A:100/CA
                            try:
                                chain_part, res_atom = item.split(':')
                                res_idx, atom_name = res_atom.split('/')
                                if chain_part not in constraints:
                                    constraints[chain_part] = []
                                constraints[chain_part].append((int(res_idx), atom_name))
                            except Exception as e:
                                print(f"[WARNING] Could not parse constraint '{item}': {e}")
                        compound_entry["constraints"] = constraints
                    if compound_entry is not None:
                        compound_data.append(compound_entry)
            parsed_data.append(compound_data)
    return parsed_data


def generate_msa(entity_name: str, sequence: Union[str, list[str]]) -> str:
    """
    Generate a Multiple Sequence Alignment (MSA) using mmseqs2.

    Parameters:
]    - entity_name (str): name of the receptor (without extension).
    - sequence (str): The sequence to use for MSA generation.

    Returns:
    - str: Path to the generated MSA file.
    """
    project_dir = os.getenv("PROJECT_DIR")

    msa_dir = os.path.join(project_dir, "input_files/msa/")
    if isinstance(sequence, str):
        msa_file = os.path.join(msa_dir, f"{entity_name}_mmseqs2.a3m")
    elif isinstance(sequence, list):
        msa_file = os.path.join(msa_dir, f"{entity_name}_mmseqs2_pair.a3m")
    # Check if MSA file already exists
    if os.path.exists(msa_file):
        print(f"MSA file already exists at {msa_file}. Skipping MSA generation.")
        return msa_file

    print("Now pregenerating MSA with mmseqs2...")
    os.makedirs(msa_dir, exist_ok=True)

    # Run mmseqs2 to generate MSA
    #print(f"Running mmseqs2 for entity '{entity_name}' with sequence:")
    #print(sequence)
    if isinstance(sequence, list):
        msa_result = run_mmseqs2(sequence, prefix=f"{msa_dir}/{entity_name}", use_pairing=True)
    else:
        msa_result = run_mmseqs2(sequence, prefix=f"{msa_dir}/{entity_name}", use_pairing=False)
    #print(msa_result)
    #print(msa_result.keys())
    with open(msa_file, 'w') as msa_output:
        msa_output.write("\n".join(msa_result))
    print(f"MSA saved to {msa_file}")
    print(f"MSA generation completed for entity '{entity_name}'.")
    # Verify that the MSA file was created successfully
    if not os.path.isfile(msa_file):
        print(f"Error: MSA file '{msa_file}' was not created successfully.")
        sys.exit(1)
    msa_file_path = msa_file
    return msa_file_path

def clean_dos_chars(input_file: str, output_file: str = None) -> dict:
    """
    Detect and clean DOS/Windows-specific characters from a text file.

    Parameters
    ----------
    input_file : str
        Path to the input file (.csv or .txt).
    output_file : str, optional
        Path to the cleaned output file. If None, overwrites input_file.

    Returns
    -------
    dict
        A dictionary summarizing findings:
        {
            "had_bom": bool,
            "had_crlf": bool,
            "had_cr": bool,
            "output_file": str
        }
    """
    with open(input_file, "rb") as f:
        raw = f.read()

    findings = {
        "had_bom": raw.startswith(b"\xef\xbb\xbf"),
        "had_crlf": b"\r\n" in raw,
        "had_cr": b"\r" in raw and b"\r\n" not in raw,  # stray Mac-style CR
        "output_file": None
    }

    # Normalize: remove BOM, convert CRLF -> LF, CR -> LF
    text = raw.decode("utf-8-sig", errors="replace")  # removes BOM if present
    text = text.replace("\r\n", "\n")  # DOS to Unix
    text = text.replace("\r", "\n")   # old Mac to Unix

    if output_file is None:
        output_file = input_file

    with open(output_file, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)

    findings["output_file"] = os.path.abspath(output_file)
    return findings