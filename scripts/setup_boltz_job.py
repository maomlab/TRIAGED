import os
import csv
import sys
import subprocess
# from pdb_to_fasta import pdb_to_fasta
from fasta_utils import build_fasta_seq
from mmseqs2 import run_mmseqs2
import argparse
from covalent_utils import get_link_atoms

from rdkit import Chem
from rdkit.Chem import AllChem

def check_smiles(smiles: str, verbose: bool = False):
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

        # Optionally do other cleanups (like kekulization or 2D coordinates)
        # AllChem.Compute2DCoords(mol)

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
    if not os.getenv("PROJECT_DIR"):
        # print("Environment variable PROJECT_DIR is not set. Running setup_enviorment.sh...")
        setup_script = os.path.join(os.path.dirname(__file__), "setup_enviorment.sh")
        subprocess.run(f"source {setup_script}", shell=True, executable="/bin/bash", check=True)
        print("Environment variables set successfully.")


def create_boltz_job(pdb_file, output_dir, csv_file=None, covalent_docking=False):
    """
    Creates directories with .yaml files based on the input CSV and PDB files.
    :param csv_file: Path to the input CSV file. Not required for covalent docking. 
    :param pdb_file: Path to the input PDB file.
    :param output_dir: Path to the output directory.
    :param covalent_docking: Whether to write yaml for covalent docking. 
    """
    # Ensure environment variables are set
    ensure_environment_variables()

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Check if the PDB file exists
    if not os.path.isfile(pdb_file):
        print(f"Error: PDB file '{pdb_file}' does not exist.")
        sys.exit(1)
    # Check if the CSV file exists only for non-covalent docking 
    if not covalent_docking:
        if not csv_file or not os.path.isfile(csv_file):
            print(f"Error: CSV file '{csv_file}' does not exist.")
            sys.exit(1)
    
    # Convert PDB to FASTA
    project_dir = os.getenv("PROJECT_DIR", "/home/ymanasa/turbo/ymanasa/boltz_benchmark")
    fasta_dir = os.path.join(project_dir, output_dir)
    pdb_name = os.path.splitext(os.path.basename(pdb_file))[0]
    fasta_file = os.path.join(fasta_dir, f"{pdb_name}.fasta")
    # pdb_to_fasta(pdb_file, fasta_file)
    sequence = build_fasta_seq(pdb_name)
    # Read the FASTA sequence
    # with open(fasta_file, 'r') as fasta:
    #     sequence = fasta.read().strip()

    if not covalent_docking:

        # Generate MSA using mmseqs2
        if os.path.exists(os.path.join(project_dir, "input_files/msa", f"{pdb_name}_mmseqs2.a3m")):
            print(f"MSA file already exists at {os.path.join(project_dir, 'input_files/msa', f'{pdb_name}_mmseqs2.a3m')}. Skipping MSA generation.")

        print("now pregenerating MSA with mmseqs2...")
        msa_dir = os.path.join(project_dir, "input_files/msa/")
        os.makedirs(msa_dir, exist_ok=True)
        msa_file = os.path.join(msa_dir, f"{pdb_name}_mmseqs2.a3m")
        msa_result = run_mmseqs2(sequence, prefix=f"{msa_dir}/{pdb_name}")
        with open(msa_file, 'w') as msa_output:
            msa_output.write("\n".join(msa_result))
        print(f"MSA saved to {msa_file}")
        if not os.path.isfile(msa_file):
            print(f"Error: MSA file '{msa_file}' was not created successfully.")
            sys.exit(1)

        # Read CSV file
        with open(csv_file, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                catalog_id = row["compound_ID"]
                smiles = row["SMILES"]
                smiles = check_smiles(smiles, verbose=True)
                if smiles is None:
                    print(f"Invalid SMILES for compound ID {catalog_id}: {row['SMILES']}")
                    continue
                # Create .yaml file in the output directory
                yaml_file = os.path.join(output_dir, f"{catalog_id}.yaml")
                with open(yaml_file, 'w') as yaml:
                    yaml.write("version: 1\n")
                    yaml.write("sequences:\n")
                    yaml.write("  - protein:\n")
                    yaml.write("      id: A\n")
                    yaml.write(f"      sequence: {sequence}\n")
                    yaml.write(f"      msa: {msa_file}\n")

                    yaml.write("  - ligand:\n")
                    yaml.write(f"      id: B\n")
                    yaml.write(f"      smiles: '{smiles}'\n")
                    yaml.write("properties:\n")
                    yaml.write("  - affinity:\n")
                    yaml.write("      binder: B\n")
    else:
       
        prot_atom, res_name, res_idx, _, ccd, lig_atom, _ = get_link_atoms(pdb_file)

        yaml_file = os.path.join(output_dir, f"{pdb_name}_{ccd}.yaml")
        with open(yaml_file, 'w') as yaml:
            yaml.write("version: 1\n")
            yaml.write("sequences:\n")
            yaml.write("  - protein:\n")
            yaml.write("      id: A\n")
            yaml.write(f"      sequence: {sequence}\n")
            # yaml.write(f"      msa: {msa_file}\n")
            yaml.write(f"      modification:\n")
            yaml.write(f"        - position: {res_idx}\n") # fix res idx 
            yaml.write(f"          ccd: {res_name}\n")

            yaml.write("  - ligand:\n")
            yaml.write(f"      id: LIG\n")
            yaml.write(f"      ccd: '{ccd}'\n")
            
            yaml.write("constraints:\n")
            yaml.write("    - bond:\n")
            yaml.write(f"        atom1: [A, {res_idx}, {prot_atom}]\n") # fix res idx 
            yaml.write(f"        atom2: [LIG, 1, {lig_atom}]\n") # lig atm idx needs to be 1 regardless of PDB idx 

            yaml.write("properties:\n")
            yaml.write("    - affinity:\n")
            yaml.write(f"        binder: LIG\n")
    print(f"YAML files created successfully: {yaml_file}.")

def main():
    parser = argparse.ArgumentParser(description="Setup Boltz job directories and YAML files.")
    parser.add_argument("-i","--input_csv_file", type=str, required=False,help="Path to the input CSV file. Required for non-covalent docking.")
    parser.add_argument("-p","--input_pdb_file", type=str, required=True,help="Path to the input PDB file.")
    parser.add_argument("-o","--output_directory", type=str, required=True,help="Path to the output directory.")
    parser.add_argument("-c", "--covalent_docking", action='store_true', default=False,help="Whether ligand must covlanetly interact with protein")

    args = parser.parse_args()
    if not args.covalent_docking and args.input_csv_file is None:
        parser.error("--input_csv_file is required when --covalent_docking is False for non-covalent docking")

    create_boltz_job(args.input_csv_file, args.input_pdb_file, args.output_directory, args.covalent_docking)

if __name__ == "__main__":
    main()
    