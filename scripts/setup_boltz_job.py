import os
import csv
import sys
from pdb_to_fasta import pdb_to_fasta
import argparse

def create_boltz_job(csv_file, pdb_file, output_dir):
    """
    Creates directories with .yaml files based on the input CSV and PDB files.
    :param csv_file: Path to the input CSV file.
    :param pdb_file: Path to the input PDB file.
    :param output_dir: Path to the output directory.
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    # Check if the PDB file exists
    if not os.path.isfile(pdb_file):
        print(f"Error: PDB file '{pdb_file}' does not exist.")
        sys.exit(1)
    # Check if the CSV file exists
    if not os.path.isfile(csv_file):
        print(f"Error: CSV file '{csv_file}' does not exist.")
        sys.exit(1)
    
    # Convert PDB to FASTA
    fasta_dir = "/home/limcaoco/turbo/limcaoco/boltz_benchmark/input_files/fastas/"
    pdb_name = os.path.splitext(os.path.basename(pdb_file))[0]
    fasta_file = os.path.join(fasta_dir, f"{pdb_name}.fasta")
    pdb_to_fasta(pdb_file, fasta_file)

    with open(fasta_file, 'r') as fasta:
        sequence = fasta.read().strip()
    # Read CSV file
    with open(csv_file, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            catalog_id = row["Catalog ID"]
            smiles = row["SMILES"]

            # Create .yaml file in the output directory
            yaml_file = os.path.join(output_dir, f"{catalog_id}.yaml")
            with open(yaml_file, 'w') as yaml:
                yaml.write("version: 1\n")
                yaml.write("sequences:\n")
                yaml.write("  - protein:\n")
                yaml.write("      id: A\n")
                
                yaml.write(f"      sequence: {sequence}\n")
                yaml.write("  - ligand:\n")
                yaml.write(f"      id: B\n")
                yaml.write(f"      smiles: '{smiles}'\n")
                yaml.write("properties:\n")
                yaml.write("  - affinity:\n")
                yaml.write("      binder: B\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Setup Boltz job directories and YAML files.")
    parser.add_argument("--input_csv_file", help="Path to the input CSV file.")
    parser.add_argument("--input_pdb_file", help="Path to the input PDB file.")
    parser.add_argument("--output_directory", help="Path to the output directory.")

    args = parser.parse_args()

    create_boltz_job(args.input_csv_file, args.input_pdb_file, args.output_directory)
