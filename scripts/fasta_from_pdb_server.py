import argparse
import requests
from Bio import SeqIO
from io import StringIO
import os
import pandas as pd

def fetch_fasta(pdb_id):
    """Fetch the FASTA sequences for all chains of a given PDB ID"""
    url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
    response = requests.get(url)
    if response.status_code != 200:
        print(f"Failed to fetch FASTA for {pdb_id}")
    return response.text

def fetch_smiles(pdb_id):
    """Fetch SMILES strings of all ligands in the given PDB entry"""
    ligands_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    entry_data = requests.get(ligands_url).json()
    
    smiles_dict = {}
    try:
        for ligand in entry_data.get("rcsb_entry_container_identifiers", {}).get("nonpolymer_entity_ids", []):
            ligand_url = f"https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{pdb_id}/{ligand}"
            ligand_data = requests.get(ligand_url).json()
            chem_id = ligand_data['chem_comp']['id']
            smiles = ligand_data['chem_comp'].get('rcsb_chem_comp_descriptor', {}).get('smiles', None)
            if smiles:
                smiles_dict[chem_id] = smiles
    except Exception as e:
        print(f"Error retrieving ligands for {pdb_id}: {e}")
    return smiles_dict

def save_fasta_to_file(fasta_str, out_path):
    with open(out_path, "w") as f:
        f.write(fasta_str)

def save_smiles_to_file(smiles_dict, out_path):
    with open(out_path, "w") as f:
        for ligand_id, smiles in smiles_dict.items():
            f.write(f"{smiles}\t{ligand_id}\n")

def main():
    parser = argparse.ArgumentParser(description="Fetch FASTA and SMILES from the PDB.")
    parser.add_argument("--pdb_id", help="PDB ID (e.g., 1A2B)", default=None)
    parser.add_argument("--receptor-name", help="Name of the receptor (optional, for output naming)", default=None)
    parser.add_argument("--fasta-out", help="Output path for FASTA file", default=None)
    parser.add_argument("--smiles-out", help="Output path for SMILES file", default=None)
    parser.add_argument("--input-csv", help="Input CSV file with Receptor_IDs, pdb_IDs (optional, for batch processing)", default=None)
    args = parser.parse_args()

    if args.input_csv:
        input_df = pd.read_csv(args.input_csv)
        for index, row in input_df.iterrows():
            pdb_id = row['pdb_ID'].upper()
            receptor_name = row['Receptor_ID'].strip() if 'Receptor_ID' in row else None
            print(f"Processing PDB ID: {pdb_id}, Receptor Name: {receptor_name}")
        
            if args.fasta_out is None:
                print("No output path for FASTA provided, skipping FASTA fetch.")
            else:
                print(f"Fetching FASTA for {pdb_id}...")
                fasta = fetch_fasta(pdb_id)
                if receptor_name:
                    fasta_out_path = os.path.join(args.fasta_out, f"{receptor_name}.fasta")
                else:
                    fasta_out_path = os.path.join(args.fasta_out, f"{pdb_id}.fasta")
                save_fasta_to_file(fasta, fasta_out_path)
                print(f"Saved FASTA to {fasta_out_path}")

            if args.smiles_out is None:
                print("No output path for SMILES provided, skipping SMILES fetch.")
            else:
                print(f"Fetching ligand SMILES for {pdb_id}...")
                smiles = fetch_smiles(pdb_id)
                smiles_out_path = os.join.path(args.smiles_out, f"{pdb_id}_ligands.smi")
                save_smiles_to_file(smiles, smiles_out_path)
                print(f"Saved SMILES to {smiles_out_path}")
    else:
        pdb_id = args.pdb_id.upper()
        receptor_name = args.receptor_name or pdb_id
        
        if args.fasta_out:
            print(f"Fetching FASTA for {pdb_id}...")
            fasta = fetch_fasta(pdb_id)
            fasta_out_path = os.path.join(args.fasta_out, f"{receptor_name}.fasta")
            save_fasta_to_file(fasta, fasta_out_path)
            print(f"Saved FASTA to {fasta_out_path}")

        if args.smiles_out:
            print(f"Fetching ligand SMILES for {pdb_id}...")
            smiles = fetch_smiles(pdb_id)
            smiles_out_path = os.path.join(args.smiles_out, f"{pdb_id}_ligands.smi")
            save_smiles_to_file(smiles, smiles_out_path)
            print(f"Saved SMILES to {smiles_out_path}")
if __name__ == "__main__":
    main()
