import os
import csv
import pandas as pd 
import argparse
import subprocess
from pdb_to_fasta import residue_to_three_letter, build_sequence
from covalent_utils import verify_covalent, residue_cov_atom, remove_leaving_group
from covalent_utils import process_covalent_smiles

def ensure_environment_variables():
    '''
    Ensures necessary environment variables are set. If not, runs setup_enviorment.sh.
    '''
    setup_script = os.path.join(os.path.dirname(__file__), "../setup_enviorment.sh")
    
    command = f"bash -c 'source {setup_script} && env'"
    proc = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True, executable="/bin/bash")
    output, _ = proc.communicate()
    
    for line in output.decode().splitlines():
        key, _, value = line.partition("=")
        os.environ[key] = value

    print("Environment variables set successfully.")

def validate_file(filename):
    ''' Validates if the file is either a PDB or a TXT file.'''
    valid_exts = {".pdb", ".txt"}
    _, ext = os.path.splitext(filename)
    ext = ext.lower()
    if ext not in valid_exts:
        raise ValueError(f"[ERROR] File '{filename}' must have one of these extensions: {valid_exts}")
    else:
        return ext

def process_protein(pdb, idx):
    '''Returns protein information.'''
    ext = validate_file(pdb)
    if ext==".pdb":
        sequence = build_sequence(pdb)
    else: # txt with sequence
        with open(pdb, 'r') as f:
            content = f.read()
            sequence = "".join(content.split())

    res_aa = sequence[idx]
    res_name = residue_to_three_letter(res_aa)

    if verify_covalent(res_name) != True: # verifies if this residue can participate in a covalent bond w the
        raise ValueError("[ERROR] The res_idx provided does NOT map to a covalent residue. " \
        "Please verify res_idx matches expected residue in sequence.")
    
    res_atom = residue_cov_atom(res_name)
    return sequence, res_name, res_atom

def generate_csv(name, prot_file, res_idx, lig_csv, out_csv, ccd_db):
    '''Generates CSV required for input into setup_cov_job.py with information required by Boltz2 for covalent docking.'''
    
    validate_file(prot_file) # check if either txt or pdb
    ensure_environment_variables()
    ccd_db = os.getenv("CCD_DB")
    
    with open(lig_csv, 'r') as lig:
        reader = csv.reader(lig)
        ligands = [(row[0], row[1]) for row in reader if row]  # gives [("lig1", "SMILES1"), ("lig2", "SMILES2"), ...]

    seq, res_name, res_atom = process_protein(prot_file, res_idx)

    expected_header = ["Compound_ID", "SMILES", "CCD", "WH_Type", "Lig_Atom", "Prot_ID", "Prot_Seq", "Res_Idx", "Res_Name", "Res_Atom"]
    # writing header if needed
    write_header = True
    if os.path.exists(out_csv):
        with open(out_csv, "r") as existing:
            reader = csv.reader(existing)
            first_row = next(reader, None)
            if first_row == expected_header:
                write_header = False  # don’t rewrite header

    with open(out_csv, 'a') as f:
        writer = csv.writer(f)

        if write_header:
            writer.writerow(expected_header)
        # for each ligand, append the protein information, assuming one protein target 
        for lig in ligands:
            id = lig[0]
            smiles_no_lg, lig_atom, wh_type = remove_leaving_group(lig[1])
            ccd = process_covalent_smiles(smiles_no_lg, ccd_db=ccd_db) # makes pkl file 
            writer.writerow([id, smiles_no_lg, ccd, wh_type, lig_atom, str(name), seq, res_name, str(res_idx), res_atom])
        
def main():
    parser = argparse.ArgumentParser(description="Generates CSV required for input into setup_cov_job.py with information required by Boltz2 for covalent docking. " \
                                     "Processes ligands and makes required files for Boltz2." \
                                     "Assumes one protein target. Consult README for input format inforamtion.")
    parser.agg_argument("-n", "--name", type=str, required=True, help="Name of protein. Used for naming output files.")
    parser.add_argument("-p","--prot_file", type=str, required=True, help="Path to either a PDB file or a TXT file with a single chain sequence.")
    parser.add_argument("-r", "--res_idx", type=int, required=True, help="Index of the residue to be covalently targeted by a covalent ligand. Starting at 1.")
    parser.add_argument("-l","--lig_csv", type=str, required=True, help="Path to CSV with Ligand info.")
    parser.add_argument("-o","--out_csv", type=str, required=True, help="Path to output CSV. Will be formatted to work with setup_cov_job.py.")
    parser.add_argument("-c","--ccd_db", type=str, required=True, help="Path to directory with covalent compound pkl files", default="/home/ymanasa/.boltz/mols")

    args = parser.parse_args()

    generate_csv(args.name, args.prot_file, args.res_idx, args.lig_csv, args.out_csv, args.ccd_db)

if __name__ == "__main__":
    main()
    