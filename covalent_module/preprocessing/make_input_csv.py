import os
import csv
import pandas as pd 
import tempfile
from dotenv import load_dotenv
from .pdb_to_fasta import residue_to_three_letter, build_sequence
from .covalent_utils import verify_covalent, residue_cov_atom, remove_leaving_group, process_covalent_smiles, lookup_compound_id, unique_ccd

# use ccd_pkl env
def validate_file(filename):
    ''' Validates if the file is either a PDB or a TXT file.'''
    valid_exts = {".pdb", ".txt"}
    _, ext = os.path.splitext(filename)
    ext = ext.lower()
    if ext not in valid_exts:
        raise ValueError(f"[ERROR] File '{filename}' must have one of these extensions: {valid_exts}")
    else:
        return ext

def process_protein(pdb, idx, lig_chain):
    '''Returns protein information.'''
    VERBOSE = os.environ.get("VERBOSE", "FALSE").upper() == "TRUE"
    ext = validate_file(pdb)
    if ext==".pdb":
        sequence = build_sequence(pdb, lig_chain)
    else: # txt with sequence
        with open(pdb, 'r') as f:
            content = f.read()
            sequence = "".join(content.split())
    if idx < 1:
        idx = 0
    elif idx > len(sequence):
        raise ValueError(f"[ERROR] res_idx {idx} exceeds sequence length {len(sequence)}.")
    else:
        res_aa = sequence[idx-1]

    res_name = residue_to_three_letter(res_aa)
    if VERBOSE: print("Boltz will now dock to this residue: ", res_name)

    if verify_covalent(res_name) != True: # verifies if this residue can participate in a covalent bond w the
        if VERBOSE: print(sequence)
        raise ValueError(f"[ERROR] res_idx {idx} does NOT map to a covalent residue. " \
        "Please verify res_idx matches expected residue in sequence.")
    
    res_atom = residue_cov_atom(res_name)
    return sequence, res_name, res_atom

def generate_csv(name, prot_file, res_idx, lig_chain, out_csv, ccd_db):
    '''Generates CSV required for input into setup_cov_job.py with information required by Boltz2 for covalent docking.'''
    VERBOSE = os.environ.get("VERBOSE", "FALSE").upper() == "TRUE"
    COMPOUND_RECORD = os.environ.get("COMPOUND_RECORD", None) # none if user does not want to pass ligand id records 
    LIG_CSV = os.environ.get("LIGAND_CSV")

    ## protein processing
    seq, res_name, res_atom = process_protein(prot_file, res_idx, lig_chain)
    
    ## ligand processing
    with open(LIG_CSV, 'r') as lig:
        reader = csv.reader(lig)
        header = next(reader)  
        name_idx = header.index('vault_id')
        smiles_idx = header.index('SMILES')
        
        ligands = [(row[name_idx], row[smiles_idx]) for row in reader if row]

    # initiate compound_records df
    if COMPOUND_RECORD is None:
        if VERBOSE: print("Compound Records not Found! Writing new one...")
        compound_rec_df = pd.DataFrame(columns=['vault_id', 'compound_id'])
        new_rows_list = [{'vault_id': name, 'compound_id': 'XXXXXXX'} for name, _ in ligands]
        compound_rec_copy = pd.concat([compound_rec_df, pd.DataFrame(new_rows_list)], ignore_index=True)
    else:
        # check cols in cmpd_rec and reads
        cmp_rec = pd.read_csv(COMPOUND_RECORD, nrows=0).columns
        missing_cols = [c for c in ['vault_id', 'compound_id'] if c not in cmp_rec]
        if missing_cols:
            raise ValueError(f"COMPOUND_RECORD is missing expected columns: {missing_cols}")
        compound_rec_df = pd.read_csv(COMPOUND_RECORD)
        compound_rec_copy = compound_rec_df.copy(deep=True)

    ## writing output csv for yaml making
    if os.path.exists(out_csv):
        if VERBOSE: print(f"[WARNING] Output CSV '{out_csv}' already exists. Deleting and rewriting.")
        os.remove(out_csv)

    expected_header = ["SMILES", "compound_id", "vault_id" ,"WH_Type", "Lig_Atom", "Prot_ID", "Prot_Seq", "Res_Idx", "Res_Name", "Res_Atom"]
    # write header once
    write_header = True
    if os.path.exists(out_csv):
        with open(out_csv, "r") as existing:
            reader = csv.reader(existing)
            first_row = next(reader, None)
            if first_row == expected_header:
                write_header = False  # won't rewrite header

    with open(out_csv, 'a') as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(expected_header)
        # for each ligand, append the protein information, assuming one protein target 
        for lig in ligands:
            vault_id = lig[0]
            # matches vault_id to unqiue compound_id
            compound_id = lookup_compound_id(vault_id, compound_rec_copy) 
            if compound_id is None or len(compound_id) > 5: 
                # get unique 5 char compound_id
                compound_id = unique_ccd(ccd_db=ccd_db, len=5)
                compound_rec_copy.loc[len(compound_rec_copy)] = {"vault_id": vault_id, "compound_id": compound_id}
            elif len(vault_id) <= 5 and compound_id is None: # case where vault_id is valid 
                compound_id = vault_id
                compound_rec_copy.loc[len(compound_rec_copy)] = {"vault_id": vault_id, "compound_id": compound_id}
            elif compound_id is not None and len(compound_id) <= 5: # case where compound_id is valid
                # no need to update the record 
                compound_id = compound_id

            smiles_no_lg, lig_atom, wh_type = remove_leaving_group(lig[1])

            # makes pkl file if dne
            process_covalent_smiles(ccd_db, smiles_no_lg, compound_id) 
            writer.writerow([smiles_no_lg, compound_id, vault_id, wh_type, lig_atom, str(name), seq, int(res_idx), res_name, res_atom])

    # need to update records 
    compound_rec_copy.to_csv(COMPOUND_RECORD, index=False)  