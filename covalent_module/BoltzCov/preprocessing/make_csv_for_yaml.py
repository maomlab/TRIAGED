import pandas as pd
import os
from datetime import date
import csv 
import string 
import random 
from BoltzCov.preprocessing import covalent_utils, pdb_to_fasta

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
        sequence = pdb_to_fasta.build_sequence(pdb, lig_chain)
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

    res_name = pdb_to_fasta.residue_to_three_letter(res_aa)
    if VERBOSE: print("Boltz will now dock to this residue: ", res_name)

    if covalent_utils.verify_covalent(res_name) != True: # verifies if this residue can participate in a covalent bond w the
        if VERBOSE: print(sequence)
        raise ValueError(f"[ERROR] res_idx {idx} does NOT map to a covalent residue. " \
        "Please verify res_idx matches expected residue in sequence.")
    
    res_atom = covalent_utils.residue_cov_atom(res_name)

    return sequence, res_name, res_atom

def unique_ccd(ccd_db, len=5, max_attempts=1000):
    '''Generates a unique alphanumeric string of specified length and ensures uniqueness.'''

    chars = string.ascii_uppercase + string.digits

    for _ in range(max_attempts):
        ccd = ''.join(random.choices(chars, k=len))
        if not os.path.exists(f"{ccd_db}/{ccd}.pkl"):
            return ccd
    raise RuntimeError("[ERROR] Could not find a unique CCD ID after max attempts.")

def lookup_compound_id(substance_id, compound_record):
    match = compound_record.loc[compound_record['substance_id'] == substance_id, 'five_char_id']
    five_char_id = match.iloc[0] if not match.empty else None
    return five_char_id

def generate_csv(pdb, ligand_df, boltz_cache, res_idx, ligand_chain, VERBOSE):
    '''
    1. Process protein to get sequence, covalent residue name and atom.
    2. Ligand processing: 
        - leaving group removal 
        - pkl making 
    3. Make CSV to make yamls for docking.
    '''
    if VERBOSE: print("- Processing protein right now.")
    seq, res_name, res_atom = process_protein(pdb, res_idx, ligand_chain)

    # writing csv for yaml 
    if VERBOSE: print("- Writing CSV to make yamls for docking.")
    today = date.today()
    
    tmp_csv = os.path.join(boltz_cache, f'ligands_{today}.csv') 

    if os.path.exists(tmp_csv):
        if VERBOSE:  print(f"[WARNING] Output CSV '{tmp_csv}' already exists. Deleting and rewriting.")
        os.remove(tmp_csv)

    # write header once
    expected_header = ["smiles", "substance_id", "five_char_id" ,"WH_Type", "Lig_Atom", "Prot_ID", "Prot_Seq", "Res_Idx", "Res_Name", "Res_Atom"]
    write_header = True
    with open(tmp_csv, "r") as existing:
        reader = csv.reader(existing)
        first_row = next(reader, None)
        if first_row == expected_header:
            write_header = False  # won't rewrite header
    
    # ligand processing 
    if VERBOSE: print("- Processig ligands right now.")
    ligands = [(row['substance_id'], row['smiles']) for _, row in ligand_df.iterrows()]

    # setup to assign 5-char codes 
    compound_rec_df = pd.DataFrame(columns=['substance_id', 'five_char_id'])
    new_rows_list = [{'substance_id': name, 'five_char_id': 'XXXXXXX'} for name, _ in ligands]
    compound_rec_copy = pd.concat([compound_rec_df, pd.DataFrame(new_rows_list)], ignore_index=True)

    boltz_cache_pkls = os.path.join(boltz_cache, 'cache_pkls')
    os.makedirs(boltz_cache_pkls, exist_ok=True)

    with open(tmp_csv, 'a') as f: 
        writer = csv.writer(f)
        if write_header:
            writer.writerow(expected_header)
        # for each ligand, append the protein information, assuming one protein target 
        for lig in ligands:
            substance_id = lig[0]
            five_char_id = lookup_compound_id(substance_id, compound_rec_copy) 
            if five_char_id is None or len(five_char_id) > 5: 
                # get unique 5 char five_char_id
                five_char_id = unique_ccd(ccd_db=boltz_cache_pkls, len=5)
                compound_rec_copy.loc[len(compound_rec_copy)] = {"substance_id": substance_id, "five_char_id": five_char_id}
            elif len(substance_id) <= 5 and five_char_id is None: # case where substance_id is valid 
                five_char_id = substance_id
                compound_rec_copy.loc[len(compound_rec_copy)] = {"substance_id": substance_id, "five_char_id": five_char_id}
            elif five_char_id is not None and len(five_char_id) <= 5: # case where five_char_id is valid
                # no need to update the record 
                five_char_id = five_char_id

            smiles_no_lg, lig_atom, wh_found = covalent_utils.remove_leaving_group(lig[1])
            if wh_found is None: 
                print('skipping ligand', lig[0])
                continue 
            
             # makes pkl file if dne
            covalent_utils.process_covalent_smiles(ccd_db=boltz_cache_pkls, smiles=smiles_no_lg, compound_id=five_char_id) 
            # update tmp_csv
            writer.writerow([smiles_no_lg, five_char_id, substance_id, wh_found, lig_atom, seq, int(res_idx), res_name, res_atom])

    return tmp_csv 

        
    
                





