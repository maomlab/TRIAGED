import os
import csv
import requests
from triage_boltz.covalent_utils import get_link_atoms

def fetch_fasta_from_pdb(pdb_id: str) -> str:
    """
    Fetches and returns the header anad FASTA sequence(s) for a given PDB ID from RCSB PDB.
    """
    url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
    response = requests.get(url)
    
    if response.status_code == 200:
        return response.text
    else:
        raise ValueError(f"Failed to fetch FASTA for PDB ID {pdb_id}: {response.status_code}")


def read_csv_pdbs(csv_path: str) -> list:
    """
    Reads a CSV file and returns a list of PDB IDs from the 'PDB' column.
    """
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        pdb_list = []
        for row in reader:
            pdb_list.append(row['PDB'])

    return pdb_list 

def build_fasta_dict(pdb_list: list) -> dict[str, str]:
    """
    Fetches FASTA sequences for a list of PDB IDs and stores them in a dictionary.

    Usage:
        pdbs = ['6ALZ', '7LZW']
        fasta_dict = build_fasta_dict(pdbs)
        Returns: 
            {'6WVO': '>6WVO_1|Chains A, B|Acetylcholinesterase...\nGREDAELLVTVRGGRLRGIRLKTPGGPVSA...\n',
            '7LZW': '>7LZW_1|Chains A, B|3C-like proteinase|... \nSNIGSGFRKMAFPSGKVEGCMVQVTCGTTTL...\n'}
    """
    fasta_dict = {}

    for pdb in pdb_list:
        try:
            fasta = fetch_fasta_from_pdb(pdb)
            fasta_dict[pdb] = fasta
            print(pdb)
        except Exception as e:
            print(f"Failed to fetch FASTA for {pdb}: {e}")

    return fasta_dict

def build_fasta_seq(pdb_id: str) -> str:
    """
    Returns fasta sequences of all chains as a single string which can be used as input for Boltz inference. 
    """
    fasta = fetch_fasta_from_pdb(pdb_id)
    list_fasta = fasta.split('\n')
    final_fasta = ''
    for i in list_fasta:
        if 'Chain' not in i: 
            final_fasta += i 

    return final_fasta

def build_covalent_chains_dict(pdb_list: list[str], pdb_dirs: str) -> dict[str, tuple]:
    """
    Extracts protein-ligand interaction information from LINK records in PDB files.
    Parameters:
        pdb_list (list of str): A list of PDB IDs (e.g., ['6ALZ', '7LZW']).
        pdb_dirs (str): Path to the directory containing subdirectories for each PDB ID.
                        pdb_dirs/
                        ├── 1XYZ/
                        │   └── 1XYZ.pdb

    Returns:
        prot_chains (dict): {PDB_ID: (prot_atom, res_name, res_num, chain_name, ccd, lig_atom, lig_idx), ...}
            Example:
              {'6ALZ': ('SG', 'CYS', 120, 'A', 'BKM', 'C01', 127), ...}
    """
    prot_chains = {}

    for pdb in pdb_list:
        pdb_dir = os.path.join(pdb_dirs, pdb) 
        pdb_file = os.path.join(pdb_dir, f"{pdb}.pdb") 
        try:
            prot_atom, res_name, res_num, chain_name, ccd, lig_atom, lig_idx = get_link_atoms(parent_file=pdb_file)
            prot_chains[pdb] = prot_atom, res_name, res_num, chain_name, ccd, lig_atom, lig_idx 
        except Exception as e: 
            print(pdb, f": {e}")

    return prot_chains


def build_fasta_dict_for_cov_inference(fasta_dict: dict, prot_chains: dict) -> tuple:
    """
    Extracts the amino acid sequence of the protein chain that covalently interacts 
    with the ligand for each PDB entry, filtering out peptide ligands.

    Parameters:
        fasta_dict (dict): A dictionary mapping PDB IDs to their full FASTA string retrieved from the RCSB PDB.
                        Example: {'6ALZ': '>6ALZ_1|Chains A|...\nMASEQUENCEHERE...', ...}

        prot_chains (dict): A dictionary mapping PDB IDs to tuples containing ligand interaction details.
                        Example: {'6ALZ': ('SG', 'CYS', 120, 'A', 'BKM', 'C01', 127), ...}

    Returns:
        tuple:
            final_fasta (dict): A dictionary mapping PDB IDs to their corresponding 
                                full protein sequence (string) of the interacting chain.
            errors (list): A list of PDB IDs that were skipped
    """

    errors = []
    final_fasta = {}
    for key, value in fasta_dict.items():
        list_fasta = value.split('\n')
        new_list_fasta = []
        for item in list_fasta:
            if item.startswith('>'):
                new_list_fasta.extend([item.split('|')])  # replace with split parts
            else:
                new_list_fasta.append(item)  
        for sub_items in new_list_fasta:
            for i in sub_items:
                if 'Chain' in i:
                    if prot_chains[key][3] in i: # get ligand interacting chain 
                        seq_idx = new_list_fasta.index(sub_items) + 1
        
                        if len(new_list_fasta[seq_idx]) > 20: # leaving out peptide ligands 
                            final_fasta[key] = new_list_fasta[seq_idx] 
                        else:
                            print(f'peptide ligand found or interacting chain not found for {key}')
                            errors.append(key)
    return final_fasta, errors


