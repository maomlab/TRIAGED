import os
import csv
import requests
from covalent_utils import get_link_atoms

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