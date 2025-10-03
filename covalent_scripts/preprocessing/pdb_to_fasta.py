import sys
import re 

def build_sequence(pdb_file):
    """
    Converts a PDB file to a sequence.

    :param pdb_file: Path to PDB file.

    :return sequence: string of the protein sequence for the ligand interacting chain.
    """

    sequence = ''
    with open(pdb_file, 'r') as pdb:
        for line in pdb:
            # Extract residue information from ATOM records
            if line.startswith("ATOM") and line[13:15].strip() == "CA":

                residue = line[17:20].strip()
                res_name = residue_to_one_letter(residue)

                res_idx = line[23:27].strip() # pdb assigned
                res_idx = int(re.sub(r'[A-Za-z]', '', res_idx))
                
                sequence += res_name

    return sequence 

def pdb_to_fasta(pdb_file, fasta_file):
    """
    Converts a PDB file to a TXT file.
    :param pdb_file: Path to the input PDB file.
    :param fasta_file: Path to the output FASTA file.
    """
    sequence = []
    with open(pdb_file, 'r') as pdb:
        for line in pdb:
            # Extract residue information from ATOM records
            if line.startswith("ATOM") and line[13:15].strip() == "CA":
                residue = line[17:20].strip()
                sequence.append(residue_to_one_letter(residue))
    
    # Write sequence to FASTA file
    with open(fasta_file, 'w') as fasta:
        fasta.write("".join(sequence) + "\n")

def residue_to_one_letter(residue):
    """
    Converts a three-letter residue code to a one-letter code.
    :param residue: Three-letter residue code.
    :return: One-letter residue code.
    """
    mapping = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
    }
    return mapping.get(residue, "X")  

def residue_to_three_letter(residue):
    """
    Converts a one-letter residue code to a three-letter code.
    :param residue: One-letter residue code.
    :return: Three-letter residue code.
    """
    mapping = {
        "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
        "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
        "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
        "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL"
    }
    return mapping.get(residue.upper(), "UNK") 

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python pdb_to_fasta.py <input_pdb_file> <output_fasta_file>")
        sys.exit(1)
    
    pdb_file = sys.argv[1]
    fasta_file = sys.argv[2]
    pdb_to_fasta(pdb_file, fasta_file)
