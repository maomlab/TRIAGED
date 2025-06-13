import sys

def pdb_to_fasta(pdb_file, fasta_file):
    """
    Converts a PDB file to a FASTA file.
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
    return mapping.get(residue, "X")  # Default to 'X' for unknown residues

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python pdb_to_fasta.py <input_pdb_file> <output_fasta_file>")
        sys.exit(1)
    
    pdb_file = sys.argv[1]
    fasta_file = sys.argv[2]
    pdb_to_fasta(pdb_file, fasta_file)
