from rdkit.Chem import CombineMols, rdMolTransforms
from Bio.PDB import PDBParser, Superimposer
from Bio import PDB
import os
from pathlib import Path

class ProteinSelect(PDB.Select):
    def accept_residue(self, residue):
        # Filter only standard amino acids
        return residue.get_resname().strip() in PDB.Polypeptide.standard_aa_names
    
class LigandSelect(PDB.Select):
    def accept_residue(self, residue):
        name = residue.get_resname().strip()
        standard_aa = PDB.Polypeptide.standard_aa_names
        standard_pdb_solvents = [
        "HOH", "DOD", "H2O", "ETH", "EOH", "IPA", "MPD", "DMS",
        "ACT", "ACE", "GOL", "PEG", "SO4", "PO4", "TRS", "MES",
        "HEP", "FMT", "TLA", "EDO"
        ]

        standard_pdb_ions = [
            "NA", "K", "CL", "CA", "MG", "ZN", "FE", "FE2", "FE3",
            "MN", "CO", "CU", "NI", "CD", "SR", "BA", "BR", "IOD",
            "CS", "RB"
        ]
        # Exclude solvents, ions and standard amino acids
        return name not in standard_aa and name not in standard_pdb_solvents and name not in standard_pdb_ions

def extract_entities(pdb_file, protein=True):
    """
    Extracts either the protein or ligand components from a PDB file and writes them to a new PDB file.

    Args:
        pdb_file (str): Path to the input PDB file containing the full protein-ligand complex.
        protein (bool): If True, extracts the protein component (standard amino acids only).
                        If False, extracts the ligand(s), excluding standard amino acids, solvents, and ions.

    Returns:
        str: Path to the output PDB file containing the extracted entity.

    Output:
        Writes a new PDB file in the same directory as the input, with suffix '_prot.pdb' or '_lig.pdb'
        based on the extracted entity.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", pdb_file)
    output_base = Path(pdb_file).stem
   
    io = PDB.PDBIO()
    io.set_structure(structure)

    if protein:

        output = os.path.join(os.path.dirname(pdb_file), f"{output_base}_prot.pdb")
        io.save(output, ProteinSelect())
        print(f"Protein written to {output}")

    else: # ligand
        output = os.path.join(os.path.dirname(pdb_file), f"{output_base}_lig.pdb")
        io.save(output, LigandSelect())
        print(f"Ligand written to {output}")

    return output 
    

def find_atom_index(mol, atom_name, res_num):
    """
    Finds the atom index in an RDKit molecule given the atom name and residue number.

    Args:
        mol (rdkit.Chem.rdchem.Mol): RDKit molecule object (typically from a PDB).
        atom_name (str): The name of the atom (e.g., "SG", "C8").
        res_num (int): The residue number in the PDB (e.g., 25 for protein, 301 for ligand).

    Returns:
        int or None: The atom index in the molecule if found, otherwise None.
    """
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info and info.GetResidueNumber() == res_num and info.GetName().strip() == atom_name:
            return atom.GetIdx() # rdkit_idx 
    return None

def get_residue_index_from_pdb(res_num, chain_id, pdb_file):
    """
    Given a residue number and chain ID from a PDB file, return its 0-based index 
    as used in the sequence passed to `parse_polymer`.

    Args:
        res_num (int): Residue number from the PDB (e.g., 25).
        chain_id (str): Chain ID (e.g., 'A').
        pdb_file (str): Path to the PDB file.

    Returns:
        int or None: The 0-based index of the residue in the chain.
    """
    seen = []
    with open(pdb_file, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                res_id = int(line[22:26])
                chain = line[21]
                if chain == chain_id and res_id not in seen:
                    seen.append(res_id)
    try:
        return seen.index(res_num)
    except ValueError:
        return None

def ccd_is_ligand(ccd):
    """
    Returns True if the CCD parsed is a ligand object and not a solvent, ion or amino acid
    """   

    standard_aa = PDB.Polypeptide.standard_aa_names

    standard_pdb_solvents = [
        "HOH", "DOD", "H2O", "ETH", "EOH", "IPA", "MPD", "DMS",
        "ACT", "ACE", "GOL", "PEG", "SO4", "PO4", "TRS", "MES",
        "HEP", "FMT", "TLA", "EDO"
        ]

    standard_pdb_ions = [
            "NA", "K", "CL", "CA", "MG", "ZN", "FE", "FE2", "FE3",
            "MN", "CO", "CU", "NI", "CD", "SR", "BA", "BR", "IOD",
            "CS", "RB", "FE2", "FE3"
        ]
    if ccd not in standard_pdb_solvents and ccd not in standard_pdb_ions and ccd not in standard_aa:
        return True
    
def get_link_atoms(parent_file):
    """
    Parses the LINK record from a PDB file to extract information about a covalent bond 
    between a protein atom and a ligand atom.

    Args:
        parent_file (str): Path to the PDB file containing a LINK record describing the covalent bond.
        In file example: LINK         SG  CYS A  25                 C8  7KH A 301     1555   1555  1.76  

    Returns:
        tuple:
            prot_atom (str): Name of the protein atom involved in the covalent bond (e.g., "SG").
            res_name (str): Residue name of the protein atom (e.g., "CYS").
            res_idx (int): Residue number of the protein atom (e.g., 25).
            chain_name (str): Chain identifier for the protein residue (e.g., "A").
            lig_atom (str): Name of the ligand atom involved in the covalent bond (e.g., "C8").
            ccd (str): Three-letter chemical component ID of the ligand (e.g., "7KH").
            lig_idx (int): Residue number (HETATM resSeq) of the ligand (e.g., 301).
    """
    import ipdb; ipdb.set_trace()
    with open(parent_file, 'r') as pdb:
        for line in pdb:
            if line.startswith("LINK"): # look at all LINK lines until we find cov ligand line 
                ccd = line[47:51].strip() 
                if ccd_is_ligand(ccd): 
                    prot_atom = line[13:17].strip()
                    res_name = line[17:21].strip() 
                    chain_name = line[21:23].strip() 

                    res_idx = int(line[23:27].strip()) # incorrect due to indexing by PDB
                    res_num = get_residue_index_from_pdb(res_idx, chain_name, parent_file) # resort to using this instead of res_idx

                    lig_atom = line[43:47].strip() 
                    lig_idx = int(line[22:26].strip()) 
                    return prot_atom, res_name, res_num, chain_name, ccd, lig_atom, lig_idx
                else: # ccd detected was a solvent, ion or aa
                    # try getting info using alternate pdb format where ccd comes before covalent residue columns in LINK record
                    ccd = line[17:21].strip()
                    if ccd_is_ligand(ccd):
                        prot_atom = line[43:47].strip()
                        res_name = line[47:51].strip() 
                        chain_name = line[51:53].strip() 

                        res_idx = line[53:57].strip()
                        res_num = get_residue_index_from_pdb(res_idx, chain_name, parent_file)

                        lig_atom = line[12:16].strip()
                        lig_idx = line[23:27].strip()
                        return prot_atom, res_name, res_num, chain_name, ccd, lig_atom, lig_idx
    print(f"Warning: LINK record not present or no valid ligand found in PDB, {parent_file}")
    return None

def combine_lig_prot(protein, ligand):
    """
    Combines a protein and ligand RDKit molecule into a single molecule and returns its conformer.

    Args:
        protein (rdkit.Chem.Mol): RDKit Mol object representing the protein.
        ligand (rdkit.Chem.Mol): RDKit Mol object representing the ligand.

    Returns:
        tuple:
            - rdkit.Chem.Mol: Combined RDKit Mol object containing both protein and ligand.
            - rdkit.Chem.rdchem.Conformer: Conformer object of the combined molecule containing 3D coordinates.
    """
    combined = CombineMols(protein, ligand)
    conf = combined.GetConformer()
    return combined, conf

def get_distance(conf, rdkit_idx1, rdkit_idx2):
    """
    Measures the Euclidean distance between two atoms in a combined protein-ligand RDKit conformer.
    
    Args:
        conf (rdkit.Chem.rdchem.Conformer): RDKit combined conformer object containing 3D coordinates.
        rdkit_prot_idx (int): RDKit atom index of the protein atom.
        rdkit_lig_idx (int): RDKit atom index of the ligand atom.

    Returns:
        float: Distance between the two atoms in Ångströms.
    
    Notes:
        - Use find_atom_index() to get rdkit idxs. 
        - Rdkit idxs are different from atom idxs listed in PDB files. 
    """
    return rdMolTransforms.GetBondLength(conf, rdkit_idx1, rdkit_idx2)

def get_ca_atoms_by_residues(structure):
    """Return a dictionary {(chain_id, resid): atom} for CA atoms"""
    ca_atoms = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    resid = residue.get_id()[1]  
                    chain_id = chain.id.strip()
                    key = (chain_id, resid)
                    ca_atoms[key] = residue['CA']
    return ca_atoms


def compute_rmsd(ref_path, pred_path):
    """
    Computes the Root-Mean-Square Deviation (RMSD) between two PDB structures 
    using C-alpha atoms of common residues.

    Args:
        ref_path (str): File path to the reference PDB structure.
        pred_path (str): File path to the predicted PDB structure.

    Returns:
        float: The RMSD value in Ångstroms if common residues are found.
        None: If no common residues are found between the structures.

    Notes:
        - Only residues with matching chain ID and residue number are compared.
        - Superimposes predicted structure onto the reference before computing RMSD.
    """
    parser = PDBParser(QUIET=True)
    ref_struct = parser.get_structure("ref", ref_path)
    pred_struct = parser.get_structure("pred", pred_path)

    ref_atoms_map = get_ca_atoms_by_residues(ref_struct)
    pred_atoms_map = get_ca_atoms_by_residues(pred_struct)

    # Get common residues
    common_keys = sorted(set(ref_atoms_map) & set(pred_atoms_map))
    print(f"Found {len(common_keys)} common residues for RMSD.")

    if len(common_keys) == 0:
        print("No overlapping residues found.")
        return

    ref_atoms = [ref_atoms_map[key] for key in common_keys]
    pred_atoms = [pred_atoms_map[key] for key in common_keys]

    # Superimpose and calculate RMSD
    sup = Superimposer()
    sup.set_atoms(ref_atoms, pred_atoms)
    print(f"RMSD: {sup.rms:.3f} Å")

    return sup.rms

def get_neighbor(atom, exclude_idx):
    """
    Returns the index of a neighboring atom, prioritizing atoms with 'C' in their PDB name,
    while excluding a specified atom index.

    Args:
        atom (rdkit.Chem.Atom): RDKit Atom object whose neighbors are to be examined.
        exclude_idx (int): Atom index to be excluded (e.g., the atom it's bonded to for a dihedral).

    Returns:
        int or None: Index of the selected neighboring atom, or None if no suitable neighbor is found.

    Notes:
        - Only neighbors with 'C' in their PDB atom name are considered (e.g., CB, CG).
        - This is a heuristic to prioritize carbon atoms in dihedral definitions.
    """
    for neighbor in atom.GetNeighbors():
        pdb_info = neighbor.GetPDBResidueInfo()
        name = pdb_info.GetName().strip()
        
        if 'C' in name:
            print(f"Neighbor name: {name}")
            if neighbor.GetIdx() != exclude_idx:
                return neighbor.GetIdx()
    return None

def get_dihedral(combined, conf, atomB_idx, atomC_idx):
    """
    Measures the dihedral (torsion) angle defined by four atoms: A-B-C-D, where B and C are specified,
    and A and D are their respective bonded neighbors.

    Args:
        combined (rdkit.Chem.Mol): RDKit molecule containing both protein and ligand.
        conf (rdkit.Chem.rdchem.Conformer): Conformer object of the combined molecule with 3D coordinates.
        atomB_idx (int): Index of atom B (typically the protein atom in the covalent bond).
        atomC_idx (int): Index of atom C (typically the ligand atom in the covalent bond).

    Returns:
        float: Dihedral angle in degrees.

    Raises:
        ValueError: If a neighboring atom (A or D) cannot be found to define the dihedral angle.
    """
    atom_B = combined.GetAtomWithIdx(atomB_idx)
    atom_C = combined.GetAtomWithIdx(atomC_idx)
    # get neighboring atoms
    atom_A_idx = get_neighbor(atom_B, atom_C.GetIdx())
    atom_D_idx = get_neighbor(atom_C, atom_B.GetIdx())

    # check
    if atom_A_idx is None or atom_D_idx is None:
        raise ValueError("Couldn't find bonded neighbor to define dihedral.")
    
    #dihedral
    return rdMolTransforms.GetDihedralDeg(conf, atom_A_idx, atomB_idx, atomC_idx, atom_D_idx)

# testing
get_link_atoms("/home/ymanasa/turbo/ymanasa/opt/boltz/covalent_testing/4EBP/4EBP.pdb")