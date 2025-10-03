
import os
import re
import csv
import random
import string
import pickle
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import pandas as pd
from pathlib import Path
from Bio import PDB
from Bio.PDB import PDBParser, Superimposer
from pdb_to_fasta import residue_to_one_letter
from rdkit.Chem import CombineMols, rdMolTransforms
from pdbeccdutils.core.component import ConformerType # might need to use ccd_pkl env to run this 

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

    :param pdb_file (str): Path to the input PDB file containing the full protein-ligand complex.
    :param protein (bool): If True, extracts the protein component (standard amino acids only).
                        If False, extracts the ligand(s), excluding standard amino acids, solvents, and ions.

    :return (str): Path to the output PDB file containing the extracted entity.
    :output: PDB file in the same directory as the input, with suffix '_prot.pdb' or '_lig.pdb' based on the extracted entity.
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

def pdb_to_map(pdb_file, records_csv):
    """
    Converts a PDB file to a sequence and gives the PDB idx of each residue as well as the literal index.

    :param pdb_file: Path to PDB file.
    :param records_csv: Path to CSV of CovalentInDB2.0 Database with information regarding positions of the covalent bond.

    :return res_info: map for indexing. list with tuples of (residue name, residue idx, literal index).
    """
    pdb_id = Path(pdb_file).stem # pdb code to look up 
    _, res_chain, _ = res_idx_chain(pdb_id, records_csv)

    index = 0
    res_info = []
    with open(pdb_file, 'r') as pdb:
        for line in pdb:
            # Extract residue information from ATOM records
            if line.startswith("ATOM") and line[13:15].strip() == "CA" and line[21:23].strip() == res_chain:

                residue = line[17:20].strip()
                res_name = residue_to_one_letter(residue)

                res_idx = line[23:27].strip() # pdb assigned
                res_idx = int(re.sub(r'[A-Za-z]', '', res_idx))
                
                index += 1 # literal counter 

                res_info.append((res_name, res_idx, index))

    return res_info 

def res_idx_chain(pdb_id, records_csv):
    """
    Get residue idx and chain id interacting with a ligand for a given PDB ID.

    :param pdb_id: PDB id of the query protein.
    :param records_csv: Path to CSV of CovalentInDB2.0 Database with information regarding positions of the covalent bond.

    :return: Tuple of covalent residue positon, ligand interacting chain and residue name. 
    """
    with open(records_csv, 'r') as f:
        reader = csv.reader(f)
        headers = next(reader)  # Get header row
    
        pdb_col = headers.index('PDB')
        resi_pos = headers.index('Resi_posi')
        chain_col = headers.index('Resi_chain')
        resi_name = headers.index('Resi_name')

        for row in reader:
            if row[pdb_col] == pdb_id:
                return int(row[resi_pos]), row[chain_col], row[resi_name]
    
def convert_pdbIDX_boltzIDX(pdb_file, records_csv):
    """
    Converts pdb indexing to boltz indexing and finds the Boltz index of the covalent residue. 
    :param pdb_file (str): Path to PDB file in question.
    :param records_csv (str): Path to CSV of CovalentInDB2.0 Database with information regarding positions of the covalent bond.

    :return mapped_idx (int): residue index based on Boltz mapping protocol.
    """
    pdb_id = Path(pdb_file).stem
    res_info = pdb_to_map(pdb_file, records_csv)
    res_idx, _, res_name = res_idx_chain(pdb_id, records_csv)
    
    mapped_idx = [item for item in res_info if item[0] == residue_to_one_letter(res_name) and item[1] == res_idx]

    if len(mapped_idx) == 1:
        return mapped_idx[0][-1]
    else: 
        print(f'more than one matching residue found for {pdb_file}')

def get_link_atoms(parent_file, records_csv):
    """
    Parses the LINK record from a PDB file to extract information about a covalent bond between a protein atom and a ligand atom.
    :param parent_file (str): Path to the PDB file containing a LINK record describing the covalent bond.
    In file example: LINK         SG  CYS A  25                 C8  7KH A 301     1555   1555  1.76  
    :records_csv (str): Path to CSV of CovalentInDB2.0 Database with information regarding positions of the covalent bond.

    :return (tuple):
            prot_atom (str): Name of the protein atom involved in the covalent bond (e.g., "SG").
            res_name (str): Residue name of the protein atom (e.g., "CYS").
            res_idx (int): Residue number of the protein atom (e.g., 25).
            chain_name (str): Chain identifier for the protein residue (e.g., "A").
            lig_atom (str): Name of the ligand atom involved in the covalent bond (e.g., "C8").
            ccd (str): Three-letter chemical component ID of the ligand (e.g., "7KH").
            lig_idx (int): Residue number (HETATM resSeq) of the ligand (e.g., 301).
    """
    with open(parent_file, 'r') as pdb:
        for line in pdb:
            if line.startswith("LINK"):
                # follows 0-indexing and end num in range not included
                ccd = line[47:51].strip() 
                if ccd_is_ligand(ccd): 
                    # protein
                    prot_atom = line[13:17].strip()
                    pdb_id =  pdb_id = Path(parent_file).stem
                    _, chain_name, res_name  = res_idx_chain(pdb_id, records_csv)

                    res_idx = convert_pdbIDX_boltzIDX(parent_file, records_csv) # fixes res_idx
                    # ligand
                    lig_atom = line[43:46].strip() 
                    lig_idx = int(line[22:26].strip()) # ?not needed?
                    return prot_atom, res_name, res_idx, chain_name, ccd, lig_atom, lig_idx
                else: # get info using alternate pdb format where ccd comes before covalent residue columns in LINK record
                    ccd = line[17:21].strip()
                    if ccd_is_ligand(ccd):
                        # protein
                        prot_atom = line[43:46].strip()
                        _, chain_name, res_name  = res_idx_chain(parent_file, records_csv)
                
                        res_idx = convert_pdbIDX_boltzIDX(parent_file, records_csv) 
                        # ligand
                        lig_atom = line[13:17].strip()
                        lig_idx = line[23:27].strip() # ?not needed?
                        return prot_atom, res_name, res_idx, chain_name, ccd, lig_atom, lig_idx
    print(f"Warning: LINK record not present or no valid ligand found in PDB, {parent_file}")
    return None

def combine_lig_prot(protein, ligand):
    """
    Combines a protein and ligand RDKit molecule into a single molecule and returns its conformer.

    :param protein (rdkit.Chem.Mol): RDKit Mol object representing the protein.
    :param ligand (rdkit.Chem.Mol): RDKit Mol object representing the ligand.

    :return (tuple):
        rdkit.Chem.Mol: Combined RDKit Mol object containing both protein and ligand.
        rdkit.Chem.rdchem.Conformer: Conformer object of the combined molecule containing 3D coordinates.
    """
    combined = CombineMols(protein, ligand)
    conf = combined.GetConformer()
    return combined, conf

def get_distance(conf, rdkit_idx1, rdkit_idx2):
    """
    Measures the Euclidean distance between two atoms in a combined protein-ligand RDKit conformer.
    
    :param conf (rdkit.Chem.rdchem.Conformer): RDKit combined conformer object containing 3D coordinates.
    :param rdkit_prot_idx (int): RDKit atom index of the protein atom.
    :param rdkit_lig_idx (int): RDKit atom index of the ligand atom.

    :return float: Distance between the two atoms in Ångströms.
    
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

def get_conformer(mol, c_type: ConformerType):
    ''' 
    Retrieve an rdkit object for a deemed conformer.
    Taken from `pdbeccdutils.core.component.Component`.
    :param mol: Mol
        The molecule to process.
    :param c_type: ConformerType
        The conformer type to extract.
    :return: Conformer
        The desired conformer, if any.
    '''
    for c in mol.GetConformers():
        try:
            if c.GetProp("name") == c_type.name:
                return c
        except KeyError:  # noqa: PERF203
            pass

    msg = f"Conformer {c_type.name} does not exist."
    raise ValueError(msg)

def compute_3d(mol) -> bool:
    '''Generate 3D coordinates using EKTDG method.
    Taken from `pdbeccdutils.core.component.Component`.
    :param mol: Mol
        The RDKit molecule to process
    :return: bool 
        Whether computation was successful.
    '''
    options = rdkit.Chem.AllChem.ETKDGv3()
    options.clearConfs = False
    conf_id = -1

    try:
        conf_id = rdkit.Chem.AllChem.EmbedMolecule(mol, options)
        rdkit.Chem.AllChem.UFFOptimizeMolecule(mol, confId=conf_id, maxIters=1000)

    except RuntimeError:
        pass  # Force field issue here
    except ValueError:
        pass  # sanitization issue here

    if conf_id != -1:
        conformer = mol.GetConformer(conf_id)
        conformer.SetProp("name", ConformerType.Computed.name)
        conformer.SetProp("coord_generation", f"ETKDGv3")
        
        return True
    return False

def unique_ccd(length=5, max_attempts=10, ccd_db="/home/ymanasa/.bolts/mols"):
    '''Generates a unique alphanumeric string of specified length and ensures uniqueness.'''
    chars = string.ascii_uppercase + string.digits

    for _ in range(max_attempts):
        ccd = ''.join(random.choices(chars, k=length))
        if not os.path.exists(f"{ccd_db}/{ccd}.pkl"):
            return ccd
    raise RuntimeError("Could not find an available CCD after max attempts.")

def process_covalent_smiles(smiles, ccd_db="/home/ymanasa/.bolts/mols"):
    '''
    Creates a covalent CCD ligand from a conanical SMILES string using Rdkit.
    :param smiles (str): Sanitized canonical SMILES string of the covalent ligand.
    :param ccd_db (str): Path to the CCD database directory where .pkl files are stored.
    :return (str): Unique CCD code for the covalent ligand.
    '''
    mol_sdf = Chem.MolFromSmiles(smiles) 

    success = compute_3d(mol_sdf)
    if not success:
        raise ValueError("3D coordinate generation failed for the provided SMILES.")
    _ = get_conformer(mol_sdf, ConformerType.Computed)

    for i, atom in enumerate(mol_sdf.GetAtoms()):
        atom_id = atom.GetSymbol()+ str(atom.GetIdx())
        atom.SetProp("name", atom_id)
    
    Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)
    ccd = unique_ccd(ccd_db=ccd_db)
    with open(f"{ccd_db}/{ccd}.pkl", "wb") as f:
        pickle.dump(mol_sdf, f)
    return ccd

def remove_lg(smiles, warhead_type, smarts_rxns="SMARTS.csv"):  

    mol = Chem.MolFromSmiles(smiles)
    with open(smarts_rxns, 'r') as csv:
        reader = csv.DictReader(csv)
        for row in reader:
            if warhead_type == row["warhead_type"]:
                rxn = AllChem.ReactionFromSmarts(row["SMARTS"])
                try:
                    products = rxn.RunReactants((mol,))

                    for prod_tuple in products:
                        for prod in prod_tuple:
                            smi_no_lg = Chem.MolToSmiles(prod)
                     
                    return smi_no_lg
                except: 
                    raise ValueError("Reaction could not be applied to the given SMILES.")

    # use SMARTS to identify warhead and leaving group 
    # make new file with SMARTS reactions and map to warheads 
    # define warheads and patterns in README 
    # identify covalent atom as well dependening on warhead type
    # lig_atom just needs to index and see which atom is covaolently interating after leaving group is removed 
    # return lig_atom

   
def get_prot_atom(res_name):
    # given res name, return which atom is covalent 
    # have a dict with standard covalent residues and their covalent atom types
    pass 

def get_lig_atom():
    pass 

# testing
# get_link_atoms('/home/ymanasa/turbo/ymanasa/opt/boltz/covalent_testing/boltz_infer/6WVO/6WVO.pdb', '/home/ymanasa/turbo/ymanasa/opt/boltz/covalent_testing/records/records_test.csv')