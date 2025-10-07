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
from .pdb_to_fasta import residue_to_one_letter
from rdkit.Chem import CombineMols, rdMolTransforms
from pdbeccdutils.core.component import ConformerType # might need to use ccd_pkl env to run this 

########################################################################
# PROTEIN+LIGAND Processing
########################################################################

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
    '''
    Extracts either the protein or ligand components from a PDB file and writes them to a new PDB file.

    :param pdb_file: str 
        Path to the input PDB file containing the full protein-ligand complex.
    :param protein: bool
        If True, extracts the protein component (standard amino acids only).
        If False, extracts the ligand(s), excluding standard amino acids, solvents, and ions.
                        
    :return: str 
        Path to the output PDB file containing the extracted entity.
        PDB file in the same directory as the input, with suffix '_prot.pdb' or '_lig.pdb' based on the extracted entity.
    '''
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", pdb_file)
    output_base = Path(pdb_file).stem
   
    io = PDB.PDBIO()
    io.set_structure(structure)

    if protein:

        output = os.path.join(os.path.dirname(pdb_file), f"{output_base}_prot.pdb")
        io.save(output, ProteinSelect())
        print(f"[SUCCESS] Protein written to {output}")

    else: # ligand
        output = os.path.join(os.path.dirname(pdb_file), f"{output_base}_lig.pdb")
        io.save(output, LigandSelect())
        print(f"[SUCCESS] Ligand written to {output}")

    return output 

def combine_lig_prot(protein, ligand):
    '''
    Combines a protein and ligand RDKit molecule into a single molecule and returns its conformer.

    :param protein: rdkit.Chem.Mol
        RDKit Mol object representing the protein.
    :param ligand: rdkit.Chem.Mol
        RDKit Mol object representing the ligand.

    :return: tuple
        rdkit.Chem.Mol: Combined RDKit Mol object containing both protein and ligand.
        rdkit.Chem.rdchem.Conformer: Conformer object of the combined molecule containing 3D coordinates.
    '''
    combined = CombineMols(protein, ligand)
    conf = combined.GetConformer()
    return combined, conf


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
    raise ValueError("[ERROR] ", msg)


########################################################################
# PROTEIN Processing
########################################################################

def find_atom_index(mol, atom_name, res_num):
    '''
    Finds the atom index in an RDKit molecule given the atom name and residue number.

    :param mol: rdkit.Chem.rdchem.Mol
        RDKit molecule object (typically from a PDB).
    :param atom_name: str
        The name of the atom (e.g., "SG", "C8").
    :param res_num: int
        The residue number in the PDB (e.g., 25 for protein, 301 for ligand).

    :returns: int or None
        The atom index in the molecule if found, otherwise None.
    '''
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info and info.GetResidueNumber() == res_num and info.GetName().strip() == atom_name:
            return atom.GetIdx() # rdkit_idx 
    return None

def pdb_to_map(pdb_file, records_csv):
    '''
    Converts a PDB file to a sequence and gives the PDB idx of each residue as well as the literal index.

    :param pdb_file: str
        Path to PDB file.
    :param records_csv: str
        Path to CSV of CovalentInDB2.0 Database with information regarding positions of the covalent bond.

    :return res_info: list
        Map for indexing. List with tuples of (residue name, residue idx, literal index).
    '''
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
    '''
    Get residue idx and chain id interacting with a ligand for a given PDB ID.

    :param pdb_id: str
        PDB id of the query protein.
    :param records_csv: str
        Path to CSV of CovalentInDB2.0 Database with information regarding positions of the covalent bond.

    :return: tuple
        Tuple of covalent residue positon, ligand interacting chain and residue name. 
    '''
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
    '''
    Converts pdb indexing to boltz indexing and finds the Boltz index of the covalent residue. 

    :param pdb_file: str
        Path to PDB file in question.
    :param records_csv: str
        Path to CSV of CovalentInDB2.0 Database with information regarding positions of the covalent bond.

    :return mapped_idx: int
        Residue index based on Boltz mapping protocol.
    '''
    pdb_id = Path(pdb_file).stem
    res_info = pdb_to_map(pdb_file, records_csv)
    res_idx, _, res_name = res_idx_chain(pdb_id, records_csv)
    
    mapped_idx = [item for item in res_info if item[0] == residue_to_one_letter(res_name) and item[1] == res_idx]

    if len(mapped_idx) == 1:
        return mapped_idx[0][-1]
    else: 
        print(f'[WARNING] More than one matching residue found for {pdb_file}')

def get_link_atoms(parent_file, records_csv):
    '''
    Parses the LINK record from a PDB file to extract information about a covalent bond between a protein atom and a ligand atom.

    :param parent_file: str
        Path to the PDB file containing a LINK record describing the covalent bond.
        In file example: LINK         SG  CYS A  25                 C8  7KH A 301     1555   1555  1.76  
    :param records_csv: str
        Path to CSV of CovalentInDB2.0 Database with information regarding positions of the covalent bond.

    :return: tuple
            prot_atom: str
                Name of the protein atom involved in the covalent bond (e.g., "SG").
            res_name: str 
                Residue name of the protein atom (e.g., "CYS").
            res_idx: int
                Residue number of the protein atom (e.g., 25).
            chain_name: str
                Chain identifier for the protein residue (e.g., "A").
            lig_atom: str
                Name of the ligand atom involved in the covalent bond (e.g., "C8").
            ccd: str
                Three-letter chemical component ID of the ligand (e.g., "7KH").
            lig_idx: int
                Residue number (HETATM resSeq) of the ligand (e.g., 301).
    '''
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
    print(f"[WARNING] LINK record not present or no valid ligand found in PDB, {parent_file}")
    return None

##################################################################################
# Covalent LIGAND functions
##################################################################################

WARHEAD_REACTIONS = { "nitrile": "[C:3][C:4]#[N:5]>>[C:3][13C:4]=[N:5]", 
"alkylhalide" : "[CX4;CH,CH2:2][I,Br,Cl:3]>>[13C:2]",
"vinyl-sulfone" : "[C:3]=[C:4][S:5](=O)=O>>[13C:3][C:4][S:5](=O)=O", # for CYS rxn; might be diff for HIS (Schneider, Grabowsky 2015)
"acrylamide" : "[C:2]=[C:3]-C(=O)-[N:4]>>[13C:3]-[C:2]-C(=O)-[N:4]",
"nitrile2": "[N:4]#[C:5]>>[N:4]=[13C:5]"
}
WARHEAD_REACTANTS = {name: smarts.split(">>")[0] for name, smarts in WARHEAD_REACTIONS.items()}
compiled_reactants = {name: Chem.MolFromSmarts(smarts) for name, smarts in WARHEAD_REACTANTS.items()}

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
    raise RuntimeError("[ERROR] Could not find an available CCD after max attempts.")

def process_covalent_smiles(smiles, ccd_db="/home/ymanasa/.boltz/mols"):
    '''
    Creates a covalent CCD ligand from a conanical SMILES string using Rdkit.

    :param smiles: str
        Sanitized canonical SMILES string of the covalent ligand.
    :param ccd_db: str
        Path to the CCD database directory where .pkl files are stored.

    :return: str
        Unique CCD code for the covalent ligand.
    '''
    mol_sdf = Chem.MolFromSmiles(smiles) 

    success = compute_3d(mol_sdf)
    if not success:
        raise ValueError("[ERROR] 3D coordinate generation failed for the provided SMILES.")
    _ = get_conformer(mol_sdf, ConformerType.Computed)

    for i, atom in enumerate(mol_sdf.GetAtoms()):
        atom_id = atom.GetSymbol()+ str(atom.GetIdx())
        atom.SetProp("name", atom_id)
    
    Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)
    ccd = unique_ccd(ccd_db=ccd_db)
    with open(f"{ccd_db}/{ccd}.pkl", "wb") as f:
        pickle.dump(mol_sdf, f)
    return ccd

def identify_warhead(smiles):
    '''
    Identifies and returns the names of covalent warheads found in a given molecule.
    
    :param smiles: str
        The SMILES representation of the molecule.

    :returns: list
        A list of warhead names found in the molecule.
    '''

    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        print("[ERROR] Invalid SMILES string.")
        return []
        
    found_warheads = []
    for name, pattern in compiled_reactants.items():
        if pattern and mol.HasSubstructMatch(pattern):
            found_warheads.append(name)

    if len(found_warheads) > 1:
        print("[WARNING] More than 1 warhead found. Choosing first match.")
    elif len(found_warheads) == 0:
        raise ValueError(f'[ERROR] No matching warhead was found for {smiles}')
        
    return found_warheads[0]

def ligand_cov_atom(no_lg_smiles):
    '''
    Finds the index of the first Carbon-13 ([13C]) atom in a SMILES string.

    :param smiles_string: str
        The SMILES string to parse.

    :returns: int
        The 0-based index of the first C13 atom found.
        Returns -1 if no C13 atom is found or the SMILES is invalid.
    '''
    mol = Chem.MolFromSmiles(no_lg_smiles)

    if not mol:
        print("[ERROR] Invalid SMILES string provided.")
        return -1

    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetIsotope() == 13:
            return f'C{atom.GetIdx()}'
    return -1

def remove_leaving_group(smiles):
    '''
    Apply a covalent warhead-specific reaction to remove the leaving group
    from a ligand SMILES string.

    :param smiles: str
        Input ligand SMILES string.

    :returns: tuple
            smi_no_c13: str
                Processed SMILES with leaving group removed and isotopic labels stripped.
            lig_atom: int
                Index/identifier of the covalent attachment atom.
            wh_type: str
                Identified warhead type. 
    '''
    mol = Chem.MolFromSmiles(smiles)
    wh_type = identify_warhead(smiles)

    rxn = AllChem.ReactionFromSmarts(WARHEAD_REACTIONS[wh_type])

    products = rxn.RunReactants((mol,))

    for prod_tuple in products: # returns ((<rdkit.Chem.rdchem.Mol object at 0x14d267c24040>,),) 
        for prod in prod_tuple:
            smi_no_lg = Chem.MolToSmiles(prod)

    lig_atom = ligand_cov_atom(smi_no_lg)
    smi_no_c13 = smi_no_lg.replace("13", "")

    print('[SUCCESS] Leaving group removed:', smi_no_c13)

    return smi_no_c13, lig_atom, wh_type
    
def verify_covalent(res_name):
    '''
    Checks if a residue is a covalent amino acid.

    :param res_name: Three-letter residue code (e.g., 'CYS').

    :return (bool): True if covalent-capable, False otherwise.
    '''
    cov_aa = {"CYS", "SER", "THR", "LYS", "HIS", "PRO", "TYR", "GLU", "ASP"}
    return res_name.upper() in cov_aa
   
def residue_cov_atom(res_name):
    '''
    Given a residue name, return the atom commonly used for covalent bonding.

    :param res_name: Three-letter residue code.
    
    :return: Atom name string.
    '''
    mapping = {
        "CYS": "SG",
        "SER": "OG",
        "THR": "OG1",
        "LYS": "NZ",
        "HIS": "NE2",  # could also be NE1, ND1
        "PRO": "N",
        "TYR": "OH",
        "GLU": "OE2",
        "ASP": "OD2",
    }
    try:
        return mapping[res_name]
    except KeyError:
        raise ValueError(f"[ERROR] Residue '{res_name}' not supported for covalent docking.")

def ccd_is_ligand(ccd):
    '''
    Returns True if the CCD parsed is a ligand object and not a solvent, ion or amino acid
    '''   
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
