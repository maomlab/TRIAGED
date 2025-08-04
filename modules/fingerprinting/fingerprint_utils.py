from triage_biomolecule.triage_biomolecule import TriageBiomolecule
from rdkit import Chem
from rdkit.Chem import AllChem
import oddt
from oddt.interactions import compute_interactions

def get_ecfp_from_biomolecule(ligand: TriageBiomolecule, fingerprint_type: str = "ecfp4"):
    """
    Calculates a fingerprint for a given biomolecule. 
    input: ligand - TriageBiomolecule object representing the ligand
    fingerprint_type - Type of fingerprint to calculate, either "ecfp4" or "ecfp6"
    output: A fingerprint object representing the ligand
    """
    if ligand.entity_type != "ligand":
        raise ValueError("Only Ligand entity types can be fingerprinted")
    
    ligand_smiles = ligand.smiles
    if not ligand_smiles:
        raise ValueError("Ligand must have a valid SMILES representation")
    mol = Chem.MolFromSmiles(ligand_smiles)
    if not mol:
        raise ValueError("Invalid SMILES representation for the ligand")
    if fingerprint_type == "ecfp4":
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, nBits=2048, radius=2)
    elif fingerprint_type == "ecfp6":
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, nBits=2048, radius=3)
    return fp 

def get_plif_from_biomolecule(ligand: TriageBiomolecule, receptor: TriageBiomolecule, fingerprint_type: str = "plif"):
    """
    Calculates a PLIF fingerprint for a given biomolecule.
    """
    if ligand.entity_type != "ligand":
        raise ValueError("Only Ligand entity types can be fingerprinted")
    if receptor.entity_type != "protein":
        raise ValueError("Receptor must be a protein for fingerprint calculation")
    
    receptor = next(oddt.toolkit.readfile("pdb", ligand.pdb_path))
    ligand = next(oddt.toolkit.read)
    # Placeholder for actual PLIF calculation logic
    return None  # Replace with actual PLIF calculation result

def get_prolif_from_biomolecule():
    pass