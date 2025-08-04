from triage_biomolecule.triage_biomolecule import TriageBiomolecule
from rdkit import Chem
from rdkit.Chem import AllChem
import oddt
from oddt.interactions import compute_interactions
import pandas as pd
import prolif
from prolif.molecule import mol2_supplier, Molecule
from prolif.utils import to_bitvectors, to_countvectors

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

def get_prolif_from_biomolecule(ligands: list[TriageBiomolecule], receptor: TriageBiomolecule, format="bit", interactions="all"):
    """
    Generate a ProLIF interaction fingerprint DataFrame for multiple docking poses.
    
    Parameters
    ----------
    ligands : list of TriageBiomolecule
        List of TriageBiomolecule objects representing ligand poses.
    receptor : TriageBiomolecule
        TriageBiomolecule object representing the receptor.
    format : str, optional
        'bit' to obtain a boolean fingerprint (default) or 'count' to obtain counts of
        each interaction type for each residue pair.
    interactions : list of str or 'all', optional
        The list of interaction names to consider. Use 'all' (default) for every available
        interaction class.
    
    Returns
    -------
    pandas.DataFrame
        A DataFrame where each row corresponds to a ligand pose and each column
        corresponds to a residue–interaction pair. Values are boolean if
        `format='bit'` or integer counts if `format='count'`.
    
    Notes
    -----
    The receptor is loaded with RDKit and converted to a ProLIF `Molecule`.
    Each TriageBiomolecule ligand can contain multiple poses; all are processed.
    """
    # Load the receptor PDB into an RDKit Mol with residue information
    rec_rdkit = Chem.MolFromPDBFile(receptor.pdb_path, removeHs=False, sanitize=False)
    if rec_rdkit is None:
        raise ValueError(f"Could not read PDB file at {receptor.pdb_path}")
    # Wrap as a ProLIF Molecule; residue info from the PDB is kept by default
    receptor_molecule = Molecule.from_rdkit(rec_rdkit)
    
    # Create a fingerprint generator; use all interactions or a custom list
    if interactions == "all":
        fp = prolif.Fingerprint()
    else:
        fp = prolif.Fingerprint(interactions)
    
    # Gather ligands from the TriageBiomolecule objects
    ligands_molecules = []
    for ligand in ligands:
        
        for ligand_pose in mol2_supplier(ligand.pose_path):
            if ligand_pose is None:
                raise ValueError(f"Invalid pose in {ligand.pose_path}")
            ligands_molecules.append(ligand_pose)

    # Compute the fingerprints for all ligand poses against the receptor
    fp.run_from_iterable(ligands_molecules, receptor_molecule)
    
    # Convert to a DataFrame. The count flag controls whether a count or bit fingerprint is returned.
        # Get appropriate output format
    df = fp.to_dataframe(count=True, dtype=int)

    if format == "dataframe":
        return df.reset_index(drop=True)
    elif format == "bitvectors":
        return to_bitvectors(df)
    elif format == "countvectors":
        return to_countvectors(df)
    else:
        raise ValueError(f"Unknown format: {format}. Use 'dataframe', 'bitvectors', or 'countvectors'.")
    
    # The DataFrame index is the pose number; reset to ensure a clean integer index
    df = df.reset_index(drop=True)

    return df