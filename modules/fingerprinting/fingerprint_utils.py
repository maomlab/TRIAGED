from triage_biomolecule.triage_biomolecule import TriageBiomolecule
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import numpy as np
import oddt
from oddt.interactions import compute_interactions
import pandas as pd
import prolif
from prolif.molecule import mol2_supplier, Molecule
from prolif.utils import to_bitvectors, to_countvectors
from tqdm import tqdm
import time

def get_ecfp_from_biomolecule(ligand: TriageBiomolecule, fingerprint_type: str = "ecfp4") -> DataStructs.ExplicitBitVect:
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
    
    ligand.ecfp = fp
    return fp 

def get_plif_from_biomolecule(ligand: TriageBiomolecule, receptor: TriageBiomolecule, fingerprint_type: str = "plif"):
    """
    Calculates a PLIF fingerprint for a given biomolecule. Not implimented yet
    """
    if ligand.entity_type != "ligand":
        raise ValueError("Only Ligand entity types can be fingerprinted")
    if receptor.entity_type != "protein":
        raise ValueError("Receptor must be a protein for fingerprint calculation")
    
    receptor = next(oddt.toolkit.readfile("pdb", ligand.pdb_path))
    ligand = next(oddt.toolkit.read)
    # Placeholder for actual PLIF calculation logic
    pass
def get_prolif_from_biomolecule(ligand: TriageBiomolecule, receptor: TriageBiomolecule,
                                 format="bitvectors", interactions="all") -> pd.DataFrame:
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
        corresponds to a residue–interaction pair. Values are boolean 
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
    ligand_poses = []
    for ligand_pose in mol2_supplier(ligand.pose_path):
        if ligand_pose is None:
            raise ValueError(f"Invalid pose in {ligand.pose_path}")
        ligand_poses.append(ligand_pose)

    # Compute the fingerprints for all ligand poses against the receptor
    fp.run_from_iterable(ligand_poses, receptor_molecule)
    
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
 

def create_fingerprint_dataframe(ligands: list, receptor: TriageBiomolecule,
                                  fingerprint_type: str = "prolif",
                                  format="bitvectors",
                                  interactions="all") -> pd.DataFrame: 
    """
    Create a DataFrame of fingerprints for a list of TriageBiomolecule ligands against a receptor.
    
    Parameters
    ----------
    ligands : list of TriageBiomolecule
        List of TriageBiomolecule objects representing ligands.
    receptor : TriageBiomolecule
        TriageBiomolecule object representing the receptor.
    fingerprint_type : str, optional
        Type of fingerprint to calculate: 'ecfp4', 'ecfp6', or 'prolif' (default).
    format : str, optional
        Format for ProLIF fingerprints: 'bitvectors' (default), 'countvectors', or 'dataframe'.
    interactions : list of str or 'all', optional
        Interaction types for ProLIF fingerprints. Use 'all' (default) for all types.
    
    Returns
    -------
    pandas.DataFrame
        A DataFrame where each row corresponds to a ligand and each column to a fingerprint feature.
    
    Notes
    -----
    For ECFP fingerprints, the DataFrame contains binary features indicating the presence of substructures.
    For ProLIF fingerprints, the DataFrame structure depends on the specified format and interactions.
    """
    fingerprint_list = []
    ligand_ids = []
    
    start_time = time.time()
    for ligand in tqdm(ligands, desc="Calculating fingerprints", unit="ligand"):
        if fingerprint_type in ["ecfp4", "ecfp6"]:
            fp = get_ecfp_from_biomolecule(ligand, fingerprint_type)
            arr = np.zeros((1,), dtype=int)
            DataStructs.ConvertToNumpyArray(fp, arr)
            fingerprint_list.append(arr)
            ligand_ids.append(ligand.entity_id)
        elif fingerprint_type == "prolif":
            df = get_prolif_from_biomolecule(ligand, receptor, format=format, interactions=interactions)
            fingerprint_list.append(df.values[0])  # grabbing the first fingerprint (pose) this is the usual case
            ligand_ids.append(ligand.entity_id)
        else:
            raise ValueError(f"Unknown fingerprint type: {fingerprint_type}. Use 'ecfp4', 'ecfp6', or 'prolif'.")
        
        elapsed_time = time.time() - start_time
        avg_time_per_ligand = elapsed_time / (ligands.index(ligand) + 1)
        remaining_time = avg_time_per_ligand * (len(ligands) - (ligands.index(ligand) + 1))
        tqdm.write(f"Estimated time remaining: {remaining_time:.2f} seconds")
    
    # Combine all fingerprints into a single DataFrame
    fingerprint_df = pd.DataFrame(fingerprint_list, index=ligand_ids)
    return fingerprint_df

def store_fingerprints_from_df(fingerprint_df: pd.DataFrame, 
                               biomolecules: list[TriageBiomolecule], 
                               fingerprint_type: str = "ecfp4") -> None:
    """
    Store fingerprints from a DataFrame into a TriageBiomolecule object.
    
    Parameters
    ----------
    fingerprint_df : pd.DataFrame
        DataFrame containing fingerprints with ligand IDs as index.
    biomolecule : TriageBiomolecule
        TriageBiomolecule object to store the fingerprints in.
    fingerprint_type : str, optional
        Type of fingerprint being stored (default is 'ecfp4').
    
    Notes
    -----
    The DataFrame index should match the entity IDs of the biomolecule ligands.
    """
    if biomolecule.entity_type != "ligand":
        raise ValueError("Only Ligand entity types can store fingerprints")
    
    for fingerprint_df_id in fingerprint_df.index:
        for ligand in biomolecules:
            if ligand.entity_id == fingerprint_df_id:
                ligand.store_fingerprint(fingerprint_df.loc[fingerprint_df_id].values, fingerprint_type)
                break
