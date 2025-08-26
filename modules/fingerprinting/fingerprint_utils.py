#from ..triage_biomolecule.triage_biomolecule import TriageBiomolecule
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from fingerprinting.mmcif_utils import extract_protein_receptor, extract_ligands_to_sdf, extract_ligand
from triage_biomolecule.triage_biomolecule import TriageBiomolecule
from misc_utils.split_mmcif import split_mmcif
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import numpy as np
import oddt
import pandas as pd
import prolif
from prolif.molecule import mol2_supplier, Molecule
from prolif.utils import to_bitvectors, to_countvectors
from tqdm import tqdm
import time
import tempfile
from pathlib import Path
from typing import Iterable, Set
from Bio.PDB import MMCIFParser
from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Geometry import Point3D
from rdkit.Chem.rdmolfiles import SDWriter
from rdkit.Chem import rdDetermineBonds


from Bio.PDB import MMCIFParser, PDBIO, Select


from prolif.molecule import sdf_supplier
from openbabel import pybel
import subprocess

def get_ecfp_from_biomolecule(ligand: TriageBiomolecule, fingerprint_type: str = "ecfp4") -> DataStructs.ExplicitBitVect:
    """
    Calculates a fingerprint for a given biomolecule. 
    input: ligand - TriageBiomolecule object representing the ligand
    fingerprint_type - Type of fingerprint to calculate, either "ecfp4" or "ecfp6"ls
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
                                 format="countvectors", interactions="all") -> pd.DataFrame:
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
    rec_rdkit = Chem.MolFromPDBFile(ligand.paired_receptor_path, removeHs=False, sanitize=False)
    if rec_rdkit is None:
        raise ValueError(f"Could not read PDB file at {ligand.paired_receptor_path}")
    # Wrap as a ProLIF Molecule; residue info from the PDB is kept by default
    receptor_molecule = Molecule.from_rdkit(rec_rdkit)
    
    # Create a fingerprint generator; use all interactions or a custom list
    if interactions == "all":
        fp = prolif.Fingerprint()
    else:
        fp = prolif.Fingerprint(interactions)

   
    ligand_poses = []
    # Use prolif.molecule.sdf_supplier for SDF files
    for ligand_pose in sdf_supplier(ligand.pose_path):
        if ligand_pose is None:
            raise ValueError(f"Invalid pose in {ligand.pose_path}")
        ligand_poses.append(ligand_pose)

    # Compute the fingerprints for all ligand poses against the receptor
    fp.run_from_iterable(ligand_poses, receptor_molecule, progress=False)
    
    # Convert to a DataFrame. The count flag controls whether a count or bit fingerprint is returned.
        # Get appropriate output format
    df = fp.to_dataframe(count=True, dtype=int)
    print("prolif fingerprint here")
    print(df)
    if format == "dataframe":
        return df.reset_index(drop=True)
    elif format == "bitvectors":
        return to_bitvectors(df)
    elif format == "countvectors":
        return to_countvectors(df)
    else:
        raise ValueError(f"Unknown format: {format}. Use 'dataframe', 'bitvectors', or 'countvectors'.")

def fps_to_dense_matrix(fps):
    n_bits = fps.GetNumBits()
    X = np.zeros((len(fps), n_bits), dtype=np.uint8)
    for i, fp in enumerate(fps):
        DataStructs.ConvertToNumpyArray(fp, X[i])
    return X  

#creates fingerprints
def create_fingerprint_dataframe(ligands: list,
                                  fingerprint_type: str = "prolif_bit",
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
        elif fingerprint_type == "prolif_bit" or fingerprint_type == "prolif_count":
            try:
                numpy_fp = get_prolif_from_biomolecule(ligand, ligand, format=format, interactions=interactions)
            except Exception as e:
                tqdm.write(f"Skipping ligand {getattr(ligand, 'entity_id', 'unknown')}: {e}")
                continue
            fingerprint_list.append(numpy_fp)  # grabbing the first fingerprint (pose) this is the usual case
            ligand_ids.append(ligand.entity_id)
        else:
            raise ValueError(f"Unknown fingerprint type: {fingerprint_type}. Use 'ecfp4', 'ecfp6', or 'prolif'.")
        
        elapsed_time = time.time() - start_time
        avg_time_per_ligand = elapsed_time / (ligands.index(ligand) + 1)
        remaining_time = avg_time_per_ligand * (len(ligands) - (ligands.index(ligand) + 1))
        tqdm.write(f"Estimated time remaining: {remaining_time:.2f} seconds")
    
    # Combine all fingerprints into a single DataFrame
    fingerprint_df = pd.Series(fingerprint_list, index=ligand_ids)
    return fingerprint_df

#pushes calculations back to the biomolecule object chain
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

    for fingerprint_df_id in fingerprint_df.index:
        for ligand in biomolecules:
            if ligand.entity_id == fingerprint_df_id:
                print("fingerprint to be stored: ")
                print(fingerprint_df.loc[fingerprint_df_id])
                print("fingerprint type: ")
                print(fingerprint_type)
                value = fingerprint_df.loc[fingerprint_df_id]

                # If it's an RDKit ExplicitBitVect (or duck-typed), convert to numpy array
                if isinstance(value, DataStructs.ExplicitBitVect) or hasattr(value, "GetNumBits"):
                    try:
                        n_bits = value.GetNumBits()
                        arr = np.zeros((n_bits,), dtype=np.uint8)
                        DataStructs.ConvertToNumpyArray(value, arr)
                        ligand.store_fingerprint(arr, fingerprint_type)
                    except Exception:
                        # Fallback: try list conversion
                        ligand.store_fingerprint(np.asarray(list(value)), fingerprint_type)
                elif isinstance(value, np.ndarray):
                    ligand.store_fingerprint(value, fingerprint_type)
                else:
                    # For pandas Series, lists, or other array-like objects
                    ligand.store_fingerprint(np.asarray(value), fingerprint_type)
                break

def calculate_fingerprints(biomolecules: list[TriageBiomolecule], fingerprint_type: str) -> None:
    """
    Calculate and store fingerprints for a list of TriageBiomolecule objects.
    
    Parameters
    ----------
    biomolecules : list of TriageBiomolecule
        List of TriageBiomolecule objects representing ligands.
    fingerprint_type : str
        Type of fingerprint to calculate: 'ecfp4', 'ecfp6', or 'prolif'.
    
    Notes
    -----
    This function calculates fingerprints for each biomolecule and stores them in the biomolecule object.
    """
    
   # if not receptor:
    #    raise ValueError("No receptor founds in the provided biomolecules.")
    
    fingerprint_df = create_fingerprint_dataframe(biomolecules, fingerprint_type=fingerprint_type)
    store_fingerprints_from_df(fingerprint_df, biomolecules, fingerprint_type=fingerprint_type)



def generate_prolif_fingerprints_from_outputs(
    outputs_dir: str,
    receptor_pdb_path: str,
    output_parquet: str,
    pose_ext: str = ".cif",
    verbose: bool = False):
    """
    Traverse the outputs directory, in the outputs directory there are subdirectories 
    with the names of compounds. in each compound subdirectory there is a .cif file that ends in the
    format {compound_name}_model_0.cif. With that .cif file as input use the extract_ligand() and
    extract_protein_receptor() to genereate parallel files iwth the suffix {compound_name}_model_ligand.cif and {compound_name}_model_protein.cif,
    then use openbabel to convert the ligand.cif file to .sdf format. With these new files, generated, create 
    a TriageBiomolecule object out of it and put that into a list. From that list of TriageBiomolecule objects, calculate prolif
    fingerprints from them using the calculate_fingerprints() funciton and then save the generated data as a parquet.
    and export as a parquet file.
    """
    biomolecules = []
    for compound_name in os.listdir(outputs_dir):
        compound_dir = os.path.join(outputs_dir, compound_name)
        if not os.path.isdir(compound_dir):
            continue
        # Find the .cif file ending with _model_0.cif
        cif_files = [f for f in os.listdir(compound_dir) if f.endswith("_model_0.cif")]
        if not cif_files:
            if verbose:
                print(f"No _model_0.cif file found in {compound_dir}")
            continue
        cif_file = cif_files[0]
        cif_path = os.path.join(compound_dir, cif_file)
        ligand_cif_path = os.path.splitext(cif_path)[0] + "_ligand.cif"
        protein_cif_path = os.path.splitext(cif_path)[0] + "_protein.cif"
        protein_pdb_path = os.path.splitext(cif_path)[0] + "_protien.pdb"
        ligand_sdf_path = os.path.splitext(ligand_cif_path)[0] + ".sdf"
        # Extract ligand and protein
        try:
            extract_ligand(cif_path, ligand_cif_path)
            extract_protein_receptor(cif_path, protein_pdb_path, out_format="pdb")
        except Exception as e:
            if verbose:
                print(f"Failed to extract ligand/protein from {cif_path}: {e}")
            continue
        # Convert ligand cif to sdf using openbabel
        try:
            mols = list(pybel.readfile("cif", ligand_cif_path))
            if not mols:
                if verbose:
                    print(f"No molecules found in {ligand_cif_path}")
            else:
                out = pybel.Outputfile("sdf", ligand_sdf_path, overwrite=True)
                for m in mols:
                    out.write(m)
                out.close()
        except Exception as e:
            if verbose:
                print(f"OpenBabel conversion failed for {ligand_cif_path}: {e}")
            continue
        # Create TriageBiomolecule object
        biomol = TriageBiomolecule(
            entity_id=compound_name,
            entity_type="ligand",
            pose_path=ligand_sdf_path,
            paired_receptor_path=protein_pdb_path
        )
        biomolecules.append(biomol)
    if verbose:
        print(f"Calculating prolif fingerprints for {len(biomolecules)} biomolecules...")
    if verbose:
        print(f"Exporting fingerprints to {output_parquet}...")
           
    # Calculate prolif fingerprints and export as parquet
    calculate_fingerprints(biomolecules, fingerprint_type="prolif_count")
    TriageBiomolecule.to_parquet(biomolecules, output_parquet)
    