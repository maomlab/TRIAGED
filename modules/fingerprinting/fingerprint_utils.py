#from ..triage_biomolecule.triage_biomolecule import TriageBiomolecule
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from triage_biomolecule.triage_biomolecule import TriageBiomolecule
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


# Try to import Biopython for mmCIF parsing; fall back to pybel/openbabel if needed
_HAS_BIOPY = True
try:
    from Bio.PDB import MMCIFParser, PDBIO, Select
except Exception:
    _HAS_BIOPY = False

_HAS_PYBEL = True
try:
    from openbabel import openbabel
    import pybel
except Exception:
    _HAS_PYBEL = False

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
 
#creates fingerprints
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
                ligand.store_fingerprint(fingerprint_df.loc[fingerprint_df_id].values, fingerprint_type)
                break

def calculate_fingerprints(biomolecules: list[TriageBiomolecule], receptor: TriageBiomolecule, fingerprint_type: str) -> None:
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
    
    if not receptor:
        raise ValueError("No receptor found in the provided biomolecules.")
    
    fingerprint_df = create_fingerprint_dataframe(biomolecules, receptor, fingerprint_type=fingerprint_type)
    store_fingerprints_from_df(fingerprint_df, biomolecules, fingerprint_type=fingerprint_type)

def extract_ligands_to_sdf(
    cif_path: str,
    out_sdf: str,
    model_index: int = 0,
    exclude_resnames: Iterable[str] = ("HOH", "WAT", "DOD", "H2O"),
    exclude_common_ions: bool = True,
    min_atoms: int = 3,
) -> int:
    """
    Extract hetero small-molecule ligands from an mmCIF and save them to an SDF.
    
    Parameters
    ----------
    cif_path : str
        Path to input .cif (or .cif.gz).
    out_sdf : str
        Path to output .sdf (multi-molecule).
    model_index : int
        Which model to use (0 = first).
    exclude_resnames : Iterable[str]
        Residue names to ignore (waters by default).
    exclude_common_ions : bool
        If True, skip common monoatomic/diatomic ions (NA, CL, MG, CA, ZN, ...).
    min_atoms : int
        Skip residues with fewer than this many atoms (avoids single ions).
    
    Returns
    -------
    int
        Number of ligand molecules written.
    """
    # --- Parse structure
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("X", cif_path)
    models = list(structure.get_models())
    if not models:
        raise ValueError("No models found in structure.")
    if model_index >= len(models):
        raise IndexError(f"Requested model_index {model_index} but only {len(models)} model(s) present.")
    model = models[model_index]

    # Sets for filtering
    exclude_resnames = set(exclude_resnames)
    common_ions: Set[str] = {
        "NA", "K", "CL", "CA", "MG", "MN", "ZN", "FE", "CU", "CO", "NI", "CD", "SR",
        "CS", "BR", "F", "I", "HG", "PB", "BA", "AL", "RB", "YB", "PT", "AG", "AU"
    }

    def is_hetero_residue(res) -> bool:
        # res.id tuple: (hetfield, resseq, icode); hetfield like "H_" for HETATM, " " for ATOM
        hetflag = res.id[0]
        return hetflag.strip() != ""  # True for HETATM

    def is_water_or_excluded(res) -> bool:
        rn = res.get_resname().strip().upper()
        if rn in exclude_resnames:
            return True
        if exclude_common_ions and (rn in common_ions or len(rn) <= 2):
            return True
        return False

    def choose_altloc(atoms):
        """
        Yield atoms with a consistent altloc choice per atom name:
        Prefer blank altloc '', else 'A', else first seen.
        """
        chosen = {}
        for a in atoms:
            alt = a.get_altloc() or ""
            key = (a.get_id(), a.get_name())
            if key not in chosen:
                chosen[key] = a
            else:
                # prefer '' over others; then 'A'
                curr = chosen[key]
                curr_alt = curr.get_altloc() or ""
                if curr_alt != "" and (alt == "" or (curr_alt != "A" and alt == "A")):
                    chosen[key] = a
        return list(chosen.values())

    # --- Collect hetero residues (ligands)
    lig_residues = []
    for chain in model.get_chains():
        for res in chain.get_residues():
            if not is_hetero_residue(res):
                continue
            if is_water_or_excluded(res):
                continue
            # Require a minimum atom count to avoid ions
            if sum(1 for _ in res.get_atoms()) < min_atoms:
                continue
            lig_residues.append((chain.id, res))

    if not lig_residues:
        # Create an empty SDF to be explicit, or raise?
        open(out_sdf, "w").close()
        return 0

    writer = SDWriter(out_sdf)
    n_written = 0

    for chain_id, res in lig_residues:
        atoms = list(res.get_atoms())
        atoms = choose_altloc(atoms)

        # Build an RDKit molecule (no bonds yet)
        rw = rdchem.RWMol()
        conf = rdchem.Conformer(len(atoms))

        idx_map = {}  # Biopython atom -> RDKit index
        for i, a in enumerate(atoms):
            # element symbol can be None in some mmCIFs; fall back to first letter of atom name
            elem = (a.element or a.get_name().strip()[0]).capitalize()
            try:
                rd_atom = rdchem.Atom(elem)
            except RuntimeError:
                # Fallback to carbon if element fails (very rare labels)
                rd_atom = rdchem.Atom("C")
            rd_idx = rw.AddAtom(rd_atom)
            x, y, z = a.coord  # numpy array
            conf.SetAtomPosition(rd_idx, Point3D(float(x), float(y), float(z)))
            idx_map[i] = rd_idx

        rw.AddConformer(conf, assignId=True)

        # Infer bonds from 3D geometry
        mol = rw.GetMol()
        try:
            rdDetermineBonds.DetermineBonds(mol)  # heuristic bond & order perception
            Chem.SanitizeMol(mol)
        except Exception:
            # As a fallback, try connectivity-only (no bond orders)
            try:
                rdDetermineBonds.DetermineConnectivity(mol)
                Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ADJUSTHS)
            except Exception:
                # If it still fails, skip this residue
                continue

        # Title & annotations
        resname = res.get_resname().strip()
        resseq = res.id[1]
        icode = res.id[2].strip() if res.id[2] else ""
        title = f"{resname}_{chain_id}_{resseq}{icode}"
        mol.SetProp("_Name", title)
        mol.SetProp("resname", resname)
        mol.SetProp("chain_id", chain_id)
        mol.SetProp("resseq", str(resseq))
        mol.SetProp("icode", icode)

        writer.write(mol)
        n_written += 1

    writer.close()
    return n_written

def generate_prolif_fingerprints_from_outputs(
    outputs_dir: str,
    receptor_pdb_path: str,
    output_parquet: str,
    pose_ext: str = ".cif",
    verbose: bool = False):
    """
    Traverse the outputs directory, generate ProLIF fingerprints for each ligand-receptor pair,
    and export as a parquet file.

    Args:
        outputs_dir (str): Path to the outputs directory.
        receptor_pdb_dir (str): Path to the directory containing receptor PDB files named as {receptor_name}.pdb.
        output_parquet (str): Path to save the resulting parquet file.
        pose_ext (str): Extension of the ligand pose files (default: ".cif").
    """

    # Determine which subdirectory to use
    combined_dir = os.path.join(outputs_dir, "combined_predictions")
    predictions_dir = os.path.join(outputs_dir, "predictions")
    print(combined_dir)
    print(predictions_dir)
    if os.path.isdir(combined_dir):
        preds_dir = combined_dir
    elif os.path.isdir(predictions_dir):
        preds_dir = predictions_dir
    else:
        raise FileNotFoundError("Neither 'combined_predictions' nor 'predictions' subdirectory found in outputs_dir.")
    
    # Walk through every subdirectory in preds_dir
    for root, dirs, files in os.walk(preds_dir):
        for subdir in dirs:
            subdir_path = os.path.join(root, subdir)
            # Now you can access the contents of each subdirectory
            subdir_files = os.listdir(subdir_path)
            if verbose:
                print(f"\033[92mSubdirectory: {subdir_path}\033[0m")
                print(f"\033[92mFiles: {subdir_files}\033[0m")
            # You can process each file in subdir_files as needed
            cif_files = [f for f in subdir_files if f.endswith(pose_ext)]
            if not cif_files:
                if verbose:
                    print(f"No .cif file found in {subdir_path}")
                continue
            cif_path = os.path.join(subdir_path, cif_files[0])
            if verbose:
                print(f"Found .cif file: {cif_path}")
            #now creating a sdf file based on the cif file
            sdf_path = os.path.splitext(cif_path)[0] + ".sdf"
            if verbose:
                print(f"Creating .sdf file: {sdf_path}")

            extract_ligands_to_sdf(cif_path, sdf_path)