"""
PLIP Batch Analysis for Boltz Predictions with Replicate Structure

Handles directory structure:
    base_dir/
    ├── tgcpl_cov_1/
    │   ├── compound_001_XXXX/
    │   │   └── boltz_results_XXX/predictions/XXX/PROT_LIG_model_0.cif
    │   └── compound_002_XXXX/
    └── tgcpl_cov_2/
        └── ...

Finds unique ligands across all replicates and performs PLIP analysis once per unique ligand.
"""

from plip.basic import config
from plip.structure.preparation import PDBComplex
import subprocess
import argparse
import fnmatch
import shutil
import re
import gemmi
import sys
import os
import tempfile
import csv
from collections import defaultdict


def get_interactions(interactions):
    """Takes a PLIP interaction object and gets the number of each interaction."""
    num_saltbridges = len(interactions.saltbridge_lneg + interactions.saltbridge_pneg)
    num_hbonds = len(interactions.hbonds_ldon + interactions.hbonds_pdon)
    num_pication = len(interactions.pication_laro + interactions.pication_paro)
    num_pistack = len(interactions.pistacking)
    num_halogen = len(interactions.halogen_bonds)
    num_waterbridges = len(interactions.water_bridges)

    interactions_vals = [num_saltbridges, num_hbonds, num_pication, num_pistack, num_halogen, num_waterbridges]
    return interactions_vals


def get_interacting_residues(interactions):
    """Takes a PLIP interaction object and gets the resids for residues in all interactions."""
    interacting_residues = {}

    for saltbridge in interactions.saltbridge_lneg + interactions.saltbridge_pneg:
        interacting_residues[f'{saltbridge.restype}{saltbridge.resnr}'] = "saltbridge"

    for hbond in interactions.hbonds_ldon + interactions.hbonds_pdon:
        interacting_residues[f'{hbond.restype}{hbond.resnr}'] = "hbond"

    for pication in interactions.pication_laro + interactions.pication_paro:
        interacting_residues[f'{pication.restype}{pication.resnr}'] = "pication"

    for pistack in interactions.pistacking:
        interacting_residues[f'{pistack.restype}{pistack.resnr}'] = "pistack"

    for halogen in interactions.halogen_bonds:
        interacting_residues[f'{halogen.restype}{halogen.resnr}'] = "halogen"

    for waterbridge in interactions.water_bridges:
        interacting_residues[f'{waterbridge.restype}{waterbridge.resnr}'] = "waterbridge"

    return interacting_residues

def find_all_cif_files_in_replicates(root_dir,protein_name, VERBOSE):
    """
    Finds all CIF files across all replicate directories.
    
    Structure expected:
        root_dir/
        ├── tgcpl_cov_1/
        │   ├── PROT_001_model_0.cif
        │   ├── PROT_002_model_0.cif
        │   └── ...
        ├── tgcpl_cov_2/
        │   ├── PROT_001_model_0.cif
        │   └── ...
        └── tgcpl_cov_3/
            └── ...
    
    :param root_dir: Base directory containing replicate subdirectories
    :return: Dictionary mapping substance_id -> list of CIF file paths
    """
    print(f"🔍 Searching for CIF files in {root_dir}...")
    
    compound_files = defaultdict(list)
    
    # Look for subdirectories (replicates)
    if not os.path.isdir(root_dir):
        print(f"Directory not found: {root_dir}")
        return compound_files
    
    # Get all subdirectories
    subdirs = [d for d in os.listdir(root_dir) 
               if os.path.isdir(os.path.join(root_dir, d))]
    
    if not subdirs:
        print(f"No subdirectories found in {root_dir}")
        return compound_files
    
    print(f"Found {len(subdirs)} replicate directories: {', '.join(subdirs[:3])}{'...' if len(subdirs) > 3 else ''}")
    
    # Search each replicate directory
    for subdir in subdirs:
        if protein_name in subdir: 
            subdir_path = os.path.join(root_dir, subdir)
            subdir_cifs = os.path.join(subdir_path,'cifs')
            
            # Find all *.cif files directly in this directory
            cif_files = [f for f in os.listdir(subdir_cifs) if f.endswith('.cif')]
            
            if VERBOSE: (f"   {subdir}: {len(cif_files)} CIF files", False)
            
            for cif_file in cif_files:
                full_path = os.path.join(subdir_cifs, cif_file)
                
                # Extract substance_id from filename
                # Format: PROT_COMPOUND_model_0.cif or similar
                match = re.search(r'_([A-Za-z0-9\-]+)\.cif$', cif_file)
                if match:
                    substance_id = match.group(1)
                
                compound_files[substance_id].append(full_path)
        
    # Print summary
    total_files = sum(len(files) for files in compound_files.values())
    print(f"Found {total_files} CIF files for {len(compound_files)} unique compounds")
    
    if compound_files:
        # Show replicate distribution
        rep_counts = [len(files) for files in compound_files.values()]
        print(f"Replicates per compound: min={min(rep_counts)}, max={max(rep_counts)}, mean={sum(rep_counts)/len(rep_counts):.1f}")
    
    return compound_files


def select_representative_structure(cif_files, selection_method='first'):
    """
    Select one representative structure from multiple replicates.
    
    :param cif_files: List of CIF file paths for same compound
    :param selection_method: Method to select representative ('first', 'random', 'best_confidence')
    :return: Path to selected CIF file
    """
    if len(cif_files) == 1:
        return cif_files[0]
    
    if selection_method == 'first':
        # Use alphabetically first (most consistent)
        return sorted(cif_files)[0]
    elif selection_method == 'random':
        import random
        return random.choice(cif_files)
    elif selection_method == 'best_confidence':
        # Would need to parse confidence from JSON files
        # For now, fall back to first
        return sorted(cif_files)[0]
    else:
        return cif_files[0]


def extract_compound_id_from_path(filepath):
    """
    Extract compound_id from filepath.
    
    Expected patterns:
    - /path/to/PROT_COMPOUND_model_0.cif
    - /path/to/compound_XXX/...
    """
    # Try filename first
    filename = os.path.basename(filepath)
    match = re.search(r'([A-Za-z0-9_]+)_model_0\.cif$', filename)
    if match:
        base_name = match.group(1)
        parts = base_name.split('_')
        if len(parts) >= 2:
            return parts[-1]
        return base_name
    
    # Try directory structure
    path_parts = filepath.split(os.sep)
    for part in reversed(path_parts):
        if part.startswith('compound_'):
            return part.replace('compound_', '').split('_')[0]
    
    return None


def extract_ligand_code(filename):
    """
    Extracts ligand code from filename for PLIP interaction set lookup.
    
    From: PROT_LIG_model_0.cif
    Returns: First 3 chars of LIG
    """
    if "/" in filename:
        mod_filename = filename.split("/")[-1]
    else:
        mod_filename = filename
    
    match = re.search(r'_([A-Za-z0-9]+)_model_0\.cif$', mod_filename)
    if match:
        lig_full = match.group(1)
        return lig_full[:3]  # Take first 3 characters
    else:
        # Fallback
        return "UNK"


def convert_cif_to_pdb(cif_path, pdb_path, chain_map=None):
    """Convert mmCIF to PDB using gemmi, fixing long chain names."""
    doc = gemmi.cif.read_file(cif_path)
    block = doc.sole_block()
    structure = gemmi.make_structure_from_block(block)

    # Fix chain names longer than 1 character
    for model in structure:
        for chain in model:
            if len(chain.name) != 1:
                new_chain = chain_map.get(chain.name, 'X') if chain_map else 'X'
                chain.name = new_chain

    structure.write_pdb(pdb_path)


def convert_pdb_to_pdb(source_pdb_path, target_pdb_path):
    """Standardize or reformat a PDB file."""
    structure = gemmi.read_structure(source_pdb_path)
    structure.write_pdb(target_pdb_path)


def main(args):
    headers = ["compound_id", "name", "smiles", "inchi", "molwt", 
               "numheavy", "numrotbonds", "numrings", 
               "hydrophobicatoms", "hbondacceptors", 
               "saltbridges", "hbonds", "pication", 
               "pistack", "halogen", "waterbridge"]
    
    res_headers = ["compound_id", "name", "residue", "interaction_type"]

    temp_dir = tempfile.mkdtemp()
    VERBOSE = args.verbose
    # Set config based on receptor type
    if args.receptor_type == "protein":
        config.DNARECEPTOR = False
    elif args.receptor_type in ["rna", "dna"]:
        config.DNARECEPTOR = True 
    else:
        print("Incorrect receptor type! Use: DNA, RNA, or PROTEIN")
        sys.exit(1)

    # Find all CIF files across replicates
    compound_files = find_all_cif_files_in_replicates(args.directory, args.protein_name, VERBOSE)
    
    if not compound_files:
        print("❌ No CIF files found in directory structure!")
        sys.exit(1)
    
    print(f"\nAnalyzing {len(compound_files)} unique compounds...")
    
    collected_data = []
    residue_data = []
    errors = []
    
    for compound_id, cif_files in sorted(compound_files.items()):
        # Select one representative structure per compound
        target_file = select_representative_structure(cif_files, args.selection_method)
        
        if VERBOSE: print(f"\nProcessing compound {compound_id} ({len(cif_files)} replicates found)")
        if VERBOSE: print(f"  Using: {target_file}")
        
        # Get name from filepath
        name = os.path.splitext(os.path.basename(target_file))[0]
        
        # Extract ligand code for PLIP
        lig_code = extract_ligand_code(target_file)
        
        # Convert to PDB
        filetype = target_file.split(".")[-1].strip()
        target_pdb = os.path.join(temp_dir, f"{compound_id}.pdb")
        
        try:
            if filetype == "cif":
                convert_cif_to_pdb(target_file, target_pdb)
            elif filetype == "pdb":
                convert_pdb_to_pdb(target_file, target_pdb)
            else:
                print(f"Unsupported file type: {filetype}")
                errors.append(target_file)
                continue

            # Run PLIP analysis
            my_mol = PDBComplex()
            my_mol.load_pdb(target_pdb)
            my_mol.analyze()
            
            # Try to find interaction set
            # PLIP uses format: "LIG:CHAIN:RESNUM"
            interaction_key = f"{lig_code}:X:1"
            
            if interaction_key not in my_mol.interaction_sets:
                # Try to find any interaction set
                if my_mol.interaction_sets:
                    interaction_key = list(my_mol.interaction_sets.keys())[0]
                    if VERBOSE: print(f"  Using interaction set: {interaction_key}")
                else:
                    print(f"⚠️  No interactions found for {compound_id}")
                    errors.append(target_file)
                    continue
            
            interactions = my_mol.interaction_sets[interaction_key]
            interactions_vals = get_interactions(interactions)
            
            # Organize data for CSV
            plip_fprint = [
                compound_id,
                name,
                interactions.ligand.smiles.strip() if interactions.ligand.smiles else "",
                interactions.ligand.inchikey.strip() if interactions.ligand.inchikey else "",
                interactions.ligand.molweight,
                interactions.ligand.heavy_atoms,
                interactions.ligand.num_rot_bonds,
                interactions.ligand.num_rings,
                len(interactions.ligand.hydroph_atoms),
                len(interactions.ligand.hbond_acc_atoms)
            ]
            
            for val in interactions_vals:
                plip_fprint.append(val)
            
            collected_data.append(plip_fprint)

            # Collect residue interactions
            residue_interactions = get_interacting_residues(interactions)
            for resid, interaction_type in residue_interactions.items():
                res_fprint = [compound_id, name, resid, interaction_type]
                residue_data.append(res_fprint)
            
            if VERBOSE: print(f"✓ {compound_id}")
            
        except Exception as e:
            print(f"❌ Error processing {compound_id}: {e}")
            errors.append(target_file)
            continue
    
    # Write results
    print(f"\n💾 Writing results...")
    
    csv_out = os.path.join(args.outdir, f"{args.csv_name}.csv")
    if VERBOSE: print(f"Writing fingerprints to {csv_out}")
    with open(csv_out, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(headers)
        writer.writerows(collected_data)
    
    residues_csv_out = os.path.join(args.outdir, f"{args.csv_name}_residues.csv")
    if VERBOSE: print(f"Writing residue interactions to {residues_csv_out}")
    with open(residues_csv_out, 'w', newline='') as res_csv:
        writer = csv.writer(res_csv)
        writer.writerow(res_headers)
        writer.writerows(residue_data)

    # Clean up
    shutil.rmtree(temp_dir)
    
    print(f"\n✅ Analysis complete!")
    print(f"   Processed: {len(collected_data)} compounds")
    print(f"   Failed: {len(errors)} structures")
    print(f"   Output: {csv_out}")

    return errors


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Analyze Boltz predictions with PLIP across replicate directories.\n"
                   "Finds unique compounds and analyzes one representative structure per compound.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument("-d", "--directory", required=True,
                       help="Root directory containing replicate subdirectories")
    parser.add_argument("-o", "--outdir", required=True,
                       help="Output directory for results")
    parser.add_argument("-rt", "--receptor_type", required=True,
                       help="Receptor type: 'dna', 'rna', or 'protein'", 
                       choices=['dna', 'rna', 'protein'])
    parser.add_argument("--verbose", "-v", action="store_true",
                       help="Enable verbose output")
    parser.add_argument("-n", "--csv_name", required=False, default="plip_fingerprints",
                       help="Base name for CSV output files (no extension)")
    parser.add_argument("-p", "--protein_name", required=True,
                       help="Protein to run PLIP on.")
    parser.add_argument("--selection_method", default="first",
                       choices=['first', 'random', 'best_confidence'],
                       help="Method to select representative from replicates (default: first)")

    args = parser.parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    
    main(args)