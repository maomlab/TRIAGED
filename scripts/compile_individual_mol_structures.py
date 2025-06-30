import argparse
import os
import glob
import tempfile
import gemmi
import mdtraj as md
import re

def parse_args():
    parser = argparse.ArgumentParser(description="Align substructures for a single " \
    "Boltz2 job to parent and save a multi-model PDB. The DIRECTORY should point " \
    "to the output 'predictions' directory. Example usage " \
    "for PDB 3JQZ with all models in the 'boltz_jobs' directory:" \
    "python compile_best_model_structures -p 3JQZ.pdb -d boltz_jobs/mol1/boltz_results_mol1/predictions/mol1"\
    " -o 3JQZ_top_models.pdb")
    parser.add_argument("-p", "--parent", 
                        help="Path to the parent file to align to. Can be CIF" \
                        " or PDB. If no structure is given, will align to " \
                        "first structure in DIRECTORY.")
    
    parser.add_argument("-d", "--directory", required=True, 
                        help="Directory of CIF substructures from Boltz2.")
    parser.add_argument("-o", "--output", default="aligned_models.pdb", 
                        help="Name for the output multi-model PDB file. Defaults to 'aligned_models.pdb'")
    parser.add_argument("--max-models", type=int, 
                        help="OPTIONAL: Maximum number of models to include in"
                        "the output PDB. Defaults to all in the directory.")
    parser.add_argument("--verbose", "-v", action="store_true", 
                        help="Enable verbose output.")
    parser.add_argument("-ep", "--exclude_parent", action='store_true', default=False,
                        help="OPTIONAL: If enabled, will exclude the parent trajectory" \
                        "from the final output aligned PDB.")
    return parser.parse_args()

def verbose_print(msg, verbose):
    """Prints message if verbose is True."""
    if verbose:
        print(msg)


def convert_pdb_to_pdb(source_pdb_path, target_pdb_path):
    """Standardize or reformat a PDB file."""
    # Read the structure from the source PDB
    structure = gemmi.read_structure(source_pdb_path)
    
    # Write the structure to the target PDB
    structure.write_pdb(target_pdb_path)


def convert_cif_to_pdb(cif_path, pdb_path):
    """Convert mmCIF to PDB using gemmi."""
    doc = gemmi.cif.read_file(cif_path)
    block = doc.sole_block()
    structure = gemmi.make_structure_from_block(block)
    structure.write_pdb(pdb_path)

def get_matching_atoms(ref, traj):
    """Find matching atoms between parent and substructure by (residue index, atom name)."""
    ref_atoms = {(a.residue.index, a.name): i for i, a in enumerate(ref.topology.atoms)}
    traj_atoms = {(a.residue.index, a.name): i for i, a in enumerate(traj.topology.atoms)}

    common_keys = sorted(set(ref_atoms) & set(traj_atoms))
    if not common_keys:
        raise ValueError("No matching atoms found between structures.")

    ref_idx = [ref_atoms[k] for k in common_keys]
    traj_idx = [traj_atoms[k] for k in common_keys]
    return ref_idx, traj_idx

def extract_key(filename):
    """Extracts the numeric key from a filename of the form NAME_KEY.cif."""
    match = re.search(r'_(\d+)\.cif$', filename)
    if match:
        return int(match.group(1))
    else:
        raise ValueError(f"Filename {filename} does not match the expected pattern.")

def save_with_model_names(trajs, filenames, output_path):
    """Write out multi-model PDB with MODEL/ENDMDL blocks and filename REMARKs for each."""
    with open(output_path, "w") as out:
        for i, (traj, filename) in enumerate(zip(trajs, filenames), 1):
            # Extract the base name without the directory and the .cif extension
            model_name = os.path.splitext(os.path.basename(filename))[0]
            
            # Create a temporary file to save the trajectory as PDB
            with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as temp_pdb:
                traj.save_pdb(temp_pdb.name)
                temp_pdb_path = temp_pdb.name

            # Read the content back from the temporary file
            with open(temp_pdb_path, "r") as pdb_file:
                pdb_str = pdb_file.read()

            out.write(f"MODEL     {i:>4}\n")
            out.write(f"TITLE     {model_name}\n")
            out.write(f"REMARK   1 Model derived from {model_name}\n")
            for line in pdb_str.splitlines():
                if not line.startswith("MODEL") and not line.startswith("ENDMDL"):
                    out.write(line + "\n")
            out.write("ENDMDL\n")

            # Clean up the temporary file
            os.remove(temp_pdb_path)

def load_structure(filename):
    """Loads a structure file which can be CIF or PDB"""
    extension = os.path.splitext(filename)[1].lower()
    if extension == '.cif':
        temp_pdb = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
        convert_cif_to_pdb(filename, temp_pdb)
        return md.load(temp_pdb)
    elif extension == '.pdb':
        return md.load(filename)
    else:
        raise ValueError(f"Unsupported file format for {filename}")

def main():
    args = parse_args()
    temp_dir = tempfile.mkdtemp()

    verbose_print(f"Verbose mode activated. Temporary directory is {temp_dir}", 
                  args.verbose)

    # Find the parent structure
    parent_file = args.parent
    cif_files = glob.glob(os.path.join(args.directory, "*.cif"))
    cif_files.sort(key=extract_key)

    # Use first CIF as parent if not provided
    if not parent_file:
        if cif_files:
            parent_file = cif_files[0]
            verbose_print(f"No parent provided. Using {parent_file} as the parent structure.", args.verbose)
        else:
            print("❌ No CIF files available to process.")
            return

    # Load parent structure
    parent_traj = load_structure(parent_file)
    if (args.exclude_parent == False):
        parent_pdb = os.path.join(temp_dir, "parent.pdb")
        if parent_file.split(".")[-1].strip() == "cif":
            convert_cif_to_pdb(parent_file, parent_pdb)
            
        elif parent_file.split(".")[-1].strip() == "pdb":
            parent_pdb = os.path.join(temp_dir, "parent.pdb")
            convert_pdb_to_pdb(parent_file, parent_pdb)

        else:
            print("Unsupported parent file type! Must be PDB or CIF!")
            exit()

        ptraj = md.load(parent_pdb)
        
        aligned_trajectories = [ptraj]
        model_names = ["parent"]
    else:
        aligned_trajectories = []
        model_names = []


    # Loop through substructure CIFs
    print(f"Found {len(cif_files)} substructure CIFs.")

    for idx, cif_path in enumerate(cif_files):
        if args.max_models is not None and idx >= args.max_models:
            break  # Stop processing more files if the maximum limit is reached

        try:
            sub_pdb = os.path.join(temp_dir, f"{os.path.splitext(os.path.basename(cif_path))[0]}.pdb")

            convert_cif_to_pdb(cif_path, sub_pdb)
            traj = md.load(sub_pdb)
            verbose_print(f"Loaded substructure {cif_path} and converted to PDB.", args.verbose)

            # Match and align
            ref_idx, traj_idx = get_matching_atoms(parent_traj, traj)
            traj.superpose(parent_traj, atom_indices=traj_idx, ref_atom_indices=ref_idx)

            aligned_trajectories.append(traj)
            model_names.append(cif_path)  # Track the filename for each trajectory

            print(f"Aligned: {cif_path}")

        except Exception as e:
            print(f"❌ Skipping {cif_path} due to error: {e}")
            verbose_print(f"Error with {cif_path}: {e}", args.verbose)

    if aligned_trajectories:
        save_with_model_names(aligned_trajectories, model_names, args.output)
        print(f"\n✅ Saved aligned multi-model PDB: {args.output}")
    else:
        print("⚠️ No models were successfully aligned.")

if __name__ == "__main__":
    main()
