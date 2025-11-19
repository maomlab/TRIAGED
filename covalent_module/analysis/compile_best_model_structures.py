import os
import tempfile
import argparse
import re
import gemmi
import string
import fnmatch
import mdtraj as md

def verbose_print(msg, verbose):
    """Prints message if verbose is True."""
    if verbose:
        print(msg)

def convert_cif_to_pdb(cif_path, pdb_path):
    """Convert mmCIF to PDB using gemmi."""
    doc = gemmi.cif.read_file(cif_path)
    block = doc.sole_block()
    structure = gemmi.make_structure_from_block(block)
    # structure.write_pdb(pdb_path)

    # rename chains if longer than 1 character
    used_ids = set()
    available_ids = list(string.ascii_uppercase + string.digits)
    
    for chain in structure[0]:  # iterate over chains in the first model
        if len(chain.name) > 1 or chain.name in used_ids:
            chain.name = available_ids[0] 
        available_ids.remove(chain.name)
        used_ids.add(chain.name)
    
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
            model_name = os.path.splitext(os.path.basename(filename))[0]
            
            with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as temp_pdb:
                traj.save_pdb(temp_pdb.name)
                temp_pdb_path = temp_pdb.name

            with open(temp_pdb_path, "r") as pdb_file:
                pdb_str = pdb_file.read()

            out.write(f"MODEL     {i:>4}\n")
            out.write(f"TITLE     {model_name}\n")
            out.write(f"REMARK   1 Model derived from {model_name}\n")
            for line in pdb_str.splitlines():
                if not line.startswith("MODEL") and not line.startswith("ENDMDL"):
                    out.write(line + "\n")
            out.write("ENDMDL\n")

            os.remove(temp_pdb_path)

def load_structure(filename):
    """Loads a structure file which can be CIF, PDB, or MOL2."""
    extension = os.path.splitext(filename)[1].lower()
    if extension == '.cif':
        temp_pdb = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
        convert_cif_to_pdb(filename, temp_pdb)
        return md.load(temp_pdb)
    elif extension == '.pdb':
        return md.load(filename)
    elif extension == '.mol2':
        return md.load(filename)
    else:
        raise ValueError(f"Unsupported file format for {filename}")

def find_first_model_files(root_dir):
    '''
    Finds the first model CIF files in 'predictions' subdirectories.
    :param root_dir: The root directory to start the search from.
        root_dir must follow this structure per ligand: 
        <root_dir>/.../predictions/../model_output_0.cif

    :return: A list of paths to the first model CIF files.
    '''
    prediction_files = []
    for dirpath, dirnames, filenames in os.walk(root_dir):

        if 'predictions/' in dirpath and fnmatch.filter(filenames, '*_0.cif'):
            prediction_files.append(os.path.join(dirpath, fnmatch.filter(filenames, '*_0.cif')[0]))

    return prediction_files

def main(args):
    temp_dir = tempfile.mkdtemp()

    verbose_print(f"Verbose mode activated. Temporary directory is {temp_dir}", args.verbose)

    # Find and collect first model _0.cif files
    collected_files = find_first_model_files(args.directory)

    # Determine parent structure
    parent_file = args.parent
    if not parent_file:
        if collected_files:
            parent_file = collected_files[0]
            verbose_print(f"No parent provided. Using {parent_file} as the parent structure.", args.verbose)
        else:
            print("❌ No CIF files available in 'predictions' subdirectories.")
            return

    # Load parent structure
    parent_traj = load_structure(parent_file)

    # Align and save models
    aligned_trajectories = []
    model_names = []
    
    for idx, cif_path in enumerate(collected_files):
        # if args.max_models is not None and idx >= args.max_models:
        #     break

        try:
            temp_pdb = os.path.join(temp_dir, f"{os.path.splitext(os.path.basename(cif_path))[0]}.pdb")
            convert_cif_to_pdb(cif_path, temp_pdb)
            traj = md.load(temp_pdb)
            verbose_print(f"Loaded {cif_path} and converted to PDB.", args.verbose)

            ref_idx, traj_idx = get_matching_atoms(parent_traj, traj)
            traj.superpose(parent_traj, atom_indices=traj_idx, ref_atom_indices=ref_idx)

            aligned_trajectories.append(traj)
            model_names.append(cif_path)

            print(f"Aligned: {cif_path}")

        except Exception as e:
            print(f"❌ Skipping {cif_path} due to error: {e}")
            verbose_print(f"Error with {cif_path}: {e}", args.verbose)

    if aligned_trajectories:
        save_with_model_names(aligned_trajectories, model_names, os.path.join(args.output, args.name))
        print(f"\n✅ Saved aligned multi-model PDB: {os.path.join(args.output, args.name)}")
    else:
        print("⚠️ No models were successfully aligned.")

    return os.path.join(args.output, args.name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Align substructures for all Boltz2 jobs to parent "
    "and save multi-model PDB. This assumes all run Boltz2 jobs are contained in"
    " a single directory with predictions located in subdirectories. Example usage " \
    "for PDB 3JQZ with all models in the 'boltz_jobs' directory:" \
    "python compile_best_model_structures -p 3JQZ.pdb -d boltz_jobs -o 3JQZ_top_models.pdb")
    parser.add_argument("-p", "--parent",
                        help="Path to the parent file to align to. Can be CIF, " \
                        "PDB, or MOL2. If no structure is given, will align to " \
                        "first._0.cif in predictions directory.")
    parser.add_argument("-d", "--directory", required=True,
                        help="Root directory to search for 'predictions' " \
                        "subdirectories containing CIF files.")
    parser.add_argument("-o", "--output", default="aligned_models.pdb",
                        help="Name for the output multi-model PDB file. " \
                        "Defaults to 'aligned_models.pdb'")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Enable verbose output.")
    parser.add_argument("--name", "-n", default="aligned_models.pdb",
                        help="Name for the output multi-model PDB file. " \
                        "Defaults to 'aligned_models.pdb'")
    args = parser.parse_args()

    main(args)