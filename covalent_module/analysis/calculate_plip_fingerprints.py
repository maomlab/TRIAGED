from plip.basic import config
from plip.structure.preparation import PDBComplex
import argparse
import fnmatch
import re
import gemmi
import os
import glob
import tempfile
import csv

def parse_args():
    parser = argparse.ArgumentParser(description="Analyzes a collection of .cif files"
                            "from a boltz2 workflow using PLIP and outputs" \
                            "ligand information as well as interaction information" \
                            "calculated by PLIP in a CSV.")
    parser.add_argument("-d", "--directory", required=True,
                        help="Root directory of .cif files for analysis.")
    parser.add_argument("-o", "--output", default="compiled_plip_fprints.csv",
                        help="Name for the output output CSV file. " \
                        "Defaults to 'compiled_plip_fprints.csv'")
    parser.add_argument("-rt", "--receptor_type", required=True,
                        help="Type of receptor - 'dna', 'rna', or 'protein.'", 
                        choices=['dna', 'rna', 'protein'])
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Enable verbose output.")

    return parser.parse_args()


def get_interactions(interactions):
    """ Takes a PLIP interaction object and gets the number of each interaction. """
    num_saltbridges = len(interactions.saltbridge_lneg + interactions.saltbridge_pneg)
    num_hbonds = len(interactions.hbonds_ldon + interactions.hbonds_pdon)
    num_pication = len(interactions.pication_laro + interactions.pication_paro)
    num_pistack = len(interactions.pistacking)
    num_halogen = len(interactions.halogen_bonds)
    num_waterbridges = len(interactions.water_bridges)

    interactions_vals = [num_saltbridges, num_hbonds, num_pication, num_pistack, num_halogen, num_waterbridges]
    return(interactions_vals)


def verbose_print(msg, verbose):
    """Prints message if verbose is True."""
    if verbose:
        print(msg)

def find_first_model_files(root_dir):
    '''
    Finds the first model CIF files in 'predictions' subdirectories.
    :param root_dir: The root directory to start the search from.
        root_dir must follow this structure: <root_dir>/<target_id>/boltz_results_<target_id>/predictions/<target_id>/<model_output_files>

    :return: A list of paths to the first model CIF files.
    '''
    prediction_files = []
    for dirpath, dirnames, filenames in os.walk(root_dir):
        if 'predictions/' in dirpath and fnmatch.filter(filenames, '*_0.cif'):
            prediction_files.append(os.path.join(dirpath, fnmatch.filter(filenames, '*_0.cif')[0]))

    return prediction_files
    
def extract_name_key(filename):
    """Extracts the name and numeric key from a filename of the form NAME_model_KEY.cif."""
    if "/" in filename:
        mod_filename = filename.split("/")[-1]
    else:
        mod_filename = filename
    match = re.search(r'_(\d+)\.cif$', mod_filename)
    if match:
        lig_name = mod_filename[0:match.span(0)[0]-6]
        return (int(match.group(1)), mod_filename[0:match.span(0)[0]-6], lig_name.split('_')[1][0:3])
    else:
        raise ValueError(f"Filename {filename} does not match the expected pattern.")

def convert_cif_to_pdb(cif_path, pdb_path, chain_map=None):
    """Convert mmCIF to PDB using gemmi, fixing long chain names."""
    doc = gemmi.cif.read_file(cif_path)
    block = doc.sole_block()
    structure = gemmi.make_structure_from_block(block)

    # Fix chain names longer than 1 character
    for model in structure:
        for chain in model:
            if len(chain.name) != 1:
                # map to a new chain name if provided, or default to 'X'
                new_chain = chain_map.get(chain.name, 'X') if chain_map else 'X'
                chain.name = new_chain

    structure.write_pdb(pdb_path)

def convert_pdb_to_pdb(source_pdb_path, target_pdb_path):
    """Standardize or reformat a PDB file."""
    # Read the structure from the source PDB
    structure = gemmi.read_structure(source_pdb_path)
    
    # Write the structure to the target PDB
    structure.write_pdb(target_pdb_path)


if __name__ == "__main__":
    # Headers for the CSV output. Includes Ligand features as well as 
    # ligand-receptor interactionns
    headers = ["name", "modelnum","smiles", "inchi", "molwt", 
               "numheavy", "numrotbonds", "numrings", 
               "hydrophobicatoms", "hbondacceptors", 
               "saltbridges", "hbonds", "pication", 
               "pistack", "halogen", "waterbridge"]
    
    args = parse_args()
    temp_dir = tempfile.mkdtemp()

    # Changes config based on receptor type
    if args.receptor_type == "protein":
        config.DNARECEPTOR = False
    elif (args.receptor_type) == "rna" or (args.receptor_type == "dna"):
        config.DNARECEPTOR = True 
    else:
        print("You have given an incorrect receptor type! DNA, RNA, or PROTEIN only!")
        exit()
    import ipdb; ipdb.set_trace()
    #Gets all the cif files in a directory.
    # cif_files = glob.glob(os.path.join(args.directory, "*.cif"))
    cif_files = find_first_model_files(args.directory)
    collected_data = []
    for target_file in cif_files:
        # Gets the name and the model number from the .cif file.
        id, name, lig = extract_name_key(target_file)

        # PLIP only works on PDBs, so converts the .cif to a temporary .pdb
        # PDB inputs, when supported, also need to be treated to be readable.
        filetype = target_file.split(".")[-1].strip()
        target_pdb = os.path.join(temp_dir, "target.pdb")
        if filetype == "cif":
            convert_cif_to_pdb(target_file, target_pdb)
        elif filetype == "pdb":
            convert_pdb_to_pdb(target_file, target_pdb)
        else:
            print("Unsupported parent file type! Must be PDB or CIF!")
            exit()

        # Loads the temporary pdb into PLIP and does the interaction analysis
        my_mol = PDBComplex()
        my_mol.load_pdb(target_pdb)
        my_mol.analyze()

        interactions = my_mol.interaction_sets[f"{lig}:X:1"]
        
        interactions_vals = get_interactions(interactions)

        #Organizes everything into one line for the CSV
        plip_fprint = [name, id,
                    interactions.ligand.smiles.strip(),
                    interactions.ligand.inchikey.strip(),
                    interactions.ligand.molweight,
                    interactions.ligand.heavy_atoms,
                    interactions.ligand.num_rot_bonds,
                    interactions.ligand.num_rings,
                    len(interactions.ligand.hydroph_atoms),
                    len(interactions.ligand.hbond_acc_atoms)]
        
        for val in interactions_vals:
            plip_fprint.append(val)

        collected_data.append(plip_fprint)

    # Writes out the info to the final CSV.
    with open(args.output, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(headers)
        writer.writerows(collected_data)