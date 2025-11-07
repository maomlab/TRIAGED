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
import compile_best_model_structures as align

def get_interactions(interactions):
    """ Takes a PLIP interaction object and gets the number of each interaction. """
    num_saltbridges = len(interactions.saltbridge_lneg + interactions.saltbridge_pneg)
    num_hbonds = len(interactions.hbonds_ldon + interactions.hbonds_pdon) # ldon: ligand H donor, pdon: protein H donor
    num_pication = len(interactions.pication_laro + interactions.pication_paro)
    num_pistack = len(interactions.pistacking)
    num_halogen = len(interactions.halogen_bonds)
    num_waterbridges = len(interactions.water_bridges)

    interactions_vals = [num_saltbridges, num_hbonds, num_pication, num_pistack, num_halogen, num_waterbridges]
    return(interactions_vals)

def get_intearcting_residues(interactions):
    """ Takes a PLIP interaction object and gets the resids for residues in all interactions. """

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

    return(interacting_residues)

def verbose_print(msg, verbose):
    """Prints message if verbose is True."""
    if verbose:
        print(msg)

def find_first_model_files(root_dir):
    '''
    Finds the first model CIF files in 'predictions' subdirectories.
    :param root_dir: The root directory to start the search from.
        root_dir must follow this structure: <root_dir>/.../predictions/../model_output_0.cif

    :return: A list of paths to the first model CIF files.
    '''
    prediction_files = []
    for dirpath, dirnames, filenames in os.walk(root_dir):
        if 'predictions/' in dirpath and fnmatch.filter(filenames, '*_0.cif'):
            prediction_files.append(os.path.join(dirpath, fnmatch.filter(filenames, '*_0.cif')[0]))

    return prediction_files
    
def extract_name_key(filename):
    """Extracts the protein id and compound id from a filename of the form PROT_LIG_model_0.cif."""
    if "/" in filename:
        mod_filename = filename.split("/")[-1]
    else:
        mod_filename = filename
    match = re.search(r'_(\d+)\.cif$', mod_filename)
    if match:
        lig_name = mod_filename[0:match.span(0)[0]-6]
        return (lig_name, lig_name.split('_')[1][0:3])
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
    
def main(args):
     # Headers for the CSV output. Includes Ligand features as well as 
    # ligand-receptor interactionns
    headers = ["name","smiles", "inchi", "molwt", 
               "numheavy", "numrotbonds", "numrings", 
               "hydrophobicatoms", "hbondacceptors", 
               "saltbridges", "hbonds", "pication", 
               "pistack", "halogen", "waterbridge"]
    
    res_headers = ["name", "residue", "interaction_type"]

    temp_dir = tempfile.mkdtemp()

    # Changes config based on receptor type
    if args.receptor_type == "protein":
        config.DNARECEPTOR = False
    elif (args.receptor_type) == "rna" or (args.receptor_type == "dna"):
        config.DNARECEPTOR = True 
    else:
        print("You have given an incorrect receptor type! DNA, RNA, or PROTEIN only!")
        exit()

    #Gets all the cif files in a directory.
    cif_files = find_first_model_files(args.directory)
    collected_data = []
    residue_data = []
    errors = []
    for target_file in cif_files:
        # Gets the name and the model number from the .cif file.
        name, lig = extract_name_key(target_file)

        filetype = target_file.split(".")[-1].strip()
        target_pdb = os.path.join(temp_dir, f"{name}.pdb")
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
        plip_fprint = [name,
                    interactions.ligand.smiles.strip(),
                    interactions.ligand.inchikey.strip(),
                    interactions.ligand.molweight,
                    interactions.ligand.heavy_atoms,
                    interactions.ligand.num_rot_bonds,
                    interactions.ligand.num_rings,
                    len(interactions.ligand.hydroph_atoms),
                    len(interactions.ligand.hbond_acc_atoms)]
        for val in interactions_vals:
            plip_fprint.append(val) # just appending the number of each interaction type to end of fprint list
        collected_data.append(plip_fprint)

        residue_interactions = get_intearcting_residues(interactions)
        for resid, interaction_type in residue_interactions.items():
            res_fprint = [name, resid, interaction_type]
            residue_data.append(res_fprint)
        
        # pse for visualization
        if args.pymol_vis and args.parent:
            verbose_print("Generating PLIP visualizations...", args.verbose)
            vis_outdir = os.path.join(args.outdir, 'pse')
            os.makedirs(vis_outdir, exist_ok=True)
            align_args = argparse.Namespace(
                parent = args.parent,
                directory=os.path.dirname(target_file),
                output=vis_outdir,
                verbose=args.verbose,
                name=f"{name}_aligned.pdb"
            )
            try:
                aligned_target = align.main(align_args)
            except Exception as e:
                print(f"Alignment failed for {target_file}: {e}")
                sys.exit(1) 
            try:
                vis_cmd = ["plip", "-f", aligned_target, "-y", "-o", vis_outdir]
                subprocess.run(vis_cmd, check=True)
            except subprocess.CalledProcessError as e:
                print(f"PLIP failed for {target_file}: {e}")
                errors.append(target_file)
                continue   
        
    csv_out = os.path.join(args.outdir, f"{args.csv_name}.csv")
    verbose_print("Writing out number of interactions...", args.verbose)
    with open(csv_out, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(headers)
        writer.writerows(collected_data)
    
    residues_csv_out = os.path.join(args.outdir, f"{args.csv_name}_residues.csv")
    verbose_print("Writing out residues of interactions...", args.verbose)
    with open(residues_csv_out, 'w', newline='') as res_csv:
        writer = csv.writer(res_csv)
        writer.writerow(res_headers)
        writer.writerows(residue_data)

    # Clean up temporary directory
    shutil.rmtree(temp_dir)

    return errors 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyzes a collection of .cif files"
                            "from a boltz2 workflow using PLIP and outputs" \
                            "ligand information as well as residue interaction information" \
                            "calculated by PLIP in a CSV.")
    parser.add_argument("-d", "--directory", required=True,
                        help="Root directory of .cif files for analysis.")
    parser.add_argument("-o", "--outdir", required=True,
                        help="Path to output directory." \
                        "Defaults to 'compiled_plip_fprints.csv'")
    parser.add_argument("-rt", "--receptor_type", required=True,
                        help="Type of receptor - 'dna', 'rna', or 'protein.'", 
                        choices=['dna', 'rna', 'protein'])
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Enable verbose output.")
    parser.add_argument("-y", "--pymol_vis", default=False, action="store_true",
                        help="If set, will generate pymol visualization of the interactions.")
    parser.add_argument("-p", "--parent", required=False, default="",
                        help="Path to the parent structure file to align to for pymol visualizations." \
                        "Can be CIF or PDB. If not given, will not generate pymol visualizations.")
    parser.add_argument("-n", "--csv_name", required=False, default="plip_fingerprints",
                        help="Name of the csv files (don't include csv extension).")

    args = parser.parse_args()
    main(args)