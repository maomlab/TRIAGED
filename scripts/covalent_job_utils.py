import os 
import shutil
import subprocess 
from pathlib import Path
import glob
import argparse
from setup_boltz_job import create_boltz_job
        
def setup_prot_dir(prot_list, db_dir, working_dir):
    """
    Given a text file with a list of PDB codes and path to PDB files, 
    this function creates a directory for the protein along with its 
    PDB file and yaml file. 
    File structure output: 
    1XYZ/
    ├── 1XYZ.pdb
    ├── 1XYZ.yaml 
    Args:
        prot_list (str): Path to txt file with list of PDB codes for inference with boltz. 
        db_dir (str): Path to CovalentInDB2.0 where all covalent protein PDB files are available. 
        working_dir (str): Path to where PDB directories should be made for boltz output. 
    """
    with open(prot_list, 'r') as prot_file: 
        prots = prot_file.readlines()

        for prot in prots: 
            prot = prot.strip() # remove trailing \n
            print(prot)
            pdb_dir = os.path.join(working_dir, prot)
            os.makedirs(pdb_dir, exist_ok=True)
            shutil.copy(os.path.join(db_dir, f"{prot}.pdb"), pdb_dir)
            os.chdir(pdb_dir)
            
            pdb_file = os.path.join(pdb_dir, f"{prot}.pdb")

            if not any(Path(pdb_dir).glob("*.y*ml")): # create the yaml file only if it dne
                create_boltz_job(pdb_file=pdb_file, covalent_docking=True, output_dir=pdb_dir)

def generate_job_list(base_dir, output_file):
    """
    Generate a two-column job input list for SLURM array jobs.
    Each line of the output file will have:
        <full_path_to_yaml> <receptor_name>
    Args:
        base_dir (str): Path to the parent directory containing receptor subdirectories.
                        Each subdirectory should have at least one .yaml file.
        output_file (str): Path to the output text file (job_input_list.txt).
    Notes:
        - Only the first .yaml file found in each receptor directory is included.
        - Receptor name is taken from the subdirectory name.
        - Output is sorted alphabetically by receptor directory name.
    """
    with open(output_file, "w") as out_f:
        for receptor_dir in sorted(os.listdir(base_dir)):
            dir_path = os.path.join(base_dir, receptor_dir)
            if os.path.isdir(dir_path):
                yaml_files = [f for f in os.listdir(dir_path) if f.endswith(".yaml")]
                if yaml_files:
                    yaml_file = yaml_files[0]  # take the first YAML file
                    yaml_path = os.path.join(dir_path, yaml_file)
                    out_f.write(f"{yaml_path} {receptor_dir}\n")

def submit_slurm_job(slurm_script, job_file):
    """
    Submit a SLURM job script with a custom job input list file.
    Args:
        slurm_script (str): Path to the SLURM script (.sh).
        job_file (str): Path to the job input list file.
    """

    result = subprocess.run(
        ["sbatch", slurm_script, job_file],
        check=True
    )

    print("Submitting SLURM jobs...")

    return result