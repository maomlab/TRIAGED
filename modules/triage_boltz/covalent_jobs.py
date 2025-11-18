import os 
import glob
from setup_boltz_job import create_boltz_job

def setup_covalent_yamls(query_dir):
    """
    Goes through a list of directories with PDB files and sets up yaml files and SLURM job files. 
    query_dir (str): path to directory with directories. 
    File structure expected: 
    query_dir/
    ├── 1XYZ/
    │   └── 1XYZ.pdb
    """
    for pdb_dir in glob.glob(os.path.join(query_dir,"*/")):
        pdb_files = glob.glob(os.path.join(pdb_dir, "*.pdb"))

        if not pdb_files:
            print(f"{pdb_files} not found")
            continue

        pdb_file = pdb_files[0]

        # will rewrite existing yamls! 
        create_boltz_job(pdb_file=pdb_file, covalent_docking=True, output_dir=pdb_dir)
        

def setup_slurm_jobs():
    pass 


# testing
setup_covalent_yamls("/home/ymanasa/turbo/ymanasa/opt/boltz/covalent_testing")