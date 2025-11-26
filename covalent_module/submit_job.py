import os 
import shutil
import subprocess
import argparse
import pandas as pd
from covalent_module.preprocessing import make_input_csv 
from covalent_module.preprocessing import setup_cov_yamls
from dotenv import load_dotenv

load_dotenv("TRIAGED/covalent_module/.env")
# load common env variables and checks if the values exist
expected_vars = ["BOLTZ_CACHE", "LIGAND_CSV", "SLURM_TEMPLATE"]
missing_vars = [var for var in expected_vars if var not in os.environ]
if missing_vars: 
    raise EnvironmentError(f"Missing expected environment variables: {missing_vars}. Please fill out .env in working directory.")

# ccd_pkl conda env 
def run(name, prot_file, res_idx, lig_chain, outdir, msa_path):
    # get environmental variables from setup_enviorment.env
    CCD_DB = os.environ.get("CCD_DB", outdir) 
    VERBOSE = os.environ.get("VERBOSE", "FALSE").upper() == "TRUE"
    SLURM_TEMPLATE= os.environ.get("SLURM_TEMPLATE")
    LIGAND_CSV = os.environ.get("LIGAND_CSV")
    DEBUG = os.environ.get("DEBUG", "FALSE").upper() == "TRUE"

    # check cols in lig csv
    lig_csv_cols = pd.read_csv(LIGAND_CSV, nrows=0).columns 
    missing_cols = [c for c in ['vault_id', 'SMILES'] if c not in lig_csv_cols]
    if missing_cols:
        raise ValueError(f"LIGAND_CSV is missing expected columns: {missing_cols}")

    if not DEBUG:
        os.makedirs(outdir, exist_ok=True)
        shutil.copy(SLURM_TEMPLATE, outdir) 

        # generate csv for yaml building
        base_name = os.path.basename(LIGAND_CSV).replace(".csv", "")
        gen_csv = os.path.join(outdir, f"generated_{base_name}_{name}.csv")
        make_input_csv.generate_csv(name, prot_file, int(res_idx), str(lig_chain), gen_csv, CCD_DB)
    
        # build yamls
        yaml_list = setup_cov_yamls.create_boltz_yamls(gen_csv, outdir, msa_path)

        job_list_file = os.path.join(outdir, "job_input_list.txt")
        if os.path.exists(job_list_file):
            os.remove(job_list_file) # because append is used later

        for yaml in yaml_list:
            base_name = os.path.basename(yaml).replace(".yaml", "")
            receptor = base_name.split("_")[0]
            receptor_lig = f'{base_name.split("_")[0]}_{base_name.split("_")[1]}'

            pred_lig_dir = os.path.join(outdir, receptor, receptor_lig)

            os.makedirs(pred_lig_dir, exist_ok=True)
            dest = os.path.join(pred_lig_dir, os.path.basename(yaml))
            if os.path.exists(dest):
                if VERBOSE: print(f"[WARNING] {os.path.basename(yaml)} exists. Deleting and rewriting.")
                os.remove(dest)  
            shutil.move(yaml, pred_lig_dir)

            yaml_path = os.path.join(pred_lig_dir, os.path.basename(yaml))

            with open(job_list_file, 'a') as f: 
                f.write(f"{yaml_path} {pred_lig_dir}\n")

        slurm_script = os.path.join(outdir, os.path.basename(SLURM_TEMPLATE))
        # submit jobs
        subprocess.run(
                ["sbatch", slurm_script, job_list_file],
                check=True
            )
        if VERBOSE: print("Submitting SLURM jobs...")
    elif DEBUG:
        print("Debug mode is ON. No files were written.")
        base_name = os.path.basename(LIGAND_CSV).replace(".csv", "")
        print(f"Using generated CSV at: {os.path.join(outdir, f'generated_{base_name}_{name}.csv')}")
        print(f"Using create yamls in: {outdir}")
        print(f"Using pkls in: {CCD_DB}")

        job_list_file = os.path.join(outdir, "job_input_list.txt")

        # submit jobs
        subprocess.run(
                ["sbatch", SLURM_TEMPLATE, job_list_file],
                check=True
            )
        
        if VERBOSE: print(f"Submitting SLURM jobs...")
    
def main(args):
    run(name=args.name, prot_file=args.prot_file, res_idx=args.res_idx, lig_chain=args.lig_chain, outdir=args.outdir, msa_path=args.msa_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="One main run script for covalent docking. Generates input CSV, builds yamls, and submits SLURM jobs. Please use ccd_pkl conda env.")
    parser.add_argument("-n", "--name", type=str, required=True, help="Name of protein. Used for naming output files.")
    parser.add_argument("--prot_file", type=str, required=True, help="Path to either a PDB file or a TXT file with a single chain sequence.")
    parser.add_argument("--res_idx", type=int, required=True, help="Index of the residue to be covalently targeted by a covalent ligand. Starting at 1.")
    parser.add_argument("--lig_chain", type=str, required=True, help="Chain interacting with ligand in PDB file. Single character.")
    parser.add_argument("--msa_path", type=str, required=False, help="Path to MSA file in csv format. If provided, will be added to yaml.", default=None)
    parser.add_argument("-o","--outdir", type=str, required=True, help="Output directory for all jobs.")

    args=parser.parse_args()
    
    main(args)


