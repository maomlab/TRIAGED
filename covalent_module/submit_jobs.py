import os 
import shutil
import subprocess
import argparse
from covalent_module.preprocessing import make_input_csv 
from covalent_module.preprocessing import setup_cov_yamls

# ccd_pkl conda env 
# last tested: 10/06/25 

def run(name, prot_file, res_idx, lig_chain, lig_csv, outdir, ccd_db, msa_path=None, slurm_template="covalent_module/covalent_slurm.sh", debug=False):
    if not debug:
        os.makedirs(outdir, exist_ok=True)

        # generate csv for yaml building
        base_name = os.path.basename(lig_csv).replace(".csv", "")
        
        gen_csv = os.path.join(outdir, f"generated_{base_name}_{name}.csv")
        make_input_csv.generate_csv(name=name,prot_file=prot_file, res_idx=int(res_idx), lig_chain=str(lig_chain), lig_csv=lig_csv, out_csv=gen_csv, ccd_db=ccd_db)
    
        # build yamls
        yaml_list = setup_cov_yamls.create_boltz_yamls(csv_file=gen_csv, output_dir=outdir, msa=msa_path)

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
                print(f"[WARNING] {os.path.basename(yaml)} exists. Deleting and rewriting.")
                os.remove(dest)  
            shutil.move(yaml, pred_lig_dir)

            yaml_path = os.path.join(pred_lig_dir, os.path.basename(yaml))

            with open(job_list_file, 'a') as f: 
                f.write(f"{yaml_path} {pred_lig_dir}\n")

        # submit jobs
        subprocess.run(
                ["sbatch", slurm_template, job_list_file],
                check=True
            )

        print("Submitting SLURM jobs...")
    else:
        print("Debug mode is ON. No files were written.")
        base_name = os.path.basename(lig_csv).replace(".csv", "")
        print(f"Using generated CSV at: {os.path.join(outdir, f'generated_{base_name}_{name}.csv')}")
        print(f"Using create yamls in: {outdir}")
        print(f"Using pkls in: {ccd_db}")

        job_list_file = os.path.join(outdir, "job_input_list.txt")

        # submit jobs
        subprocess.run(
                ["sbatch", slurm_template, job_list_file],
                check=True
            )
        
        print(f"Submitting SLURM jobs...")
    
def main():
    parser = argparse.ArgumentParser(description="One main run script for covalent docking. Generates input CSV, builds yamls, and submits SLURM jobs. Please use ccd_pkl conda env.")
    parser.add_argument("-n", "--name", type=str, required=True, help="Name of protein. Used for naming output files.")
    parser.add_argument("-p","--prot_file", type=str, required=True, help="Path to either a PDB file or a TXT file with a single chain sequence.")
    parser.add_argument("-r", "--res_idx", type=int, required=True, help="Index of the residue to be covalently targeted by a covalent ligand. Starting at 1.")
    parser.add_argument("-c", "--lig_chain", type=str, required=True, help="Chain interacting with ligand in PDB file. Single character.")
    parser.add_argument("-l","--lig_csv", type=str, required=True, help="Path to CSV with Ligand info.")
    parser.add_argument("-o","--outdir", type=str, required=True, help="Output directory for all jobs.")
    parser.add_argument("-s","--slurm", type=str, required=False, help="Path to SLURM template file.", default="covalent_module/covalent_slurm.sh")
    parser.add_argument("-m", "--msa_path", type=str, required=False, help="Path to MSA file in a3m format. If provided, will be added to yaml.", default=None)
    parser.add_argument("-d", "--debug", action="store_true", help="If set, files will not be rewritten, just submitting slurm boltz job.")

    args=parser.parse_args()

    setup_cov_yamls.ensure_environment_variables()
    ccd_db = os.getenv("CCD_DB") # set in setup_enviorment.sh

    run(name=args.name, prot_file=args.prot_file, res_idx=args.res_idx, lig_chain=args.lig_chain, lig_csv=args.lig_csv, outdir=args.outdir, ccd_db=ccd_db, slurm_template=args.slurm, msa_path=args.msa_path, debug=args.debug)

if __name__ == "__main__":
    main()


