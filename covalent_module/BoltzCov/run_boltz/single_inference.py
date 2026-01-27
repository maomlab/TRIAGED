# need it to be such that smiles and chain info can be provided in a single command and docking should happen 
import argparse
import torch
import os
import subprocess
import yaml
from covalent_module.preprocessing import setup_cov_yamls
from covalent_module.preprocessing import covalent_utils
from covalent_module.preprocessing import make_input_csv

from dotenv import load_dotenv
load_dotenv("TRIAGED/covalent_module/.env")

# load common env variables and checks if the values exist
expected_vars = ["BOLTZ_CACHE", "CCD_DB"]
missing_vars = [var for var in expected_vars if var not in os.environ]
if missing_vars: 
    raise EnvironmentError(f"Missing expected environment variables: {missing_vars}. Please fill out .env in working directory.")

def write_yaml(outdir, smiles, compound_id, pdb_file, chain_id, res_idx, msa=None):
    CCD_DB = os.environ.get("CCD_DB", outdir)
    smiles = setup_cov_yamls.check_smiles(smiles) # returns conancial smiles or None
    
    if smiles is None:
        print(f"[ERROR] Invalid SMILES.")
        return False
        
    # remove leaving group and make pkl file
    smiles_no_lg, lig_atom, _ = covalent_utils.remove_leaving_group(smiles)
    covalent_utils.process_covalent_smiles(CCD_DB, smiles_no_lg, compound_id=compound_id)

    # get protein info 
    sequence, res_name, res_atom = make_input_csv.process_protein(pdb=pdb_file, idx=res_idx, lig_chain=chain_id)
    pdb_name = os.path.basename(pdb_file).replace(".pdb","")

    protein_data = {
        "id": "A",
        "sequence": sequence,
        "modification": [{"position": int(res_idx), "ccd": res_name}],
    }

    if msa is not None:
        protein_data["msa"] = msa

    data = {
        "sequences": [
            {"protein": protein_data},
            {"ligand": {"id": "LIG", "ccd": str(compound_id)}},
        ],
        "constraints": [
            {
                "bond": {
                    "atom1": ["A", int(res_idx), res_atom],
                    "atom2": ["LIG", int(1), lig_atom],
                }
            }
        ],
        "properties": [
            {"affinity": {"binder": "LIG"}}
        ],
    }
    yaml_file = os.path.join(outdir, f"{pdb_name}_{compound_id}.yaml") # should be unique for each ligand
    with open(yaml_file, "w") as f:
        yaml.safe_dump(
            data, 
            f,
            sort_keys=False,
            indent=4,
            width=4096,  
            default_flow_style=False
        )
    return yaml_file

def run(pdb_file,res_idx,chain_id,smiles,compound_id,outdir, msa=None):
    '''
    Script runs code needed to make yaml for a single ligand and runs inference with boltz.
    '''

    if len(compound_id) > 5: 
        raise ValueError("Compound ID must be 5 characters or less!!")
    
    os.makedirs(outdir, exist_ok=True)

    try:
        yaml_path = write_yaml(outdir,smiles,compound_id,pdb_file,chain_id,res_idx,msa)
    except Exception as e:
        print(f"[ERROR] Error in YAML creation: {e}")
        return
    
    boltz_job = os.path.join(outdir, "boltz_job.sh")
    # write a bash script and run it 
    with open(boltz_job, "w") as boltz:
        boltz.write("#!/bin/bash\n")
        # if slurm not None, slurm args are added later over here 
        boltz.write("source  ${HOME}/opt/miniconda3/etc/profile.d/conda.sh\n")       
        boltz.write("conda activate boltz2\n")
        boltz.write("module load cuda cudnn\n")
        boltz.write("\n")
        
        out_stdout = os.path.join(outdir, "boltz.out")
        err_stdout = os.path.join(outdir, "boltz.err")

        # Run boltz predict
        boltz.write(f"boltz predict {yaml_path} --out_dir {outdir} --num_workers 8 --use_msa_server 1> {out_stdout} 2> {err_stdout}\n")
        boltz.write("\n")

        # File existence check (quote {outdir} to handle spaces safely)
        boltz.write(f"FILES=$(compgen -G '{outdir}/boltz_results_*/predictions/*/*_model_0.cif')\n")
        boltz.write("\n")
        boltz.write('if [ -z "$FILES" ]; then\n')
        boltz.write('    echo "Job failed or output files not found"\n')
        boltz.write('else\n')
        boltz.write('    echo "Boltz prediction completed successfully."\n')
        boltz.write('fi\n')

    SLURM = os.environ.get("SLURM_TEMPLATE", None)
    if SLURM is not None: 
        with open(SLURM, "r") as f:
            slurm_lines = f.read()

        with open(boltz_job, "r") as f:
            content = f.readlines()

        # Insert slurm lines after the first line in boltz_job
        content.insert(1, slurm_lines + "\n")
        with open(boltz_job, "w") as f:
            f.writelines(content)
        subprocess.run(["sbatch", boltz_job], check=True)

    else: # if not running with slurm 
        if torch.cuda.is_available():
            print(torch.cuda.is_available())
            print(torch.cuda.device_count())
            # in interactive shell
            subprocess.run(["bash", boltz_job], check=True)
        else: 
            print(torch.cuda.is_available())
            print(torch.cuda.device_count())
            print("No GPU found! Please use SLURM to submit to cluster or login to a GPU node.")

def main(args):
    run(pdb_file=args.prot_file, res_idx=args.res_idx, chain_id=args.lig_chain, smiles=args.smiles, compound_id=args.compound_id, msa=args.msa_path,outdir=args.outdir) 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Will create yaml and script to run Boltz2 for covalent docking. Please use ccd_pkl conda env.")
    parser.add_argument("--prot_file", type=str, required=True, help="Path to a PDB file.")
    parser.add_argument("--res_idx", type=int, required=True, help="Residue index to which the ligand needs to be docked. Starting at 1.")
    parser.add_argument("--lig_chain", type=str, required=False, help="Chain interacting with ligand in PDB file. Single character.", default='A')
    parser.add_argument("--smiles", type=str, required=True, help=" SMILES of ligand.")
    parser.add_argument("--compound_id", type=str, required=True, help="Compound ID for the ligand. Must be 5 characters or less.")
    parser.add_argument("--msa_path", type=str, required=False, help="Path to MSA file in csv format. If provided, will be added to yaml.", default=None)
    parser.add_argument("-o","--outdir", type=str, required=True, help="Output directory for Boltz.")
    
    args=parser.parse_args()
    
    main(args)
    
