import os 
import shutil
import subprocess
from datetime import datetime
from BoltzCov.preprocessing import make_csv_for_yaml, setup_cov_yamls

def run_boltz_cov(prot_file, ligand_df, boltz_cache, res_idx, ligand_chain, VERBOSE, output_dir, slurm_template, msa_path=None):
    '''
    1. Generates a temp CSV for yaml building 
    2. Builds yamls
    3. write temp file for job list with yaml path and output dir
    4. submit jobs to cluster 
    ligand_df (DataFrame): A Pandas.DataFrame that has, at least, substance_id, smiles, and inchi_key columns. 
    boltz_cache (Path): Path to directory for storing temporary files. Will be deleted at end of pipeline run. 

    output_dir (Path): Path to output directory where final cleaned up predictions will be deposited. 
    '''
    protein_name = os.path.splitext(os.path.basename(prot_file))[0]
    missing_cols = [c for c in ['substance_id', 'smiles', 'inchi_key'] if c not in ligand_df.columns]
    if missing_cols:
        raise ValueError(f"LIGAND_CSV is missing expected columns: {missing_cols}")

    boltz_cache_prot = os.path.join(boltz_cache, protein_name)
    if VERBOSE: print("1. Generating temporary CSV with ligand and protein info for yaml building.")
    tmp_docking_csv = make_csv_for_yaml.generate_csv(prot_file, ligand_df, boltz_cache_prot, res_idx, ligand_chain, VERBOSE)
    
    yaml_cache = os.path.join(boltz_cache_prot, 'yaml_cache')
    # make yamls and get list of them
    yaml_list = setup_cov_yamls.create_boltz_yamls(csv_file=tmp_docking_csv, msa_path=msa_path, boltz_cache=boltz_cache_prot, output_dir=yaml_cache)

    job_list_file = os.path.join(boltz_cache_prot, "job_input_list.txt")
    if os.path.exists(job_list_file):
        os.remove(job_list_file) 

    for yaml in yaml_list:
            pkl_id = os.path.basename(yaml).replace(".yaml", "")
            
            now = datetime.now()
            datetime_str = now.strftime("%Y-%m-%d_%H-%M-%S") 

            pkl_datetime = "_".join([datetime_str, pkl_id])  
            pred_lig_dir = os.path.join(boltz_cache_prot, pkl_datetime)
            
            os.makedirs(pred_lig_dir, exist_ok=True)
            
            # moves yaml regardless of prediction/yaml existing
            shutil.move(yaml, pred_lig_dir)

            yaml_path = os.path.join(pred_lig_dir, os.path.basename(yaml))

            with open(job_list_file, 'a') as f: 
                f.write(f"{yaml_path} {pred_lig_dir}\n")
    
    
    slurm_script = os.path.join(boltz_cache_prot, os.path.basename(slurm_template))
    shutil.copy(slurm_template, slurm_script)
    # submit jobs
    subprocess.run(
         ["sbatch", slurm_script, job_list_file],
                check=True
        )
    
    if VERBOSE: print("Submitting jobs to cluster...")  

def check_pred():
    pass

def reorg_preds():
    # uses check_pred 
    pass

def update_preds():
    # update predictions in predictions.csv and update errored.csv for ligands that didnt work
    # should use both above functions 
    pass 

