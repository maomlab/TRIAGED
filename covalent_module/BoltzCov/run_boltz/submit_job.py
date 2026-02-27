import os 
import shutil
import subprocess
from datetime import datetime
import time
from BoltzCov.preprocessing import make_csv_for_yaml, setup_cov_yamls, setup_noncov_yamls

def wait_for_slurm_job(job_id, check_interval=30):
    '''
    Polls SLURM queue until job completes.
    
    :param job_id: SLURM job ID to monitor
    :param check_interval: Seconds between status checks (default: 30)
    :return: Final job status (COMPLETED, FAILED, etc.)
    '''
    import subprocess
    import time
    
    print(f"Monitoring job {job_id}...")
    
    while True:
        check = subprocess.run(
            ["squeue", "-j", job_id, "-h", "--array"],
            capture_output=True,
            text=True
        )
        
        if not check.stdout.strip():
            print(f"Job {job_id} finished")
            break
        
        print(f"Job {job_id} still running...")
        time.sleep(check_interval)
    
    status = subprocess.run(
        ["sacct", "-j", job_id, "--format=State", "--noheader"],
        capture_output=True,
        text=True
    )
    
    final_status = status.stdout.strip().split('\n')[0].strip()
    print(f"Final status: {final_status}")
    
    return final_status

def run_boltz_cov(prot_file, ligand_df, boltz_cache, run_cache, res_idx, ligand_chain, VERBOSE, slurm_template, msa_path=None):
    '''
    1. Generates a temp CSV for yaml building 
    2. Builds yamls
    3. write temp file for job list with yaml path and output dir
    4. submit jobs to cluster 
    ligand_df (DataFrame): A Pandas.DataFrame that has, at least, substance_id, smiles, and inchi_key columns. 
    run_cache (Path): Path to directory for storing temporary files. Will be deleted at end of pipeline run. 
    boltz_cache(Path): Path to directorty for Boltz to download weights for models and PDB pkl files. 
    '''
    protein_name = os.path.splitext(os.path.basename(prot_file))[0]
    missing_cols = [c for c in ['substance_id', 'smiles', 'inchi_key'] if c not in ligand_df.columns]
    if missing_cols:
        raise ValueError(f"LIGAND_CSV is missing expected columns: {missing_cols}")

    run_cache_prot = os.path.join(run_cache, protein_name)
    if VERBOSE: print("1. Generating temporary CSV with ligand and protein info for yaml building.")
    tmp_docking_csv = make_csv_for_yaml.generate_csv(prot_file=prot_file, ligand_df=ligand_df, run_cache=run_cache_prot, boltz_cache=boltz_cache, res_idx=res_idx, ligand_chain=ligand_chain, VERBOSE=VERBOSE)

    yaml_cache = os.path.join(run_cache_prot, 'yaml_cache')
    # make yamls and get list of them
    yaml_list = setup_cov_yamls.create_boltz_yamls(csv_file=tmp_docking_csv, msa_path=msa_path, boltz_cache=boltz_cache, output_dir=yaml_cache)

    job_list_file = os.path.join(run_cache_prot, "job_input_list.txt")
    if os.path.exists(job_list_file):
        os.remove(job_list_file) 

    for yaml in yaml_list:
            pkl_id = os.path.basename(yaml).replace(".yaml", "")
            
            now = datetime.now()
            datetime_str = now.strftime("%Y%m%d_%H%M%S") 

            pkl_datetime = "_".join([datetime_str, pkl_id])  
            pred_lig_dir = os.path.join(run_cache_prot, pkl_datetime)
            
            os.makedirs(pred_lig_dir, exist_ok=True)
            
            # copy yaml regardless of prediction/yaml existing
            shutil.copy(yaml, pred_lig_dir)

            yaml_path = os.path.join(pred_lig_dir, os.path.basename(yaml))

            with open(job_list_file, 'a') as f: 
                f.write(f"{yaml_path} {pred_lig_dir}\n")
    
    if os.path.exists(yaml_cache): shutil.rmtree(yaml_cache, ignore_errors=True)
    slurm_script = os.path.join(run_cache_prot, os.path.basename(slurm_template))
    shutil.copy(slurm_template, slurm_script)
    # submit jobs
    print("Submitting jobs to cluster in 5 seconds...") 
    time.sleep(5)
    result = subprocess.run(["sbatch", slurm_script, job_list_file, boltz_cache], capture_output=True, text=True, check=True)
    job_id = result.stdout.strip().split()[-1]
    final_status = wait_for_slurm_job(job_id, check_interval=120)

    return final_status 

def run_boltz_noncov(prot_file, ligand_df, boltz_cache, run_cache, VERBOSE, slurm_template, msa_path=None):
    '''
    1. Builds non-cov yamls
    2. write temp file for job list with yaml path and output dir
    3. submit jobs to cluster 
    ligand_df (DataFrame): A Pandas.DataFrame that has, at least, substance_id, smiles, and inchi_key columns. 
    run_cache (Path): Path to directory for storing temporary files. Will be deleted at end of pipeline run. 
    boltz_cache(Path): Path to directorty for Boltz to download weights for models and PDB pkl files. 
    '''
    protein_name = os.path.splitext(os.path.basename(prot_file))[0]
    missing_cols = [c for c in ['substance_id', 'smiles', 'inchi_key'] if c not in ligand_df.columns]
    if missing_cols:
        raise ValueError(f"LIGAND_CSV is missing expected columns: {missing_cols}")

    run_cache_prot = os.path.join(run_cache, protein_name)
    yaml_cache = os.path.join(run_cache_prot, 'yaml_cache')
    # make yamls and get list of them
    if VERBOSE: print("1. Generating non-covalent yamls.")
    yaml_list = setup_noncov_yamls.create_boltz_yamls(prot_file=prot_file,ligand_df=ligand_df, msa_path=msa_path, output_dir=yaml_cache)

    job_list_file = os.path.join(run_cache_prot, "job_input_list.txt")
    if os.path.exists(job_list_file):
        os.remove(job_list_file) 

    for yaml in yaml_list:
            pkl_id = os.path.basename(yaml).replace(".yaml", "")
            
            now = datetime.now()
            datetime_str = now.strftime("%Y%m%d_%H%M%S") 

            pkl_datetime = "_".join([datetime_str, pkl_id])  
            pred_lig_dir = os.path.join(run_cache_prot, pkl_datetime)
            
            os.makedirs(pred_lig_dir, exist_ok=True)
            
            # copy yaml regardless of prediction/yaml existing
            shutil.copy(yaml, pred_lig_dir)

            yaml_path = os.path.join(pred_lig_dir, os.path.basename(yaml))

            with open(job_list_file, 'a') as f: 
                f.write(f"{yaml_path} {pred_lig_dir}\n")
    
    if os.path.exists(yaml_cache): shutil.rmtree(yaml_cache, ignore_errors=True)
    slurm_script = os.path.join(run_cache_prot, os.path.basename(slurm_template))
    shutil.copy(slurm_template, slurm_script)
    # submit jobs
    print("Submitting jobs to cluster in 5 seconds...") 
    time.sleep(5)
    result = subprocess.run(["sbatch", slurm_script, job_list_file, boltz_cache], capture_output=True, text=True, check=True)
    job_id = result.stdout.strip().split()[-1]
    final_status = wait_for_slurm_job(job_id, check_interval=120)

    return final_status 