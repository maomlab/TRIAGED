import os 
from BoltzCov.preprocessing import setup_cov_yamls, make_csv_for_yaml


def run_boltz_cov(ligand_df, boltz_cache, output_dir, VERBOSE):
    '''
    1. Generates a temp CSV for yaml building 
    2. Builds yamls
    3. Checks if .cif file exists in outdir 
    4. write temp file for job list with yaml path and output dir ? (can put both tgcpl and hscpl i think)
    5. submit jobs to cluster 
    ligand_df (DataFrame): A Pandas.DataFrame that has, at least, substance_id, smiles, and inchi_key columns. 
    boltz_cache (Path): Path to directory for storing temporary files. Will be deleted at end of pipeline run. 

    output_dir (Path): Path to output directory where final cleaned up predictions will be deposited. 
    '''
    missing_cols = [c for c in ['substance_id', 'smiles', 'inchi_key'] if c not in ligand_df.columns]
    if missing_cols:
        raise ValueError(f"LIGAND_CSV is missing expected columns: {missing_cols}")
    
    os.makedirs(boltz_cache, exist_ok=True)

    if VERBOSE: print("1. Generating temporary CSV for yaml building.")
    make_csv_for_yaml.generate_csv()

    yaml_list = setup_cov_yamls.create_boltz_yamls(g)






    
    pass


def check_pred():
    pass

def reorg_preds():
    # uses check_pred 
    pass

def update_preds():
    # update predictions in predictions.csv and update errored.csv for ligands that didnt work
    # should use both above functions 
    pass 

