# use boltz_analysis_env 

import pandas as pd

def fetch_new(pred_df, metadata_df, protein):
    '''
    Checks if predictions of compounds are in predictions.csv records.
    Does NOT check if .cif file of final prediction is present.
    pred_df: df with names of molecules, confidence values and predictions
    metadata_df: df with all compound information
    protein: any name for the protein. taken from pdb file name. 
    Returns list of new compounds that have no predictions yet. 
    '''
    no_pred_prot = []
    
    for _, row in metadata_df.iterrows():
        substance_id = row['substance_id']
        inchi_key = row['inchi_key']
        smiles = row['smiles']
        
        pred_row = pred_df[(pred_df['substance_id'] == substance_id) & 
                                (pred_df['inchi_key'] == inchi_key)]
        
        if pred_row.empty or pd.isna(pred_row[f'pred_log_ic50_{protein}'].iloc[0]):
            no_pred_prot.append((substance_id, inchi_key, smiles))
        
    return no_pred_prot

def check_attempted(error_df, no_pred_prot):
    '''
    Checks errored.csv records to skip attempted ligands that failed in past runs.
    Returns updated list of tuples for tgcpl and hscpl seperately. 
    '''
    for _, row in error_df.iterrows():
        substance_id = row['substance_id']
        inchi_key = row['inchi_key']

        no_pred_prot = [lig for lig in no_pred_prot if not (lig[0] == substance_id and lig[1] == inchi_key)]
    
    return no_pred_prot