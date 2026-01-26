# use boltz_analysis_env 
import requests
import math
import sys
import os
sys.path.append(os.path.expandvars("${TRIAGGED}/covalent_module/analysis"))
from analyize_boltz_preds import mean_metrics
from analysis_utils import read_boltz_predictions
import numpy as np
import pandas as pd

def cdd_query(API_KEY, VAULT_ID, readout_query={}, mol_query={}):
    '''
    Queries CDD Vault and returns query results.
    Eg. 
    API_KEY = "NTc3Mnw1WFFxeC9mOE1SZ09NSnZjY2VmWVFWVFpIUGhJUUJLWUVrcG1rOEM2WnRWTkZQK0pSQT09"
    VAULT_ID = 7171 
    
    params = {
    "modified_after": "2023-09-01T00:00:00Z",  # Filter molecules created after sept 2025
    "page_size": 1000,
    "protocols": 89683
    }
    params = {
    "modified_after": "2023-09-01T00:00:00Z",  # Filter molecules created after sept 2025
    "page_size": 1000
    }
    '''
    headers = {"X-CDD-Token": API_KEY}

    readout_rows_url = f"https://app.collaborativedrug.com/api/v1/vaults/{VAULT_ID}/readout_rows"
    molecules_url = f"https://app.collaborativedrug.com/api/v1/vaults/{VAULT_ID}/molecules"

    readout_response = requests.get(readout_rows_url, headers=headers, params=readout_query)
    readout_response.raise_for_status()
    readouts = readout_response.json()

    molecules_response = requests.get(molecules_url, headers=headers, params=mol_query)
    molecules_response.raise_for_status()
    molecules = molecules_response.json()

    return readouts, molecules

def reorg_query(readouts, molecules, filter=False):
    '''
    Combines, matches and reorganizes readouts and molecule IDs from CDD vault. 
    
    readouts (dict): json response from CDD vault query for readout rows. has 'objects' attribute. 
    eg. {'count': 1161,
        'offset': 0,
        'page_size': 1000,
        'objects': [{'id': 1198645303,...}]}
    molecules (dict): json reponse from CDD vault query for molecules in vault. has 'objects' attribute. 
    eg. {'count': 353,
        'offset': 0,
        'page_size': 1000,
        'objects': [{'id': 159552986, ...}]} 
    filter (bool): if True, will filter out molecules without 'TJB' in synonym. 

    returns (dict): molecule ids matched to ic50s (nM) for hscpl and tgcpl. includes log(ic50 uM) converted readouts. 
    eg. {'VMCC-0000135': {'tgcpl_ic50': 2500.0,
        'tgcpl_ic50_log': 0.3979400086720376,
        'hscpl_ic50': 2300.0,
        'hscpl_ic50_log': 0.36172783601759284,
        'smiles': 'CCNc1nc(C#N)nc(N2CCOCC2)n1'},
        'VMCC-0000147': {'tgcpl_ic50': 28.8433333333333...}}
    '''
    
    molecule_ids = {}
    no_tjb = []
    for row in molecules['objects']:
        if not filter: 
            molecule_ids[row['id']] = {row['name']: row['smiles']}

        if filter:
            for syn_name in row.get('synonyms', []):
                if 'TJB' in syn_name:
                    molecule_ids[row['id']] = {row['name']: row['smiles']}
                    break
                else:
                    no_tjb.append(row['name'])

    if filter: print("no TJB in synonym name:", no_tjb)

    readouts_and_ids = {}
    for row in readouts["objects"]:
        mol_read = {}

        for readout_id, readout_data in row["readouts"].items():
            readout_val = readout_data.get("value")

            # Map readout ID → experiment name
            if readout_id == "1125023":
                exp_name = "tgcpl_ic50"
            elif readout_id == "1125024":
                exp_name = "hscpl_ic50"
            else:
                continue  # skip unknown readouts

            # Store raw IC50 (nM)
            mol_read[exp_name] = readout_val
            if readout_val and readout_val > 0:
                mol_read[f"{exp_name}_log"] = math.log10(readout_val / 1000)
            else:
                mol_read[f"{exp_name}_log"] = None

        try:
            mol_name = next(iter(molecule_ids[row["molecule"]]))
            smiles = molecule_ids[row["molecule"]][mol_name]

            mol_read["smiles"] = smiles
            readouts_and_ids[mol_name] = mol_read

        except KeyError:
            pass

    return readouts_and_ids

# updating comound_records 
def append_comment(df, idx, col, msg, add_comment=False):
    '''Adds comments to the spread sheet based on whether a 
    prediction and experimental readout or not was present.
    '''
    if add_comment:
        current = df.at[idx, col]
        if not isinstance(current, str):
            current = ''
        df.at[idx, col] = current + msg

def lookup_compound_id(vault_id, compound_record):
    '''Matches vault_id from CDD vault with 5-char unique .pkl identifying key'''
    if 'vault_id' not in compound_record.columns or 'compound_id' not in compound_record.columns:
        print('[ERROR] add vault_id and compound_id columns')
        return None
    match = compound_record.loc[compound_record['vault_id'] == vault_id, 'compound_id']
    compound_id = match.iloc[0] if not match.empty else None
    return compound_id

# need to pull experimental info from cdd vault 
def add_experiment_reads(cdd_query, compound_rec_df): 
    '''
        Adding missing TgCPL/HsCPL info to records.
        NaN/Empty values if data not present. 
        
        cdd_query (dict): organized query results 
        eg. {'VMCC-0000135': {'tgcpl_ic50': 2500.0,
                'tgcpl_ic50_log': 0.3979400086720376,
                'hscpl_ic50': 2300.0,
                'hscpl_ic50_log': 0.36172783601759284,
                'smiles': 'CCNc1nc(C#N)nc(N2CCOCC2)n1'}, ...}
    '''
    # update regardless of whether the value is there or not, based on query date
    counter = 0
    for row in compound_rec_df.itertuples(index=True):
        vault_id = row.vault_id
        if vault_id not in cdd_query:
            append_comment(
                compound_rec_df,
                row.Index,
                'comments',
                'IC50 not in TJB TgCPL set, IC50 not in TJB HsCPL set,'
            )
            continue

        entry = cdd_query[vault_id]
        if 'tgcpl_ic50_log' in entry:
            compound_rec_df.at[row.Index, 'tgcpl_log_ic50'] = entry['tgcpl_ic50_log']
            compound_rec_df.at[row.Index, 'tgcpl_ic50'] = entry['tgcpl_ic50']
            counter += 1
        else:
            append_comment(
                compound_rec_df,
                row.Index,
                'comments',
                'IC50 not in TJB TgCPL set,'
            )

        if 'hscpl_ic50_log' in entry:
            compound_rec_df.at[row.Index, 'hscpl_log_ic50'] = entry['hscpl_ic50_log']
            compound_rec_df.at[row.Index, 'hscpl_ic50'] = entry['hscpl_ic50']
            counter += 1
        else:
            append_comment(
                compound_rec_df,
                row.Index,
                'comments',
                ' IC50 not in TJB HsCPL set,'
            )
            
    print(counter, " experimental readouts updated")

    return compound_rec_df

def add_smiles(cdd_query, compound_rec_df):
    '''
    Adds smiles from ccd query. 
    '''
    for row in compound_rec_df.itertuples(index=True):
        if pd.isna(row.smiles) and row.vault_id in cdd_query:
            compound_rec_df.at[row.Index, 'smiles'] = cdd_query[row.vault_id]['smiles']
    return compound_rec_df

def add_predictions(predictions_dir, compound_rec_df, target):
    '''
    Adding missing TgCPL/HsCPL predicted ic50s to records.
    '''
    all_reps = [
                os.path.join(predictions_dir, f)
                for f in os.listdir(predictions_dir)
                if os.path.isdir(os.path.join(predictions_dir, f))
            ]
    
    if 'tgcpl' in target:
        counter = 0
        for row in compound_rec_df.itertuples(index=True):
            compound_id = lookup_compound_id(row.vault_id, compound_rec_df) 
            if np.isnan(row.tgcpl_pred_log_ic50) and compound_id is not None and len(all_reps) > 1:
                _, mean_preds  = mean_metrics(all_reps[0], score_col='Pred log10(IC50)')
                try:
                    pred_ic50 = mean_preds.loc[mean_preds['compound_id'] == compound_id, 'mean'].values[0]
                    compound_rec_df.at[row.Index, 'tgcpl_pred_log_ic50'] = pred_ic50
                    counter += 1
                except IndexError:
                    pred_ic50 = None
                    append_comment(compound_rec_df, row.Index, 'comments', ' missing in TgCPL Pred,')

            elif np.isnan(row.tgcpl_pred_log_ic50) and compound_id is not None and len(all_reps) == 1: # only 1 rep present
                preds = read_boltz_predictions(all_reps[0], reps=False)
                try:
                    pred_ic50 = preds.loc[preds['compound_id'] == compound_id, 'Pred log10(IC50)'].values[0]
                    compound_rec_df.at[row.Index, 'tgcpl_pred_log_ic50'] = pred_ic50
                    counter += 1
                except IndexError:
                    pred_ic50 = None
                    append_comment(compound_rec_df, row.Index, 'comments', ' missing in TgCPL Pred,')
            else:
                append_comment(compound_rec_df, row.Index, 'comments', ' missing in TgCPL Pred,')

        print(f"{counter} predictions were updated for tgcpl")  

    elif 'hscpl' in target:
        counter = 0 
        for row in compound_rec_df.itertuples(index=True):
            compound_id = lookup_compound_id(row.vault_id, compound_rec_df) 
            if np.isnan(row.hscpl_pred_log_ic50) and compound_id is not None and len(all_reps) > 1:  # check if more than 1 rep present 
                    _, mean_preds  = mean_metrics(all_reps[0], score_col='Pred log10(IC50)')
                    try:
                        pred_ic50 = mean_preds.loc[mean_preds['compound_id'] == compound_id, 'mean'].values[0]
                        compound_rec_df.at[row.Index, 'hscpl_pred_log_ic50'] = pred_ic50
                        counter +=1 
                    except IndexError:
                        pred_ic50 = None
                        append_comment(compound_rec_df, row.Index, 'comments', ' missing in HsCPL Pred,')

            elif np.isnan(row.tgcpl_pred_log_ic50) and compound_id is not None and len(all_reps) == 1: # only 1 rep present
                preds = read_boltz_predictions(all_reps[0], reps=False)
                try:
                    pred_ic50 = preds.loc[preds['compound_id'] == compound_id, 'Pred log10(IC50)'].values[0]
                    compound_rec_df.at[row.Index, 'hscpl_pred_log_ic50'] = pred_ic50
                    counter+=1 
                except IndexError:
                    pred_ic50 = None
                    append_comment(compound_rec_df, row.Index, 'comments', ' missing in HsCPL Pred,')

            else:
                append_comment(compound_rec_df, row.Index, 'comments', ' missing in HsCPL Pred,')   

        print(f"{counter} predictions were updated for hscpl")  

    return compound_rec_df


