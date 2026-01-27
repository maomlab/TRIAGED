# use boltz_analysis_env 
import requests
import time
import requests
import math 

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

    readout_query_clean = {k: v for k, v in readout_query.items() if k != 'page_size'}
    mol_query_clean = {k: v for k, v in mol_query.items() if k != 'page_size'}
    
    readout_query_clean["async"] = "true"
    mol_query_clean["async"] = "true"

    readout_response = requests.get(readout_rows_url, headers=headers, params=readout_query_clean)
    readout_response.raise_for_status()
    readout_export = readout_response.json()
    
    molecules_response = requests.get(molecules_url, headers=headers, params=mol_query_clean)
    molecules_response.raise_for_status()
    molecules_export = molecules_response.json()

    readout_export_id = readout_export.get("id")
    molecules_export_id = molecules_export.get("id")
    
    if not readout_export_id or not molecules_export_id:
        raise Exception("Failed to get export IDs")

    def get_export_data(export_id):
        export_url = f"https://app.collaborativedrug.com/api/v1/vaults/{VAULT_ID}/exports/{export_id}"
        max_attempts = 120
        attempt = 0
        
        while attempt < max_attempts:
            time.sleep(5)
            attempt += 1
            
            response = requests.get(export_url, headers=headers)
            response.raise_for_status()
            data = response.json()
            
            if "objects" in data:
                return data
            
            status = data.get("status")
            
            if status == "finished":
                return data
            elif status in ["new", "started", "pending", None]:
                continue
            elif status == "failed":
                raise Exception(f"Export failed: {data.get('error', 'Unknown')}")
            else:
                print(f"Debug - status: {status}, full response: {data}")
                continue
        
        raise Exception("Export timed out")

    readouts = get_export_data(readout_export_id)
    molecules = get_export_data(molecules_export_id)
    
    return readouts, molecules

def get_ic50s(readouts, molecules):
    '''
    Loads MEAN ic50 values for both TgCPL and HsCPL.
    Uses two DataFrames where 'molecule' from readouts matches in 'id' in molecules and
    Returns experimental data and metadata of compounds from user CDD-Vault. 
    '''
    metadata_df = []
    exp_readouts = []
    for row in molecules['objects']: # per molecule 
        molecule_id = row['id']
        substance_id = row['name']
        inchi_key = row['inchi_key']
        inchi = row['inchi']
        smiles = row['smiles']
        syn = row["synonyms"]
        syn_str = (',').join(syn)

        # one per molecule
        metadata_df.append({    
                    "vault_mol_id" : molecule_id,
                    "substance_id" : substance_id,
                    "inchi_key": inchi_key,
                    "inchi": inchi, 
                    "smiles": smiles,
                    "synonyms": syn_str
                })
        
        # must find molecule match in readout
        tgcpl_log_ic50 = None
        hscpl_log_ic50 = None
        for r in readouts['objects']: 
            if r['molecule'] == str(molecule_id) or r['molecule'] == int(molecule_id): 
                for key, val in r['readouts'].items():
                    # print(molecule_id, key, val['value'])
                    if str(key) == "1125023":
                        tgcpl_mean_ic50 = val['value'] 
                        tgcpl_log_ic50 = math.log10(tgcpl_mean_ic50 / 1000)
                    elif str(key) == "1125024": 
                        hscpl_mean_ic50 = val['value']
                        hscpl_log_ic50 = math.log10(hscpl_mean_ic50 / 1000)
                    else: 
                        continue

        exp_readouts.append({
            "vault_mol_id" : molecule_id,
            "substance_id": row['name'],
            "inchi_key": row['inchi_key'],
            "mean_tgcpl_log_ic50 (uM)": tgcpl_log_ic50,
            "mean_hscpl_log_ic50 (uM)": hscpl_log_ic50 # to do: pulling reps info
        })
    return exp_readouts, metadata_df

def update_local_data(old_metadata, old_readouts, new_metadata, new_readouts):
    '''
    All inputs must be Pandas.DataFrames. 
    Will overwrite existing CSVs!!! 
    '''
    # Metadata: update existing and add new
    merged_metadata = old_metadata.merge(new_metadata, on='vault_mol_id', how='outer', suffixes=('_old', '_new'))

    for col in old_metadata.columns:
        if col == 'vault_mol_id':
            continue
        if f'{col}_new' in merged_metadata.columns:
            merged_metadata[col] = merged_metadata[f'{col}_new'].fillna(merged_metadata[f'{col}_old'])
            merged_metadata.drop([f'{col}_old', f'{col}_new'], axis=1, inplace=True)

    # Readouts: update existing and add new
    merged_readouts = old_readouts.merge(new_readouts, on='vault_mol_id', how='outer', suffixes=('_old', '_new'))

    for col in old_readouts.columns:
        if col == 'vault_mol_id':
            continue
        if f'{col}_new' in merged_readouts.columns:
            merged_readouts[col] = merged_readouts[f'{col}_new'].fillna(merged_readouts[f'{col}_old'])
            merged_readouts.drop([f'{col}_old', f'{col}_new'], axis=1, inplace=True)

    return merged_metadata, merged_readouts
