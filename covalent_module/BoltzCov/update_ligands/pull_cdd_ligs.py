# use boltz_analysis_env 
import requests
import time
import requests
import math 
import pandas as pd 

def cdd_query(API_KEY, VAULT_ID, readout_query={}, mol_query={}):
    '''
    Queries CDD Vault and returns query results.
    Supports both sync (small queries) and async (large queries) modes.
    
    Eg. 
    API_KEY = "NTc3Mnw1WFFxeC9mOE1SZ09NSnZjY2VmWVFWVFpIUGhJUUJLWUVrcG1rOEM2WnRWTkZQK0pSQT09"
    VAULT_ID = 7171 
    
    readout_query = {
        "modified_after": "2023-09-01T00:00:00Z",
        "page_size": 10,
        "protocols": 89683
    }
    mol_query = {
        "modified_after": "2023-09-01T00:00:00Z",
        "page_size": 10
    }
    '''
    headers = {"X-CDD-Token": API_KEY}

    readout_rows_url = f"https://app.collaborativedrug.com/api/v1/vaults/{VAULT_ID}/readout_rows"
    molecules_url = f"https://app.collaborativedrug.com/api/v1/vaults/{VAULT_ID}/molecules"

    # Check if page_size > 1000, use async mode
    use_async = readout_query.get('page_size', 0) >= 1000 or mol_query.get('page_size', 0) >= 1000
    
    if use_async:
        readout_query_clean = {k: v for k, v in readout_query.items() if k != 'page_size'}
        mol_query_clean = {k: v for k, v in mol_query.items() if k != 'page_size'}
        readout_query_clean["async"] = "true"
        mol_query_clean["async"] = "true"
    else:
        readout_query_clean = readout_query.copy()
        mol_query_clean = mol_query.copy()

    readout_response = requests.get(readout_rows_url, headers=headers, params=readout_query_clean)
    readout_response.raise_for_status()
    readout_export = readout_response.json()
    
    molecules_response = requests.get(molecules_url, headers=headers, params=mol_query_clean)
    molecules_response.raise_for_status()
    molecules_export = molecules_response.json()

    if use_async:
        readout_export_id = readout_export.get("id")
        molecules_export_id = molecules_export.get("id")
        
        if not readout_export_id or not molecules_export_id:
            raise Exception("Failed to get export IDs")

        def get_export_data(export_id):
            export_url = f"https://app.collaborativedrug.com/api/v1/vaults/{VAULT_ID}/exports/{export_id}"
            max_attempts = 120
            attempt = 0
            time.sleep(3)  # wait before first poll
            while attempt < max_attempts:
                time.sleep(5)
                attempt += 1
                
                response = requests.get(export_url, headers=headers)
                if response.status_code == 403:
                    continue
                response.raise_for_status()
                data = response.json()
                
                if "objects" in data:
                    return data
                
                status = data.get("status")
                
                if status == "finished":
                    response = requests.get(export_url, headers=headers)
                    return response.json()
                elif status in ["new", "started", "pending", None]:
                    print(f"Attempt {attempt}: status={status}")
                    continue
                elif status == "failed":
                    raise Exception(f"Export failed: {data.get('error', 'Unknown')}")
                else:
                    print(f"Debug - status: {status}, full response: {data}")
                    continue
            
            raise Exception("Export timed out")

        readouts = get_export_data(readout_export_id)
        molecules = get_export_data(molecules_export_id)
    else:
        readouts = readout_export 
        molecules = molecules_export 
    
    return readouts, molecules

def get_ic50s(readouts, molecules):
    '''
    Loads MEAN ic50 values for both TgCPL and HsCPL.
    Uses two DataFrames where 'molecule' from readouts matches in 'id' in molecules and
    Returns experimental data and metadata of compounds from user CDD-Vault. 
    Now also extracts all numerical properties from molecules.
    '''
    metadata_list = []
    exp_readouts_list = []
    
    # Fields to exclude
    excluded_fields = {'registration_form_id', 'exact_mass', 'cdd_registry_number'}
    
    for row in molecules['objects']: # per molecule 
        molecule_id = row['id']
        substance_id = row['name']
        inchi_key = row['inchi_key']
        inchi = row['inchi']
        smiles = row['smiles']
        syn = row["synonyms"]
        syn_str = (',').join(syn)

        # Start with basic metadata
        metadata_entry = {    
            "vault_mol_id": molecule_id,
            "substance_id": substance_id,
            "inchi_key": inchi_key,
            "inchi": inchi, 
            "smiles": smiles,
            "synonyms": syn_str
        }
        
        # Extract molecular weight if present
        if 'molecular_weight' in row and row['molecular_weight'] is not None:
            metadata_entry['molecular_weight'] = float(row['molecular_weight'])
        
        # Extract custom fields (numerical properties)
        if 'fields' in row and isinstance(row['fields'], dict):
            for field_id, field_data in row['fields'].items():
                if isinstance(field_data, dict) and 'value' in field_data:
                    value = field_data['value']
                    # Add numerical values
                    if isinstance(value, (int, float)) and value is not None:
                        # Use field name if available, otherwise use field_id
                        field_name = field_data.get('name', f'field_{field_id}')
                        # Clean field name for use as column name
                        field_name_clean = field_name.replace(' ', '_').replace('(', '').replace(')', '').lower()
                        
                        # Skip excluded fields
                        if field_name_clean not in excluded_fields and field_name.lower() not in excluded_fields:
                            metadata_entry[field_name_clean] = value
        
        # Extract any other top-level numerical fields
        for key, value in row.items():
            # Skip non-numerical or already processed fields
            if key in ['id', 'name', 'inchi_key', 'inchi', 'smiles', 'synonyms', 
                      'created_at', 'modified_at', 'projects', 'batches', 'fields',
                      'molecular_weight']:  # Already handled
                continue
            
            # Skip excluded fields
            if key.lower() in excluded_fields:
                continue
            
            # Add numerical values (int, float)
            if isinstance(value, (int, float)) and value is not None:
                metadata_entry[key] = value
        
        metadata_list.append(metadata_entry)
        
        # must find molecule match in readout
        tgcpl_log_ic50 = None
        hscpl_log_ic50 = None
        for r in readouts['objects']: 
            if str(r['molecule']) == str(molecule_id):
                print(f"Molecule {molecule_id} readout keys: {list(r['readouts'].keys())}")
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

        exp_readouts_list.append({
            "vault_mol_id": molecule_id,
            "substance_id": row['name'],
            "inchi_key": row['inchi_key'],
            "mean_tgcpl_log_ic50 (uM)": tgcpl_log_ic50,
            "mean_hscpl_log_ic50 (uM)": hscpl_log_ic50 # to do: pulling reps info
        })
    
    return pd.DataFrame(exp_readouts_list), pd.DataFrame(metadata_list)

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
            merged_readouts[col] = merged_readouts[f'{col}_new'].fillna(merged_readouts[f'{col}_old']).infer_objects(copy=False)
            merged_readouts.drop([f'{col}_old', f'{col}_new'], axis=1, inplace=True)

    return merged_metadata, merged_readouts
