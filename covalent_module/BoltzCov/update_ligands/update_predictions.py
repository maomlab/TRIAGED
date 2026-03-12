# use boltz_analysis_env 

import pandas as pd
import os 
import glob     
import shutil
import json
import math 

def fetch_new(pred_df, metadata_df, protein, min_replicates=3):
    '''
    Checks if predictions of compounds are in predictions.csv records.
    Works with long-format predictions.csv (with 'protein' column).
    Does NOT check if .cif file of final prediction is present.
    
    pred_df: df with names of molecules, confidence values and predictions (long format)
    metadata_df: df with all compound information
    protein: any name for the protein. taken from pdb file name.
    min_replicates: minimum number of non-NaN predictions required (default: 3)
    Returns DataFrame of new compounds that have no predictions yet OR have fewer than min_replicates. 
    '''
    
    # Filter predictions for this specific protein
    prot_preds = pred_df[pred_df['protein'] == protein].copy()
    pred_col = 'pred_log10ic50'  # Column name without protein suffix
    
    # If pred_col doesn't exist, return all compounds as new
    if pred_col not in prot_preds.columns:
        return metadata_df[['substance_id', 'inchi_key', 'smiles']].copy()
    
    # Count non-NaN predictions per compound for this protein
    pred_counts = prot_preds.groupby(['substance_id', 'inchi_key'])[pred_col].apply(
        lambda x: x.notna().sum()
    ).reset_index(name='pred_count')
    
    # Merge metadata with prediction counts
    merged = metadata_df.merge(
        pred_counts,
        on=['substance_id', 'inchi_key'], 
        how='left'
    )
    
    # Fill NaN counts with 0 (compounds with no predictions)
    merged['pred_count'] = merged['pred_count'].fillna(0)
    
    # Keep compounds with no predictions OR fewer than min_replicates
    no_pred_prot = merged[merged['pred_count'] <= min_replicates][['substance_id', 'inchi_key', 'smiles']]
    
    return no_pred_prot

def check_attempted(error_df, dock_compounds):
    '''
    Checks errored.csv records to skip attempted ligands that failed in past runs.
    Returns updated DataFrame with errored compounds removed.
    '''
    # Anti-join: keep only rows NOT in error_df
    dock_compounds = dock_compounds.merge(
        error_df[['substance_id', 'inchi_key']], 
        on=['substance_id', 'inchi_key'], 
        how='left', 
        indicator=True
    )
    dock_compounds = dock_compounds[dock_compounds['_merge'] == 'left_only']
    dock_compounds = dock_compounds.drop('_merge', axis=1)
    
    return dock_compounds

def get_identifier(pkl_id, boltz_cache_prot, record_path):
    '''
    Gets unique identifier information about molecule.
    
    :param pkl_id: ID generated during yaml and pkl creation before docking. 
    :param boltz_cache_prot: Path to raw predictions for a given protein.
    :param record_path: Path to metadata.csv record.
    '''

    metadata_record = os.path.join(record_path, 'metadata.csv')
    metadata_df = pd.read_csv(metadata_record)
    csv_files = glob.glob(os.path.join(boltz_cache_prot, '*.csv'))
    if csv_files: ligand_cache_csv = csv_files[0] # temp csv that was generating before boltz run to make yamls 
    ligand_csv_df = pd.read_csv(ligand_cache_csv)

    # Get compound info
    if len(pkl_id) <= 5: 
        substance_id = ligand_csv_df[ligand_csv_df['pkl_id'] == pkl_id]['substance_id'].iloc[0]
    else:
        substance_id = pkl_id
    
    meta_row = metadata_df[metadata_df['substance_id'] == substance_id].iloc[0]
    inchi_key = meta_row['inchi_key']
    vault_mol_id = meta_row['vault_mol_id']

    return substance_id, inchi_key, vault_mol_id

def convert_IC_to_energy(IC):
    """
    Convert log10(IC50) measured in uM to kcal/mol.
    :param IC: IC50.
    :return: converted kcal/mol estimate.
    """
    return (6 - IC) * 1.364

def read_boltz_predictions(predictions_dir):
    """
    Reads prediction JSON files from subdirectories and compiles them into a pandas DataFrame.
    :param predictions_dir: Path to the directory containing subdirectories with JSON files.
    :return: Pandas DataFrame with compiled data.
    """
    data = []
    list_compound_dirs = [d.name for d in os.scandir(predictions_dir) if d.is_dir()]
    for compound_name in list_compound_dirs:
        compound_id = compound_name.split('_')[-1]
        compound_dir = os.path.join(predictions_dir, compound_name)
        results = [compound_dir,  f"boltz_results_{compound_id}", "predictions", f"{compound_id}"]
        compound_result =  "/".join(results)
        if not os.path.isdir(compound_result):
            continue

        affinity_file = os.path.join(compound_result, f"affinity_{compound_id}.json")
        confidence_file = os.path.join(compound_result, f"confidence_{compound_id}_model_0.json")

        if os.path.exists(affinity_file) and os.path.exists(confidence_file):
            with open(affinity_file, 'r') as af:
                affinity_data = json.load(af)
                affinity_pred_value = affinity_data.get("affinity_pred_value", None)
                ic50_nm = (10 ** affinity_pred_value) * 1000 
                pred_pic50 = -math.log10((10 ** affinity_pred_value) * 1e-6)
                affinity_probability_binary = affinity_data.get("affinity_probability_binary", None)

            with open(confidence_file, 'r') as cf:
                confidence_data = json.load(cf)
                confidence_score = confidence_data.get("confidence_score", None)
                ptm = confidence_data.get("ptm", None)
                iptm = confidence_data.get("iptm", None)
                ligand_iptm = confidence_data.get("ligand_iptm", None)
                protein_iptm = confidence_data.get("protein_iptm", None)
                complex_plddt = confidence_data.get("complex_plddt", None)
                complex_iplddt = confidence_data.get("complex_iplddt", None)
                complex_pde = confidence_data.get("complex_pde", None)
                complex_ipde = confidence_data.get("complex_ipde", None)

            energy_value = convert_IC_to_energy(affinity_pred_value) if affinity_pred_value is not None else None

            data.append({
                "compound_id": compound_id,
                "Pred log10(IC50)": affinity_pred_value,
                "Pred pIC50": pred_pic50,
                "Pred Label (IC50-like)": True if ic50_nm < 1000 else False,
                "Binding Probability": affinity_probability_binary,
                "Pred Label":  True if affinity_probability_binary > 0.5 else False,
                "Confidence Score": confidence_score,
                "kcal/mol": energy_value,
                "PTM": ptm,
                "IPTM": iptm,
                "Ligand IPTM": ligand_iptm,
                "Protein IPTM": protein_iptm,
                "Complex pLDDT": complex_plddt,
                "Complex iPLDDT": complex_iplddt,
                "Complex PDE": complex_pde,
                "Complex iPDE": complex_ipde
            })
    return pd.DataFrame(data)

def update_preds(pkl_id, boltz_cache_prot, record_path, protein_name, VERBOSE=False):
    ''' 
    Update predictions in predictions.csv. Each run creates a unique row (no averaging).
    Now saves data in long format with a 'protein' column.
    '''
    boltz_preds = read_boltz_predictions(predictions_dir=boltz_cache_prot)
    
    new_pred = boltz_preds[boltz_preds['compound_id'] == pkl_id]
    if new_pred.empty:
        return
    
    new_pred = new_pred.iloc[0]
    
    substance_id, inchi_key, vault_mol_id = get_identifier(pkl_id, boltz_cache_prot, record_path)
    
    pred_path = glob.glob(f'{boltz_cache_prot}/*/boltz_results_*/predictions/{pkl_id}*')
    if pred_path:
        boltz_runID = pred_path[0].split('/')[-4]
    else:
        boltz_runID = None
    
    if len(pkl_id) <= 5: # covalent
        predictions_csv = os.path.join(record_path, 'predictions.csv')
    else: 
        predictions_csv = os.path.join(record_path, 'noncov_predictions.csv')

    numeric_cols = ['Pred log10(IC50)', 'Pred pIC50', 'Binding Probability', 'Confidence Score',
                    'PTM', 'IPTM', 'Ligand IPTM', 'Protein IPTM', 'Complex pLDDT', 
                    'Complex iPLDDT', 'Complex PDE', 'Complex iPDE']
    
    # Create new row WITHOUT protein suffix on metric columns
    # Add 'protein' as a separate column instead
    prediction_data = pd.DataFrame([{
        'substance_id': substance_id,
        'inchi_key': inchi_key,
        'vault_mol_id': vault_mol_id,
        'boltz_runID': boltz_runID,
        'protein': protein_name,  # Add protein as a separate column
        **{col.lower().replace(" ", "_").replace("(", "").replace(")", ""): new_pred[col] 
           for col in numeric_cols if col in new_pred.index}
    }])
    
    # Append new row (no averaging, boltz_runID makes each row unique)
    if os.path.exists(predictions_csv):
        old_preds = pd.read_csv(predictions_csv)
        if not bool((old_preds['boltz_runID'] == boltz_runID).any()):
            updated_preds = pd.concat([old_preds, prediction_data], ignore_index=True)
            updated_preds.to_csv(predictions_csv, index=False)
            if VERBOSE: print(f'Updated predictions for {substance_id}, {boltz_runID}')
        else: 
            if VERBOSE: print(f'Prediction present for {substance_id}, {boltz_runID}')
    else:
        prediction_data.to_csv(predictions_csv, index=False)

def check_pred(run_cache_prot, record_path, VERBOSE, protein_name, COVALENT=True):
    '''
    Given the path to raw predictions, check if .cif was produced or not. 
    Updates errored.csv and predictions.csv.
    '''
    error_record = os.path.join(record_path, 'errored.csv')

    csv_files = glob.glob(os.path.join(run_cache_prot, '*.csv'))
    if csv_files: 
        ligand_cache_csv = csv_files[0]
    else:
        print(f"No ligand CSV found in {run_cache_prot}")
        return
        
    ligand_csv_df = pd.read_csv(ligand_cache_csv)

    if COVALENT:
        for pkl_id in ligand_csv_df['pkl_id']:
            cif_files = glob.glob(f'{run_cache_prot}/*/boltz_results_*/predictions/{pkl_id}*/*.cif')
            json_files = glob.glob(f'{run_cache_prot}/*/boltz_results_*/predictions/{pkl_id}*/*.json')
            has_both = len(cif_files) > 0 and len(json_files) > 0

            # Get compound info
            substance_id, inchi_key, vault_mol_id = get_identifier(pkl_id, run_cache_prot, record_path)

            if has_both:
                update_preds(pkl_id, run_cache_prot, record_path, protein_name, VERBOSE)  
            else: 
                error_data = pd.DataFrame([{
                    'substance_id': substance_id,
                    'inchi_key': inchi_key,
                    'vault_mol_id': vault_mol_id
                }])
                if os.path.exists(error_record): 
                    error_data.to_csv(error_record, mode='a', header=False, index=False)
                else: 
                    error_data.to_csv(error_record, mode='w', header=True, index=False)
                if VERBOSE: print(f"Added {substance_id} to errored.csv")
    else: # noncovalent
        for substance_id in ligand_csv_df['pkl_id']:
            cif_files = glob.glob(f'{run_cache_prot}/*/boltz_results_*/predictions/{substance_id}*/*.cif')
            json_files = glob.glob(f'{run_cache_prot}/*/boltz_results_*/predictions/{substance_id}*/*.json')
            has_both = len(cif_files) > 0 and len(json_files) > 0

            if has_both:
                update_preds(substance_id, run_cache_prot, record_path, protein_name, VERBOSE)  
            else: 
                error_data = pd.DataFrame([{
                    'substance_id': substance_id,
                    'inchi_key': inchi_key,
                    'vault_mol_id': vault_mol_id
                }])
                if os.path.exists(error_record): 
                    error_data.to_csv(error_record, mode='a', header=False, index=False)
                else: 
                    error_data.to_csv(error_record, mode='w', header=True, index=False)
                if VERBOSE: print(f"Added {substance_id} to errored.csv")
    print("Record updates complete.")

def reorg_preds(run_cache_prot, ligand_df, output_temp, VERBOSE):
    '''
    Copy .cif and .yaml files to output directory for successful predictions.
    Files are renamed to boltz_runID_pklID format.
    output_temp: /run_cache_prot/temp/protein_name
    '''
    # Create output subdirectories
    cif_dir = os.path.join(output_temp, 'cifs') # /run_cache_prot/temp/protein_name/cifs
    yaml_dir = os.path.join(output_temp, 'yamls') # /run_cache_prot/temp/protein_name/yamls
    os.makedirs(cif_dir, exist_ok=True)
    os.makedirs(yaml_dir, exist_ok=True)

    for _, row in ligand_df.iterrows():
        pkl_id = row['pkl_id']
        substance_id = row.get('substance_id', pkl_id)

        pred_path = glob.glob(f'{run_cache_prot}/*/boltz_results_*/predictions/{pkl_id}*')
    
        if pred_path:
            boltz_runID = pred_path[0].split('/')[-4]
        
        # Filter cif files specific to this compound
        cif_files = glob.glob(f'{run_cache_prot}/*/boltz_results_*/predictions/*/*.cif')
        if not cif_files:
            # try by substance_id 
            cif_files = glob.glob(f'{run_cache_prot}/*/boltz_results_*/predictions/*{substance_id}*/*.cif')
        
        if cif_files:
            cif_file = cif_files[0]
            # Copy .cif file with new name
            cif_dest = os.path.join(cif_dir, f"{boltz_runID}_{substance_id}.cif")
            if substance_id == boltz_runID.split('_')[-1]: 
                cif_dest = os.path.join(cif_dir, f"{boltz_runID}.cif")
            shutil.copy(cif_file, cif_dest)
            if VERBOSE: print(f"Copied {substance_id} cif for {boltz_runID}")
            
            # Find and copy yaml file
            # Extract parent directory: run_cache_prot/YYMMDD_HHMMSS_pklID/
            parent_dir = cif_file.split('/boltz_results_')[0]
            yaml_files = glob.glob(f'{parent_dir}/*.yaml')
            if yaml_files:
                yaml_file = yaml_files[0]
                yaml_dest = os.path.join(yaml_dir, f"{boltz_runID}_{substance_id}.yaml")
                if substance_id == boltz_runID.split('_')[-1]: 
                    yaml_dest = os.path.join(yaml_dir, f"{boltz_runID}.yaml")
                shutil.copy(yaml_file, yaml_dest)
                if VERBOSE: print(f"Copied {substance_id} yaml for {boltz_runID}")
        else:
            print(f"[WARNING] No .cif file found for {substance_id}")
    if VERBOSE: print(f"Files copied to {output_temp}.")

def reorg_reps(output_temp, output_dir, protein_name):
    '''
    Reorgs predictions from temp into rep directories in output_dir.
    output_temp: /run_cache_prot/temp/protein_name
    output_dir:  /output/user_given_output_name 
    Files in output_temp are named: datetime_substanceid.cif/yaml
    Creates new rep directories dynamically if needed.
    '''
    output_parent = os.path.dirname(output_dir) # may have tgcpl or hscpl with rep1, rep2.. 
    output_parent_prot = os.path.join(output_parent, protein_name) 

    existing_reps = sorted(glob.glob(f'{output_parent_prot}/rep*'))
    
    # Create rep1 if no reps exist
    if not existing_reps:
        rep1 = os.path.join(output_dir, 'rep1')
        os.makedirs(os.path.join(rep1, 'cifs'), exist_ok=True)
        os.makedirs(os.path.join(rep1, 'yamls'), exist_ok=True)
        existing_reps = [rep1]
        print(f"Created {rep1}")

    temp_cifs = glob.glob(f'{output_temp}/temp/*/cifs/*.cif')
    if not temp_cifs:
        print(f"[WARNING] No cif files found in {output_temp}/cifs/")
        return

    for cif in temp_cifs:
        cif_name = os.path.basename(cif)
        substance_id = cif_name.split('_')[-1].split('.')[0]  # VMCC-0000740

        # Find first rep that doesn't already have this ligand
        placed = False
        for rep_dir in existing_reps:
            existing = glob.glob(f'{rep_dir}/cifs/*{substance_id}*')
            if not existing:
                shutil.copy(cif, os.path.join(rep_dir, 'cifs', cif_name))
                print(f"Copied {cif_name} → {rep_dir}/cifs/")

                yaml_path = cif.replace('/cifs/', '/yamls/').replace('.cif', '.yaml')
                if os.path.exists(yaml_path):
                    shutil.copy(yaml_path, os.path.join(rep_dir, 'yamls', os.path.basename(yaml_path)))
                    print(f"Copied {os.path.basename(yaml_path)} → {rep_dir}/yamls/")
                else:
                    print(f"[WARNING] No yaml found for {cif_name}")
                placed = True
                break

        # If all reps already have this ligand, create a new rep
        if not placed:
            new_rep = os.path.join(output_dir, f'rep{len(existing_reps)+1}')
            os.makedirs(os.path.join(new_rep, 'cifs'), exist_ok=True)
            os.makedirs(os.path.join(new_rep, 'yamls'), exist_ok=True)
            existing_reps.append(new_rep)
            print(f"Created {new_rep}")

            shutil.copy(cif, os.path.join(new_rep, 'cifs', cif_name))
            print(f"Copied {cif_name} → {new_rep}/cifs/")

            yaml_path = cif.replace('/cifs/', '/yamls/').replace('.cif', '.yaml')
            if os.path.exists(yaml_path):
                shutil.copy(yaml_path, os.path.join(new_rep, 'yamls', os.path.basename(yaml_path)))
                print(f"Copied {os.path.basename(yaml_path)} → {new_rep}/yamls/")
            else:
                print(f"[WARNING] No yaml found for {cif_name}")          

def remove_pkls(boltz_cache, run_cache_prot, VERBOSE):

    run_cache_csv = glob.glob(f'{run_cache_prot}/*.csv')

    run_cache_csv_df = pd.read_csv(run_cache_csv[0])

    for pkl_id in run_cache_csv_df['pkl_id']:
        pkl_file = os.path.join(boltz_cache, 'mols', f'{pkl_id}.pkl')
        if os.path.exists(pkl_file):
            os.remove(pkl_file)
    
    if VERBOSE: print(f"Deleted all generated pkls.")


def get_protein_stats(pred_df, protein=None):
    """
    Get statistics about predictions for each protein.
    Useful for monitoring progress across different proteins.
    
    Args:
        pred_df: predictions DataFrame (long format with 'protein' column)
        protein: specific protein to analyze (None = all proteins)
    
    Returns:
        DataFrame with stats per protein
    """
    if 'protein' not in pred_df.columns:
        print("Warning: 'protein' column not found. Cannot generate protein stats.")
        return None
    
    if protein:
        pred_df = pred_df[pred_df['protein'] == protein]
    
    stats = pred_df.groupby('protein').agg({
        'substance_id': 'count',  # Total rows
        'pred_log10ic50': lambda x: x.notna().sum(),  # Non-null predictions
        'boltz_runID': 'nunique'  # Unique runs
    }).rename(columns={
        'substance_id': 'total_rows',
        'pred_log10ic50': 'valid_predictions',
        'boltz_runID': 'unique_runs'
    })
    
    return stats