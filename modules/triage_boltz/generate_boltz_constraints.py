import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings
import numpy as np
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
import xgboost as xgb
import shap

#currently untested! -miguel 

def get_residue_coords(pdb_path):
    """
    input: pdb_path (str): Path to the PDB file.
    output: pd.DataFrame: DataFrame with columns ['chain', 'residue', 'x', 'y', 'z']
    note: xyz coordinates are computed using the center of mass of each residue.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", pdb_path)
    records = []
    for model in structure:
        for chain in model:
            for residue in chain:
                # Only process standard residues (ignore water, ligands)
                hetfield, resseq, icode = residue.get_id()
                if hetfield.strip():  # skip hetero-residues
                    continue
                # Use built-in COM (mass-weighted)
                com = residue.center_of_mass(geometric=False)
                records.append({
                    "chain": chain.id,
                    "residue": resseq,
                    "x": com[0],
                    "y": com[1],
                    "z": com[2],
                })
    return pd.DataFrame(records, columns=["chain", "residue", "x", "y", "z"])

def generate_distogram(pdb_dataframe):
    """
    Generates a distogram (distance matrix) based on XYZ coordinates from the input DataFrame.
    :param pdb_dataframe: DataFrame with columns ['chain', 'residue', 'x', 'y', 'z'].
    :return: Numpy array representing the distogram.
    """
    # Extract XYZ coordinates
    coords = pdb_dataframe[['x', 'y', 'z']].values

    # Compute pairwise distances using scipy
    distogram = squareform(pdist(coords))

    return distogram

def flatten_upper_triangle(dm):
    iu = np.triu_indices_from(dm, k=1)
    return dm[iu], iu


def analyze_distogram_shap(X, y, pairs, n_estimators=100, max_depth=4, learning_rate=0.1):
    """
    Analyzes residue-pair importance using SHAP values from an XGBClassifier.
    
    Parameters:
    - X: numpy array or DataFrame of shape (n_samples, n_features), residue-residue distances.
    - y: 1D array of class labels (e.g., 0/1 for two protein states).
    - pairs: List of tuples (i, j) mapping feature indices to residue-pair identifiers.
    - n_estimators: Number of trees in the XGBClassifier (default: 100).
    - max_depth: Maximum depth of trees in the XGBClassifier (default: 4).
    - learning_rate: Learning rate for the XGBClassifier (default: 0.1).
    
    Returns:
    - pandas.DataFrame with columns:
      - Residue1: residue index i
      - Residue2: residue index j
      - shap_score: mean(|SHAP|) for the residue pair
    """
    print("training XGBClassifier to compute SHAP values...")
    # Train XGBClassifier
    model = xgb.XGBClassifier(
        n_estimators=n_estimators,
        max_depth=max_depth,
        learning_rate=learning_rate,
        use_label_encoder=False,
        eval_metric="logloss"
    )
    model.fit(X, y)
    print("XGBClassifier trained successfully.")
    # Compute SHAP values
    print("Calculating SHAP values...")
    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(X)
    print("SHAP values calculated successfully.")
    # Aggregate mean absolute SHAP values per feature
    mean_abs_shap = np.mean(np.abs(shap_values), axis=0)

    # Create DataFrame with residue-pair identifiers and SHAP scores
    shap_df = pd.DataFrame({
        "Residue1": [pair[0] for pair in pairs],
        "Residue2": [pair[1] for pair in pairs],
        "shap_score": mean_abs_shap
    })

    return shap_df

def generate_boltz_constraints(state_one_pdb_path, state_two_pdb_path, output_path, n=25, export_csv=True, use_SHAP=False):
    """
    Identifies top residue pairs for binding constraints using PCA on distograms of two protein states.
    input:
        state_one_pdb_path (str): Path to the first PDB file (state 1).
        state_two_pdb_path (str): Path to the second PDB file (state 2).
        output_path (str): Path to save the constraints CSV file.
        n (int): Number of top residue pairs to return.
        export_csv (bool): Whether to save the results to a CSV file.
    output:
        pd.DataFrame: DataFrame containing top residue pairs and PCA scores.
    """
    # Suppress warnings from Bio.PDB
    warnings.simplefilter('ignore', PDBConstructionWarning)
    
    # Get residue coordinates for both states
    df1 = get_residue_coords(state_one_pdb_path)
    df2 = get_residue_coords(state_two_pdb_path)
    
    # Generate distograms for both states
    distogram1 = generate_distogram(df1)
    distogram2 = generate_distogram(df2)
    
    # Compute the difference between the two distograms
    distogram_diff = np.abs(distogram1 - distogram2)
    
    # Flatten the upper triangle of the distogram difference
    flattened_diff, indices = flatten_upper_triangle(distogram_diff)
    
    # Perform PCA analysis
    pca = PCA(n_components=1)
    pca_scores = pca.fit_transform(flattened_diff.reshape(-1, 1)).flatten()
    
    # Get top residue pairs by PCA1
    top_indices = np.argsort(-pca_scores)[:n]  # Select top `n` pairs
    top_pairs = [(indices[0][i], indices[1][i], pca_scores[i]) for i in top_indices]
    top_pairs_df = pd.DataFrame(top_pairs, columns=["Residue 1", "Residue 2", "PCA Score"])
    
    if use_SHAP: 
         print("Calculating SHAP values for top residue pairs...")
         shap_top_pairs_df = analyze_distogram_shap(
            flattened_diff.reshape(1, -1),  # Reshape for single sample
            np.array([1]),  # Dummy label for binary classification
            [(indices[0][i], indices[1][i]) for i in top_indices]
        )
         top_pairs_df = top_pairs_df.merge(shap_top_pairs_df, on=["Residue 1", "Residue 2"], how="left")
    
    # Save top pairs to CSV if required
    if export_csv:
        top_pairs_df.to_csv(output_path, index=False)
        print(f"Top residue pairs saved to {output_path}")
    
    return top_pairs_df