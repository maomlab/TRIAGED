import pandas as pd
import numpy as np

def compute_logauc_from_dataframe(df, score_col, label_col):
    """
    Compute logAUC for a DataFrame using the same logic as bootstrap_tldr.py.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing scores and binary labels.
    score_col : str
        Name of the column containing the scores (lower is better).
    label_col : str
        Name of the column containing the binary labels (1 = true, 0 = false).
    
    Returns
    -------
    float
        logAUC value (as a percentage, e.g., 25.0 for 25%).
    """
    # Sort by score (ascending: lower is better)
    df_sorted = df.sort_values(by=score_col, ascending=True).reset_index(drop=True)
    labels = df_sorted[label_col].values
    num_lig = np.sum(labels == 1)
    num_dec = np.sum(labels == 0)
    found_ligand = 0
    results = []
    for i in range(len(labels)):
        if i % max(1, int(len(labels)/10000)) == 0:
            results.append([i-found_ligand, found_ligand])
        if labels[i] == 1:
            found_ligand += 1
    results.append([len(labels) - found_ligand, found_ligand])
    results.append([num_dec, num_lig])
    points = []
    for x in results:
        fpr = x[0]*100.0/num_dec if num_dec > 0 else 0
        tpr = x[1]*100.0/num_lig if num_lig > 0 else 0
        points.append([fpr, tpr])
    i = 0
    while i < len(points) and points[i][0] < 0.1:
        i += 1
    if i > 0 and i < len(points):
        slope = (points[i][1] - points[i-1][1])/(points[i][0] - points[i-1][0])
        intercept = points[i][1] - slope * points[i][0]
        point_one =  [0.100001, (slope * 0.100001 + intercept)]
        points.insert(i, point_one)
    LOGAUC_MAX = 1.0
    LOGAUC_MIN = 0.001
    RANDOM_LOGAUC = (LOGAUC_MAX-LOGAUC_MIN)/np.log(10)/np.log10(LOGAUC_MAX/LOGAUC_MIN)
    npoints = []
    for x in points:
        if (x[0] >= LOGAUC_MIN*100) and (x[0] <= LOGAUC_MAX*100):
            npoints.append([x[0]/100, x[1]/100])
    area = 0.0
    for point2, point1 in zip(npoints[1:], npoints[:-1]):
        if point2[0] - point1[0] < 0.000001:
            continue
        dx = point2[0] - point1[0]
        dy = point2[1] - point1[1]
        intercept = point2[1] - (dy)/(dx) * point2[0]
        area += dy/np.log(10) + intercept*(np.log10(point2[0])-np.log10(point1[0]))
    logauc = area/np.log10(LOGAUC_MAX/LOGAUC_MIN) - RANDOM_LOGAUC
    return logauc * 100


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Compute logAUC from a DataFrame CSV file.")
    parser.add_argument('-f', '--file', type=str, default=None, help='CSV file with data (default: use toy example)')
    parser.add_argument('-s', '--score_col', type=str, default='score', help='Column name for scores (default: score)')
    parser.add_argument('-l', '--label_col', type=str, default='is_binder', help='Column name for labels (default: label)')
    args = parser.parse_args()

    if args.file is None:
        raise ValueError("No file provided!")
    else:
        df = pd.read_csv(args.file)

    logauc_func = compute_logauc_from_dataframe(df, args.score_col, args.label_col)
    print(f"logAUC: {logauc_func:.6f}")
