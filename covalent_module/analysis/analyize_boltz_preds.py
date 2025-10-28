import os
import pandas as pd
import argparse
import json
from analysis_utils import read_boltz_predictions, process_invitro, compute_vscreen_metrics

# this is just for end user usage/ update later 

def main():
    parser = argparse.ArgumentParser(description="Process ground truth and Boltz prediction data. Outputs metrics and plots.")
    parser.add_argument("-i","--input_directory", type=str, required=True, help="Directory containing boltz predictions")
    parser.add_argument("-o","--output_directory", type=str, required=True, help="Directory to save the processed output files.")    
    # parser.add_argument("-m", "--compute_metrics", action="store_true", help="Flag to compute metrics.")
    # parser.add_argument("-b", "--bootstrap", action="store_true", help="Flag to compute bootstrap metrics.")
    parser.add_argument("-v","--in-vitro-data", type=str, required=False, default=None, help="Path to the in vitro data file, requires a column stating if compound is a binder or not. Default is None.")
    args = parser.parse_args()

    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)

    pred_df = read_boltz_predictions(args.input_directory)
    