import os
import pandas as pd
import argparse
from scoring_utils import read_boltz_predictions

def format_for_tldr(input_dataframe,output_directory):
    pass

def main():
    parser = argparse.ArgumentParser(description="Process and format TLDR data.")
    parser.add_argument("-i","--input_directory", type=str, required=True, help="Directory containing boltz predictions")
    parser.add_argument("-o","--output_directory", type=str, required=True, help="Directory to save the processed output files.")    
    parser.add_argument("-m","--compute_metrics", type=bool, default=False, help="Whether to compute metrics or not. Default is False.")
    parser.add_argument("-p","--positives", type=str, required=False, default=None, help="Path to the confirmed postives file. Default is None.")
    parser.add_argument("-n","--negatives", type=str, required=False, default=None, help="Path to the confirmed negatives file. Default is None.")
    args = parser.parse_args()

    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)

    # Read and process the predictions
    df = read_boltz_predictions(args.input_directory)
    #print(df)
    # Save the processed DataFrame to a CSV file
    output_file = os.path.join(args.output_directory, "processed_tldr_data_full.csv")
    df.to_csv(output_file, index=False)
    print(f"Processed data saved to {output_file}")

    if args.compute_metrics:
        # Compute metrics if required
        print("Computing metrics...")
        
if __name__ == "__main__":
    main()
