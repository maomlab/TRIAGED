import os
import pandas as pd
import argparse
import json
from scoring_utils import read_boltz_predictions, compute_vscreen_metrics, bootstrap_single_keep_ratio

def format_for_tldr(input_dataframe,output_directory):
    pass

def main():
    parser = argparse.ArgumentParser(description="Process and format TLDR data.")
    parser.add_argument("-i","--input_directory", type=str, required=True, help="Directory containing boltz predictions")
    parser.add_argument("-o","--output_directory", type=str, required=True, help="Directory to save the processed output files.")    
    parser.add_argument("-m","--compute_metrics", type=bool, default=False, help="Whether to compute metrics or not. Default is False.")
    parser.add_argument("-b","--bootstrap", type=bool, default=False, help="Whether to compute bootstrap metrics or not. Default is False.")
    parser.add_argument("-v","--in-vitro-data", type=str, required=False, default=None, help="Path to the in vitro data file, requires a column stating if compound is a binder or not. Default is None.")
    args = parser.parse_args()

    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)

    # Read and process the predictions
    df = read_boltz_predictions(args.input_directory)
    #print(df)
    # Save the processed DataFrame to a CSV file
    output_file = os.path.join(args.output_directory, "processed_boltz_data_full.csv")
    df.to_csv(output_file, index=False)
    print(f"Processed data saved to {output_file}")

    if args.compute_metrics:
        # Compute metrics if required
        print("Computing metrics...")
        if args.in_vitro_data is None:
            print("No in vitro data provided, skipping metric computation.")
        else:
            in_vitro_df = pd.read_csv(args.in_vitro_data)
            if 'is_binder' not in in_vitro_df.columns:
                raise ValueError("In vitro data must contain a column named 'is_binder' to indicate if the compound is a binder.")
            else:
                positive_df = in_vitro_df[in_vitro_df['is_binder'] == True]
                negative_df = in_vitro_df[in_vitro_df['is_binder'] == False]
                print(f"True Positive compounds: {len(positive_df)}, True Negative compounds: {len(negative_df)}")
                
                metrics_output = {}
                for col in ["Affinity Pred Value",
                                 "kcal/mol",
                                 "Ligand IPTM",
                                 "Complex pLDDT",
                                 "Complex iPLDDT",
                                 "Complex PDE",
                                 "Complex iPDE"]:
                        metrics = compute_vscreen_metrics(positive_df, negative_df, df, col)
                        print(f"Metrics for {col}: {metrics}")
                        metrics_output[col] = {"Metrics": metrics}  # Include metrics in the JSON structure
                        if args.bootstrap:
                            print("Computing bootstrap metrics...")
                            bootstrap_metrics = compute_vscreen_metrics(positive_df, negative_df, df, col, args.output_directory,bootstrap=True)
                            print(f"Bootstrap metrics for {col}: {bootstrap_metrics}")
                            metrics_output[col]["Bootstrap_Metrics"] = bootstrap_metrics  # Add bootstrap metrics to JSON structure
                
                # Export metrics to JSON
                metrics_json_file = os.path.join(args.output_directory, "computed_metrics.json")
                with open(metrics_json_file, 'w') as json_file:
                    json.dump(metrics_output, json_file, indent=4)
                print(f"Metrics saved to {metrics_json_file}")

if __name__ == "__main__":
    main()
