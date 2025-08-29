import os
import pandas as pd
import argparse
import json
from analysis_utils import read_boltz_predictions, compute_vscreen_metrics
import shutil
import numbers


def main():
    parser = argparse.ArgumentParser(description="Process and format TLDR data.")
    parser.add_argument("-i","--input-directory", type=str, required=True, help="Directory containing boltz predictions")
    parser.add_argument("-o","--output-directory", type=str, required=True, help="Directory to save the processed output files.")    
    parser.add_argument("-m", "--compute-metrics", action="store_true", help="Flag to compute metrics.")
    parser.add_argument("-b", "--bootstrap", action="store_true", help="Flag to compute bootstrap metrics.")
    parser.add_argument("-r","--reference-data", type=str, required=False, default=None, help="Path to the reference file, requires a column stating if compound is a binder or not. Default is None.")
    args = parser.parse_args()

    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)

    # Read and process the predictions
    # Helper function to check if a directory contains only prediction files
    def is_predictions_dir(path):
        return os.path.basename(os.path.normpath(path)) == "predictions"

    input_dir = args.input_directory
    parent_dir = os.path.dirname(input_dir.rstrip("/"))

    if is_predictions_dir(input_dir):
        print("predictions directory detected.")
        df = read_boltz_predictions(input_dir)
    else:
        print("combining prediction directory")
        # Find all subdirectories matching boltz_results_{receptor_name}_{num}
        subdirs = [
            os.path.join(input_dir, d)
            for d in os.listdir(input_dir)
            if os.path.isdir(os.path.join(input_dir, d)) and d.startswith("boltz_results_")
        ]
        print(f"Found {len(subdirs)} boltz_results subdirectories.")

        combined_predictions_dir = os.path.join(input_dir, "combined_predictions")

        if not os.path.exists(combined_predictions_dir):
            os.makedirs(combined_predictions_dir)
        
        for subdir in subdirs:
            pred_dir = os.path.join(subdir, "predictions")
            if os.path.exists(pred_dir):
                for compound_dir in os.listdir(pred_dir):
                    compound_path = os.path.join(pred_dir, compound_dir)
                    if os.path.isdir(compound_path):
                        dst = os.path.join(combined_predictions_dir, f"{compound_dir}")
                        shutil.move(compound_path, dst)
        df = read_boltz_predictions(combined_predictions_dir)
    #print(df)
    # Save the processed DataFrame to a CSV file
    output_file = os.path.join(args.output_directory, "processed_boltz_data_full.csv")
    df.to_csv(output_file, index=False)
    print(f"Processed data saved to {output_file}")

    if args.compute_metrics:
        # Compute metrics if required
        print("Computing metrics...")
        if args.reference_data is None:
            print("No reference data provided, skipping metric computation.")
        else:
            reference_df = pd.read_csv(args.reference_data)
            print(reference_df)
            if 'is_binder' not in reference_df.columns:
                raise ValueError("Reference data must contain a column named 'is_binder' to indicate if the compound is a binder.")
            else:
                def _to_bool(v):
                    if pd.isna(v):
                        return None
                    if isinstance(v, numbers.Number):
                        if v == 1:
                            return True
                        if v == 0:
                            return False
                    s = str(v).strip().lower()
                    if s in ('true', 't', '1', 'yes', 'y', 'TRUE', 'True'):
                        return True
                    if s in ('false', 'f', '0', 'no', 'n', 'FALSE', 'False'):
                        return False
                    return None

                bool_series = reference_df['is_binder'].apply(_to_bool)
                if bool_series.isnull().any():
                    print("[WARN] Some 'is_binder' values could not be interpreted and will be ignored.")
                positive_df = reference_df[bool_series == True].copy()
                negative_df = reference_df[bool_series == False].copy()
                print(f"True Positive compounds: {len(positive_df)}, True Negative compounds: {len(negative_df)}")
                
                metrics_output = {}
                #computing for "DOCK score" or "score"
                print(reference_df.columns)
                if 'DOCK score' in reference_df.columns:
                    try:
                        metrics = compute_vscreen_metrics(positive_df, negative_df, reference_df, 'DOCK score', args.output_directory)
                        print(f"Metrics for DOCK score: {metrics}")
                        metrics_output['DOCK score'] = {"Metrics": metrics}
                        if args.bootstrap:
                            bootstrap_metrics = compute_vscreen_metrics(positive_df, negative_df, reference_df, 'DOCK score', args.output_directory, bootstrap=True)
                            print(f"Bootstrap metrics for DOCK score: {bootstrap_metrics}")
                            metrics_output['DOCK score']["Bootstrap_Metrics"] = bootstrap_metrics
                    except Exception as e:
                        print(f"[WARN] Failed to compute DOCK score metrics: {e}")
                elif 'score' in reference_df.columns:    
                    try:
                        metrics = compute_vscreen_metrics(positive_df, negative_df, reference_df, 'score', args.output_directory)
                        print(f"Metrics for DOCK score: {metrics}")
                        metrics_output['DOCK score'] = {"Metrics": metrics}
                        if args.bootstrap:
                            bootstrap_metrics = compute_vscreen_metrics(positive_df, negative_df, reference_df, 'score', args.output_directory, bootstrap=True)
                            print(f"Bootstrap metrics for DOCK score: {bootstrap_metrics}")
                            metrics_output['DOCK score']["Bootstrap_Metrics"] = bootstrap_metrics
                    except Exception as e:
                        print(f"[WARN] Failed to compute score metrics: {e}")
                else:
                    print("No 'DOCK score' or 'score' column found in reference data, skipping metric computation for these columns.")


                for col in ["Affinity Pred Value",
                                 "Affinity Probability Binary",
                                 "kcal/mol",
                                 "Ligand IPTM",
                                 "Complex pLDDT",
                                 "Complex iPLDDT",
                                 "Complex PDE",
                                 "Complex iPDE"]:
                    # Ensure the predictions dataframe contains the score column
                    if col not in df.columns:
                        print(f"No '{col}' column found in predictions, skipping metrics for this column.")
                        continue
                    try:
                        
                        metrics = compute_vscreen_metrics(positive_df, negative_df, df, col, args.output_directory)
                        print(f"Metrics for {col}: {metrics}")
                        metrics_output[col] = {"Metrics": metrics}  # Include metrics in the JSON structure
                        if args.bootstrap:
                            print("Computing bootstrap metrics...")
                            bootstrap_metrics = compute_vscreen_metrics(positive_df, negative_df, df, col, args.output_directory, bootstrap=True)
                            print(f"Bootstrap metrics for {col}: {bootstrap_metrics}")
                            metrics_output[col]["Bootstrap_Metrics"] = bootstrap_metrics  # Add bootstrap metrics to JSON structure
                    except ValueError as e:
                        print(f"[WARN] Skipping metrics for {col} due to insufficient data: {e}")
                    except Exception as e:
                        print(f"[WARN] Failed to compute metrics for {col}: {e}")
                
                # Export metrics to JSON
                metrics_json_file = os.path.join(args.output_directory, "computed_metrics.json")
                with open(metrics_json_file, 'w') as json_file:
                    json.dump(metrics_output, json_file, indent=4)
                print(f"Metrics saved to {metrics_json_file}")

if __name__ == "__main__":
    main()
