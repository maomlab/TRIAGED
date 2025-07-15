#!/bin/bash

# Check if DUDEZ_list is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <DUDEZ_list>"
    exit 1
fi

DUDEZ_LIST=$1

# Check if the file exists
if [ ! -f "$DUDEZ_LIST" ]; then
    echo "Error: File '$DUDEZ_LIST' not found!"
    exit 1
fi
# Iterate through each receptor in the DUDEZ_list
while IFS= read -r receptor; do
    receptor=$(echo "$receptor" | xargs) # Remove leading and trailing whitespaces
    if [ -n "$receptor" ]; then
        echo "Setting up Boltz job for receptor: $receptor"
        python setup_boltz_job.py --input_csv_file "../input_files/DUDEZ_benchmark/${receptor}_combined_ids.csv" --input_fasta_file "../input_files/DUDEZ_fastas/${receptor}.fasta" --output_directory "../boltz_inputs/${receptor}" -n 5
    fi
done < "$DUDEZ_LIST"

echo "Boltz job setup completed for all receptors in $DUDEZ_LIST."
