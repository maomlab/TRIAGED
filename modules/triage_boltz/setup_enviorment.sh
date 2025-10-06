#!/bin/bash

# Set the environment variable PROJECT_DIR
export PROJECT_DIR="/PATH/TO/PROJECT"  # Replace with your actual project path
export PYTHONPATH="/PATH/TO/MODULES/modules:$PYTHONPATH"
export SLURM_ACCOUNT="SLURM_ACCOUNT_NAME"  # Replace with your actual SLURM account name
# Set the environment variable for SLURM email notifications
export SLURM_EMAIL="SLURM_EMAIL_ADDRESS"  # Replace with your actual SLURM email address