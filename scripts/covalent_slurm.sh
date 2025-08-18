#!/bin/bash
#SBATCH --job-name=covalent_boltz
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:00:00
#SBATCH --account=tromeara99
#SBATCH --partition=spgpu
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-gpu=40000m
#SBATCH --gres=gpu:1
#SBATCH --array=0-4%2
#SBATCH --output=../../boltz/covalent_benchmark/predictions/logs/slurm_%A_%a.out
#SBATCH --error=../../boltz/covalent_benchmark/predictions/logs/slurm_%A_%a.err

# source  /nfs/turbo/umms-maom/ymanasa/miniconda3/etc/profile.d/conda.sh

JOB_FILE="job_input_list.txt"

YAML_PATH=$(awk -v id="$SLURM_ARRAY_TASK_ID" 'NR==id+1 {print $1}' $JOB_FILE)
RECEPTOR_PATH=$(awk -v id="$SLURM_ARRAY_TASK_ID" 'NR==id+1 {print $2}' $JOB_FILE)

conda activate boltz2
module load cuda cudnn

echo $RECEPTOR_PATH

boltz predict ${YAML_PATH} --out_dir ../../boltz/covalent_benchmark/predictions/${RECEPTOR_PATH} --num_workers 8  --use_msa_server 