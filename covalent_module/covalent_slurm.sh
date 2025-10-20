#!/bin/bash
#SBATCH --job-name=hscpl_200_samples_all_ligs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:10:00
#SBATCH --account=tromeara99
#SBATCH --partition=gpu
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-gpu=20GB
#SBATCH --gres=gpu:1
#SBATCH --array=0-277%10
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

source  /nfs/turbo/umms-maom/ymanasa/miniconda3/etc/profile.d/conda.sh

JOB_FILE="$1" 

YAML_PATH=$(awk -v id="$SLURM_ARRAY_TASK_ID" 'NR==id+1 {print $1}' "$JOB_FILE")
RECEPTOR_LIG_DIR=$(awk -v id="$SLURM_ARRAY_TASK_ID" 'NR==id+1 {print $2}' "$JOB_FILE")

conda activate boltz2
module load cuda cudnn

boltz predict "${YAML_PATH}" --out_dir "${RECEPTOR_LIG_DIR}" --num_workers 8 --use_msa_server 1> "${RECEPTOR_LIG_DIR}/slurm_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out" 2> "${RECEPTOR_LIG_DIR}/slurm_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"

# checking if run successfuly, otherwise logging err
FILES=$(compgen -G ${RECEPTOR_LIG_DIR}/boltz_results_*/predictions/*/*_model_0.cif)
if [ -z "$FILES" ]; then
    echo "Job failed or output files not found"
    FAILED_JOBS_FILE="$(dirname "$JOB_FILE")/failed_${SLURM_ARRAY_JOB_ID}_jobs.txt"
    echo "$RECEPTOR_LIG_DIR" >> "$FAILED_JOBS_FILE"
    exit 1
else
    echo "Job completed successfully"
fi

