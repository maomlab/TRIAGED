#!/bin/bash
#SBATCH --job-name=test_env_tgcpl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=6:00:00
#SBATCH --account=tromeara99
#SBATCH --partition=gpu
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-gpu=30GB
#SBATCH --gres=gpu:1
#SBATCH --array=0
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null


source /nfs/turbo/umms-maom/ymanasa/miniconda3/etc/profile.d/conda.sh
conda activate boltz2
module load cuda cudnn

JOB_FILE="$1"
LINES_PER_TASK=10   # number of YAMLs per array task

# Calculate starting line for this array task
START=$(( SLURM_ARRAY_TASK_ID * LINES_PER_TASK + 1 ))
END=$(( START + LINES_PER_TASK - 1 ))

# Loop through the lines for this task
sed -n "${START},${END}p" "$JOB_FILE" | while read -r YAML_PATH RECEPTOR_LIG_DIR; do
    echo "Running boltz predict for: $YAML_PATH -> $RECEPTOR_LIG_DIR"

    boltz predict "${YAML_PATH}" \
        --out_dir "${RECEPTOR_LIG_DIR}" \
        --num_workers 8 \
        --sampling_steps 500 \
        --sampling_steps_affinity 500 \
        1> "${RECEPTOR_LIG_DIR}/slurm_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out" \
        2> "${RECEPTOR_LIG_DIR}/slurm_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"
done