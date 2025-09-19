#!/bin/bash
#SBATCH --job-name=covalent_boltz
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:30:00
#SBATCH --account=tromeara99
#SBATCH --partition=spgpu
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-gpu=40000m
#SBATCH --gres=gpu:1
#SBATCH --array=0-9%3

source  /nfs/turbo/umms-maom/ymanasa/miniconda3/etc/profile.d/conda.sh

JOB_FILE="$1" 
OUTDIR="$2"

YAML_PATH=$(awk -v id="$SLURM_ARRAY_TASK_ID" 'NR==id+1 {print $1}' "$JOB_FILE")
RECEPTOR_PATH=$(awk -v id="$SLURM_ARRAY_TASK_ID" 'NR==id+1 {print $2}' "$JOB_FILE")

FINAL_OUTDIR="${OUTDIR}/${RECEPTOR_PATH}"
mkdir -p "$FINAL_OUTDIR" # in case not made/erred 

conda activate boltz2
module load cuda cudnn

boltz predict "${YAML_PATH}" --out_dir "${FINAL_OUTDIR}" --num_workers 8 --use_msa_server 1> "${FINAL_OUTDIR}/slurm_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out" 2> "${FINAL_OUTDIR}/slurm_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"

# checking if run successfuly, otherwise logging err
EXIT_CODE=$?

if [ $EXIT_CODE -ne 0 ] || [ -z "$(compgen -G "${FINAL_OUTDIR}/boltz_results_*/predictions/*/*model_0.cif")" ]; then
    echo "Job failed or output files not found for receptor: $RECEPTOR_PATH"
    
    # Log to failed jobs file
    FAILED_JOBS_FILE="${OUTDIR}/failed_jobs.txt"
    echo "$RECEPTOR_PATH" >> "$FAILED_JOBS_FILE"
    
    exit 1
else
    echo "Job completed successfully for receptor: $RECEPTOR_PATH"
    echo "Found output files: $(compgen -G "${FINAL_OUTDIR}/boltz_results_*/predictions/*/*model_0.cif" | wc -l)"
fi