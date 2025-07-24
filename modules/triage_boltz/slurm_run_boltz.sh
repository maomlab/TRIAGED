#!/bin/bash
#SBATCH --job-name=ampc_boltz_screen
#SBATCH --mail-type=FAIL,END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:30:00
#SBATCH --account=maom99
#SBATCH --partition=spgpu
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-gpu=48000m
#SBATCH --mail-user=limcaoco@umich.edu
#SBATCH --output=ampc_boltz_slurm.log


WORK_DIR=${1}
RECEPTOR=${2}
BOLTZ_PATH=/home/limcaoco/opt/boltz
mkdir -p ../outputs/${RECEPTOR}
#conda activate boltz_env 
module load cuda cudnn
boltz predict ${WORK_DIR} --out_dir ../outputs/${RECEPTOR} --num_workers 8  --use_msa_server
