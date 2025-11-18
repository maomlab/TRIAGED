#!/bin/bash
source  /nfs/turbo/umms-maom/ymanasa/miniconda3/etc/profile.d/conda.sh
conda activate boltz2
module load cuda cudnn

boltz predict /home/ymanasa/turbo/ymanasa/opt/maom_boltz/covalent_module/test/single/5MAJ_Y7FCI.yaml --out_dir /home/ymanasa/turbo/ymanasa/opt/maom_boltz/covalent_module/test/single --num_workers 8 --use_msa_server 1> /home/ymanasa/turbo/ymanasa/opt/maom_boltz/covalent_module/test/single/boltz.out 2> /home/ymanasa/turbo/ymanasa/opt/maom_boltz/covalent_module/test/single/boltz.err

FILES=$(compgen -G '/home/ymanasa/turbo/ymanasa/opt/maom_boltz/covalent_module/test/single/boltz_results_*/predictions/*/*_model_0.cif')

if [ -z "$FILES" ]; then
    echo "Job failed or output files not found"
else
    echo "Boltz prediction completed successfully."
fi
