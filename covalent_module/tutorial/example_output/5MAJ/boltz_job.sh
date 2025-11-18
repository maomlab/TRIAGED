#!/bin/bash
#!/bin/bash
#SBATCH --job-name=test_single
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:30:00
#SBATCH --account=tromeara99
#SBATCH --partition=gpu
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-gpu=10GB
#SBATCH --gres=gpu:1
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

source  /nfs/turbo/umms-maom/ymanasa/miniconda3/etc/profile.d/conda.sh
conda activate boltz2
module load cuda cudnn

boltz predict /home/ymanasa/turbo/ymanasa/opt/maom_boltz/covalent_module/tutorial/example_output/5MAJ/5MAJ_VM834.yaml --out_dir /home/ymanasa/turbo/ymanasa/opt/maom_boltz/covalent_module/tutorial/example_output/5MAJ --num_workers 8 --use_msa_server 1> /home/ymanasa/turbo/ymanasa/opt/maom_boltz/covalent_module/tutorial/example_output/5MAJ/boltz.out 2> /home/ymanasa/turbo/ymanasa/opt/maom_boltz/covalent_module/tutorial/example_output/5MAJ/boltz.err

FILES=$(compgen -G '/home/ymanasa/turbo/ymanasa/opt/maom_boltz/covalent_module/tutorial/example_output/5MAJ/boltz_results_*/predictions/*/*_model_0.cif')

if [ -z "$FILES" ]; then
    echo "Job failed or output files not found"
else
    echo "Boltz prediction completed successfully."
fi
