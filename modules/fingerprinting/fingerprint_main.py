# Example receptor lookup (update with your actual PDB paths)
from fingerprint_utils import generate_prolif_fingerprints_from_outputs

generate_prolif_fingerprints_from_outputs(
    outputs_dir="/home/limcaoco/turbo/limcaoco/boltz_benchmark/outputs/test",
    receptor_pdb_path="/home/limcaoco/turbo/limcaoco/boltz_benchmark/input_files/DUDEZ_pdbs/SRC.pdb",
    output_parquet="/home/limcaoco/turbo/limcaoco/boltz_benchmark/outputs/SRC/prolif_fingerprints.parquet",
    verbose=True
)

