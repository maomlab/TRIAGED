from fingerprint_utils import generate_prolif_fingerprints_from_outputs

# Minimal smoke test: point to one receptor folder if available
outputs_dir = '/home/limcaoco/turbo/limcaoco/boltz_benchmark/outputs'
receptor_pdb_dir = '/home/limcaoco/turbo/limcaoco/boltz_benchmark/input_files/DUDEZ_pdbs'
output_parquet = '/home/limcaoco/turbo/limcaoco/boltz_benchmark/outputs/test_prolif.parquet'

# Run only on a small subset: the function currently iterates all, so rely on existing small dataset
generate_prolif_fingerprints_from_outputs(outputs_dir, receptor_pdb_dir, output_parquet, pose_ext='.cif', verbose=True)
