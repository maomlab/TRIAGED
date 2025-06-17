import os
import csv
import sys
import subprocess
from pdb_to_fasta import pdb_to_fasta
from mmseqs2 import run_mmseqs2
import argparse

from rdkit import Chem
from rdkit.Chem import AllChem


def check_smiles(smiles: str, verbose: bool = True):
    """
    Attempts to load and sanitize a SMILES string using RDKit.
    Returns a canonicalized SMILES string if successful, otherwise None.
    
    Parameters:
    - smiles (str): The input SMILES string.
    - verbose (bool): If True, print debug messages on failure.

    Returns:
    - str or None: A valid, canonical SMILES or None if the molecule is invalid.
    """
    try:
        # Attempt to parse without sanitizing
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            if verbose:
                print(f"[ERROR] MolFromSmiles failed for: {smiles}")
            return None
        # Attempt sanitization (includes valence check, aromaticity, Hs)
        Chem.SanitizeMol(mol)
        # Return the canonical SMILES
        return Chem.MolToSmiles(mol, canonical=True)
    except Exception as e:
        if verbose:
            print(f"[ERROR] Sanitization failed for SMILES: {smiles}\n{e}")
        return None
    
def ensure_environment_variables():
    """
    Ensures necessary environment variables are set. If not, runs setup_enviorment.sh.
    """
    if not os.getenv("PROJECT_DIR"):
        print("Environment variable PROJECT_DIR is not set. Running setup_enviorment.sh...")
        setup_script = os.path.join(os.path.dirname(__file__), "setup_enviorment.sh")
        subprocess.run(f"source {setup_script}", shell=True, executable="/bin/bash", check=True)
        print("Environment variables set successfully.")


def create_boltz_job(csv_file, pdb_file, output_dir, num_jobs, covalent_docking=False):
    """
    Creates directories with .yaml files based on the input CSV and PDB files.
    :param csv_file: Path to the input CSV file.
    :param pdb_file: Path to the input PDB file.
    :param output_dir: Path to the output directory.
    :param covalent_docking: Whether to write yaml for covalent docking. 
    """
    # Ensure environment variables are set
    ensure_environment_variables()

    # Create directories based on num_jobs
    if num_jobs > 1:
        job_dirs = [f"{output_dir}_{i}" for i in range(1, num_jobs + 1)]
        for job_dir in job_dirs:
            os.makedirs(job_dir, exist_ok=True)
    else:
        job_dirs = [output_dir]
        os.makedirs(output_dir, exist_ok=True)

    # Check if the PDB file exists
    if not os.path.isfile(pdb_file):
        print(f"Error: PDB file '{pdb_file}' does not exist.")
        sys.exit(1)
    # Check if the CSV file exists only for non-covalent docking 
    if not covalent_docking:
        if not os.path.isfile(csv_file):
            print(f"Error: CSV file '{csv_file}' does not exist.")
            sys.exit(1)
    
    # Convert PDB to FASTA
    project_dir = os.getenv("PROJECT_DIR", "/home/ymanasa/turbo/ymanasa/boltz_benchmark")
    fasta_dir = os.path.join(project_dir, "/input_files/fasta/")
    pdb_name = os.path.splitext(os.path.basename(pdb_file))[0]
    fasta_file = os.path.join(fasta_dir, f"{pdb_name}.fasta")
    pdb_to_fasta(pdb_file, fasta_file)
    print(f"FASTA file created at {fasta_file}")
    # Read the FASTA sequence
    with open(fasta_file, 'r') as fasta:
        sequence = fasta.read().strip()

    if not covalent_docking:

        # Generate MSA using mmseqs2
        if os.path.exists(os.path.join(project_dir, "input_files/msa", f"{pdb_name}_mmseqs2.a3m")):
            print(f"MSA file already exists at {os.path.join(project_dir, 'input_files/msa', f'{pdb_name}_mmseqs2.a3m')}. Skipping MSA generation.")

        print("now pregenerating MSA with mmseqs2...")
        msa_dir = os.path.join(project_dir, "input_files/msa/")
        os.makedirs(msa_dir, exist_ok=True)
        msa_file = os.path.join(msa_dir, f"{pdb_name}_mmseqs2.a3m")
        msa_result = run_mmseqs2(sequence, prefix=f"{msa_dir}/{pdb_name}")
        with open(msa_file, 'w') as msa_output:
            msa_output.write("\n".join(msa_result))
        print(f"MSA saved to {msa_file}")
        if not os.path.isfile(msa_file):
            print(f"Error: MSA file '{msa_file}' was not created successfully.")
            sys.exit(1)

        # Read CSV file
        with open(csv_file, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            rows = list(reader)
            # Distribute rows evenly across job directories
            chunk_size = len(rows) // num_jobs
            chunks = [rows[i * chunk_size:(i + 1) * chunk_size] for i in range(num_jobs)]
            if len(rows) % num_jobs != 0:
                chunks[-1].extend(rows[num_jobs * chunk_size:])

            for job_dir, chunk in zip(job_dirs, chunks):
                for row in chunk:
                    catalog_id = row["compound_ID"]
                    smiles = row["SMILES"]
                    smiles = check_smiles(smiles, verbose=True)
                    if smiles is None:
                        print(f"Invalid SMILES for compound ID {catalog_id}: {row['SMILES']}")
                        continue
                    # Create .yaml file in the job directory
                    yaml_file = os.path.join(job_dir, f"{catalog_id}.yaml")
                    with open(yaml_file, 'w') as yaml:
                        yaml.write("version: 1\n")
                        yaml.write("sequences:\n")
                        yaml.write("  - protein:\n")
                        yaml.write("      id: A\n")
                        yaml.write(f"      sequence: {sequence}\n")
                        yaml.write(f"      msa: {msa_file}\n")

                        yaml.write("  - ligand:\n")
                        yaml.write(f"      id: B\n")
                        yaml.write(f"      smiles: '{smiles}'\n")
                        yaml.write("properties:\n")
                        yaml.write("  - affinity:\n")
                        yaml.write("      binder: B\n")
    else:
       
        with open(pdb_file, 'r') as pdb:
            for line in pdb:
                if line.startswith("LINK"):
                    atom1_name = line[13:17].strip()
                    atom2_name = line[43:47].strip()

                    res_idx = line[23:27].strip()
                    res_name = line[17:21].strip()

                    ccd = line[47:50].strip() 
                    break 
        yaml_file = os.path.join(output_dir, f"{pdb_name}_{ccd}.yaml")
        with open(yaml_file, 'w') as yaml:
            yaml.write("version: 1\n")
            yaml.write("sequences:\n")
            yaml.write("  - protein:\n")
            yaml.write("      id: A\n")
            yaml.write(f"      sequence: {sequence}\n")
            # yaml.write(f"      msa: {msa_file}\n")
            yaml.write(f"      modification:\n")
            yaml.write(f"        - position: {res_idx}\n") 
            yaml.write(f"          ccd: {res_name}\n")

            yaml.write("  - ligand:\n")
            yaml.write(f"      id: LIG\n")
            yaml.write(f"      ccd: '{ccd}'\n")
            
            yaml.write("constraints:\n")
            yaml.write("    - bond:\n")
            yaml.write(f"        atom1: [A, {res_idx}, {atom1_name}]\n")
            yaml.write(f"        atom2: [LIG, 1, {atom2_name}]\n")

            yaml.write("properties:\n")
            yaml.write("    - affinity:\n")
            yaml.write(f"        binder: LIG\n")
    print(f"YAML files created successfully in directories: {', '.join(job_dirs)}.")


def create_slurm_submit_script(work_dir, receptor, output_dir, job_name="boltz_screen"):
    """
    Dynamically creates a SLURM submit script for Boltz jobs.
    :param work_dir: Path to the working directory.
    :param receptor: Name of the receptor.
    :param output_dir: Path to the output directory.
    :param job_name: Name of the SLURM job.
    """
    email = os.getenv("SLURM_EMAIL", "default_email@example.com")  # Use environment variable for email
    slurm_script_path = os.path.join(output_dir, f"slurm_{receptor}.sh")
    with open(slurm_script_path, "w") as slurm_script:
        slurm_script.write(f"""#!/bin/bash
#SBATCH --job-name={job_name}_{receptor}
#SBATCH --mail-type=FAIL,END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --account=maom99
#SBATCH --partition=spgpu
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-gpu=48000m
#SBATCH --mail-user={email}
#SBATCH --output={job_name}_{receptor}_slurm.log

WORK_DIR={work_dir}
RECEPTOR={receptor}
mkdir -p ../outputs/${{RECEPTOR}}
module load cuda cudnn
boltz predict ${{WORK_DIR}} --out_dir ../outputs/${{RECEPTOR}} --num_workers 8
""")
    print(f"SLURM submit script created at {slurm_script_path}")

def main():
    parser = argparse.ArgumentParser(description="Setup Boltz job directories and YAML files.")
    parser.add_argument("-i","--input_csv_file", type=str, required=False,help="Path to the input CSV file. Required for non-covalent docking.")
    parser.add_argument("-p","--input_pdb_file", type=str, required=True,help="Path to the input PDB file.")
    parser.add_argument("-o","--output_directory", type=str, required=True,help="Path to the output directory.")
    parser.add_argument("-n", "--num_jobs", type=int, required=False, default=1, help="Number of jobs to create. Default is 1.")
    parser.add_argument("-c", "--covalent_docking", action='store_true', default=False,help="Whether ligand must covlanetly interact with protein")

    args = parser.parse_args()
    if not args.covalent_docking and args.input_csv_file is None:
        parser.error("--input_csv_file is required when --covalent_docking is False for non-covalent docking")

    create_boltz_job(args.input_csv_file, args.input_pdb_file, args.output_directory, args.num_jobs, args.covalent_docking)
    # Create SLURM submit script
    if args.num_jobs > 1:
        pdb_name= os.path.splitext(os.path.basename(args.input_pdb_file))[0]
        for i in range(args.num_jobs):
            create_slurm_submit_script(
                work_dir=os.path.join(args.output_directory, f"{args.output_directory}_{i+1}"),
                receptor=pdb_name,
                output_dir=f"{pdb_name}_slurm_submit_{i+1}.sh",
                job_name=f"boltz_screen"
            )
    else:
        create_slurm_submit_script(
                work_dir=os.path.join(args.output_directory, f"{args.output_directory}"),
                receptor=pdb_name,
                output_dir=f"{pdb_name}_slurm_submit_{i+1}.sh",
                job_name=f"boltz_screen"
            )
if __name__ == "__main__":
    main()
