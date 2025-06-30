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
    if not os.getenv("PROJECT_DIR") or os.getenv(""):
        print("Environment variable PROJECT_DIR is not set. Running setup_enviorment.sh...")
        setup_script = os.path.join(os.path.dirname(__file__), "setup_enviorment.sh")
        subprocess.run(f"source {setup_script}", shell=True, executable="/bin/bash", check=True)
        print("Environment variables set successfully.")

def sanitize_compound_id(compound_ID):
    """
    Adjust the compound_ID to make it file-system friendly.
    Parameters:
    - compound_ID (str): The compound ID to sanitize.
    Returns:
    - str: The sanitized compound ID.
    """
    if any(char in compound_ID for char in [' ', '"', "'", ",", "(", ")", "/", "\\"]):
        sanitized_compound_id = (
            compound_ID.replace(" ", "_")  # Replace spaces with underscores
                       .replace('"', '')    # Remove double quotes
                       .replace("'", '')    # Remove single quotes
                       .replace(",", '')    # Remove commas
                       .replace("(", '')    # Remove parentheses
                       .replace(")", '')    # Remove parentheses
                       .replace("/", '_')   # Replace slashes with underscores
                       .replace("\\", '_')  # Replace backslashes with underscores
        )
        print(f"Found invalid characters in compound ID '{compound_ID}'. Sanitized to '{sanitized_compound_id}'.")
    else:
        sanitized_compound_id = compound_ID
    return sanitized_compound_id


def parse_input_csv(csv_file):
    """
    parse a CSV file with varying numbers of compound_ID, SMILES, and optional num columns.
    Also supports optional cofactor_ID_X, cofactor_sequence_X, and cofactor_num_X columns.
    If num_X or cofactor_num_X is missing or empty, it will be set to 1.

    Parameters:
    - csv_file: Path to the CSV file.

    Returns:
    - A list of dictionaries, where each dictionary contains compound_ID, SMILES, num, and optional cofactor data.
    """
    parsed_data = []

    with open(csv_file, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            compound_data = []
            for key in row.keys():
                if key.startswith("compound_ID") or key.startswith("peptide_ID") and row[key]:  
                    
                    if (len(key) > len("compound_ID")) or (len(key) > len("peptide_ID")):
                        if key.startswith("compound_ID"):
                            suffix = key[len("compound_ID"):]
                        else:
                            suffix = key[len("peptide_ID"):]  
                        smiles_key = f"SMILES{suffix}" 
                        
                        num_key = f"num{suffix}"  
                        peptide_id_key = f"peptide_ID{suffix}" 
                        peptide_sequence_key = f"peptide_sequence{suffix}"  
                        peptide_num_key = f"peptide_num{suffix}" 
                    else:
                        smiles_key = "SMILES"  
                        num_key = "num"  
                        peptide_id_key = "peptide_ID"  
                        peptide_sequence_key = "peptide_sequence"  
                        peptide_num_key = "peptide_num"  

                    
                    if key.startswith("compound_ID") and row[key]:
                        if smiles_key in row and row[smiles_key]:
                            compound_entry = {
                                "compound_ID": row[key],
                                "SMILES": row[smiles_key],
                                "num": row[num_key] if num_key in row and row[num_key] else "1"  
                            }
                        compound_data.append(compound_entry)
                    if key.startswith("peptide_ID") and row[key]:    
                        if peptide_id_key in row and row[peptide_id_key]:
                            compound_entry = {
                                "peptide_ID": row[peptide_id_key],
                                "peptide_sequence": row[peptide_sequence_key] if peptide_sequence_key in row else None,
                                "peptide_num": row[peptide_num_key] if peptide_num_key in row and row[peptide_num_key] else "1" 
                            }
                        compound_data.append(compound_entry)
            parsed_data.append(compound_data)

    return parsed_data


def generate_msa(project_dir, pdb_name, sequence):
    """
    Generate a Multiple Sequence Alignment (MSA) using mmseqs2.

    Parameters:
    - project_dir (str): The base project directory.
    - pdb_name (str): The name of the PDB file (without extension).
    - sequence (str): The sequence to use for MSA generation.

    Returns:
    - str: Path to the generated MSA file.
    """
    msa_dir = os.path.join(project_dir, "input_files/msa/")
    msa_file = os.path.join(msa_dir, f"{pdb_name}_mmseqs2.a3m")

    # Check if MSA file already exists
    if os.path.exists(msa_file):
        print(f"MSA file already exists at {msa_file}. Skipping MSA generation.")
        

    print("Now pregenerating MSA with mmseqs2...")
    os.makedirs(msa_dir, exist_ok=True)

    # Run mmseqs2 to generate MSA
    msa_result = run_mmseqs2(sequence, prefix=f"{msa_dir}/{pdb_name}")
    with open(msa_file, 'w') as msa_output:
        msa_output.write("\n".join(msa_result))
    print(f"MSA saved to {msa_file}")

    # Verify that the MSA file was created successfully
    if not os.path.isfile(msa_file):
        print(f"Error: MSA file '{msa_file}' was not created successfully.")
        sys.exit(1)


def create_boltz_job(csv_file, pdb_file, fasta_file, output_dir, num_jobs, covalent_docking=False, protein_nmers=1):
    """
    Creates directories with .yaml files based on the input CSV and PDB files.
    :param csv_file: Path to the input CSV file.
    :param pdb_file: Path to the input PDB file.
    :param output_dir: Path to the output directory.
    :param covalent_docking: Whether to write yaml for covalent docking. 
    """
    


    # Create directories based on num_jobs
    if num_jobs > 1:
        job_dirs = [f"{output_dir}_{i}" for i in range(1, num_jobs + 1)]
        for job_dir in job_dirs:
            os.makedirs(job_dir, exist_ok=True)
    else:
        job_dirs = [output_dir]
        os.makedirs(output_dir, exist_ok=True)


    # Check if the PDB file exists
    if pdb_file is not None:
        print(f"Input PDB file: {pdb_file}")
        print(f"generating FASTA file from PDB...")
        if not os.path.isfile(pdb_file):
            print(f"Error: PDB file '{pdb_file}' does not exist.")
            sys.exit(1)
        # Convert PDB to FASTA
        project_dir = os.getenv("PROJECT_DIR", "/home/ymanasa/turbo/ymanasa/boltz_benchmark")
        fasta_dir = os.path.join(project_dir, "/input_files/fasta/")
        pdb_name = os.path.splitext(os.path.basename(pdb_file))[0]
        fasta_file = os.path.join(fasta_dir, f"{pdb_name}.fasta")
        pdb_to_fasta(args.input_pdb_file, fasta_file)
        print(f"FASTA file created at {fasta_file}")
    else:
        project_dir = os.getenv("PROJECT_DIR", "/home/ymanasa/turbo/ymanasa/boltz_benchmark")
        pdb_name = os.path.splitext(os.path.basename(fasta_file))[0]
        print(f"Using provided FASTA file: {fasta_file}")

    # Check if the CSV file exists only for non-covalent docking 
    if not covalent_docking:
        if not os.path.isfile(csv_file):
            print(f"Error: CSV file '{csv_file}' does not exist.")
            sys.exit(1)
    
    
    # Read the FASTA sequence
    with open(fasta_file, 'r') as fasta:
        sequence = fasta.read().strip()

    if not covalent_docking:

        # Generate MSA using mmseqs2
        msa_dir = os.path.join(project_dir, "input_files/msa/")
        msa_file = os.path.join(msa_dir, f"{pdb_name}_mmseqs2.a3m")
        generate_msa(project_dir, pdb_name, sequence)
        

        # Read CSV file
        parsed_data = parse_input_csv(csv_file)
        
        #generating MSA for peptides if peptide_ID is present in the CSV
        for row in parsed_data:
            for entry in row:
                if "peptide_ID" in entry:
                    peptide_msa = os.path.join(msa_dir, f"{entry['peptide_ID']}_peptide_mmseqs2.a3m")
                    peptide_sequence = sequence
                    if not os.path.isfile(peptide_msa):
                        print(f"Generating MSA for peptides...")
                        pdb_name = entry["peptide_ID"]
                        peptide_msa = generate_msa(project_dir, pdb_name, peptide_sequence)          

        chunk_size = len(parsed_data) // num_jobs
        chunks = [parsed_data[i * chunk_size:(i + 1) * chunk_size] for i in range(num_jobs)]
        if len(parsed_data) % num_jobs != 0:
            chunks[-1].extend(parsed_data[num_jobs * chunk_size:])
        for job_dir, chunk in zip(job_dirs, chunks):
            for row in chunk:
                
                # Find the first compound_ID in the row to name the yaml_file
                first_compound_id = None
                for entry in row:
                    if "compound_ID" in entry:
                        first_compound_id = sanitize_compound_id(entry["compound_ID"])
                        break
                if first_compound_id:
                    yaml_file = os.path.join(job_dir, f"{first_compound_id}.yaml")
                
                with open(yaml_file, 'w') as yaml:
                    yaml.write(f"version: 1\n")
                    yaml.write("sequences:\n")
                    if protein_nmers > 1:
                        for protein_nmer in range(1, protein_nmers + 1):
                            chain_id = chr(ord('A') + protein_nmer - 1)
                            yaml.write("  - protein:\n")
                            yaml.write(f"      id: {chain_id}\n")
                            yaml.write(f"      sequence: {sequence}\n")
                            yaml.write(f"      msa: {msa_file}\n")
                    else:
                        yaml.write("  - protein:\n")
                        yaml.write(f"      id: A\n")
                        yaml.write(f"      sequence: {sequence}\n")
                        yaml.write(f"      msa: {msa_file}\n")
                    
                    total_compound_num = 0
                    for entry in row:   
                        if "compound_ID" in entry:
                            catalog_id = entry["compound_ID"]
                            smiles = entry["SMILES"]
                            num = entry.get("num", "1")  # Default to 1 if not provided
                            # Sanitize and validate SMILES
                            catalog_id = sanitize_compound_id(catalog_id)
                            smiles = check_smiles(smiles, verbose=True)
                            if smiles is None:
                                print(f"Invalid SMILES for compound ID {catalog_id}, skipping...")
                                continue
                            for compound_num in range(1, int(num) + 1):
                                ligand_chain_id = chr(ord('A') + (protein_nmers + total_compound_num))
                                yaml.write("  - ligand:\n")
                                yaml.write(f"      id: {ligand_chain_id}\n")
                                yaml.write(f"      smiles: '{smiles}'\n")
                                total_compound_num += 1
                        elif "peptide_ID" in entry:
                            peptide_id = entry["peptide_ID"]
                            peptide_sequence = entry["peptide_sequence"]
                            peptide_num = entry.get("peptide_num", "1")  # Default to 1 if not provided
                            # Create .yaml file for the peptide
                            for compound_num in range(1, int(peptide_num) + 1):
                                for protein_nmer in range(1, protein_nmers + 1):
                                    peptide_chain_id = chr(ord('A') + (protein_nmers + total_compound_num))
                                    yaml.write("  - protein:\n")
                                    yaml.write(f"      id: {peptide_chain_id}\n")
                                    yaml.write(f"      sequence: {peptide_sequence}\n")
                                    yaml.write(f"      msa: {peptide_msa}\n")
                                    total_compound_num += 1

                    yaml.write
                    yaml.write("properties:\n")
                    yaml.write("  - affinity:\n")
                    yaml.write(f"      binder: {chr(ord('A') + (protein_nmers + compound_num - 1))}\n")

    else: #do covelant method
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
    slurm_script_path = output_dir
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

mkdir -p ../outputs/{receptor}
module load cuda cudnn
boltz predict {work_dir} --out_dir ../outputs/{receptor} --num_workers 8 --method "md"
""")
    print(f"SLURM submit script created at {slurm_script_path}")

def main():
    parser = argparse.ArgumentParser(description="Setup Boltz job directories and YAML files.")
    parser.add_argument("-i","--input_csv_file", type=str, required=False,help="Path to the input CSV file. Required for non-covalent docking.")
    parser.add_argument("-p","--input_pdb_file", type=str, required=False,help="Path to the input PDB file.")
    parser.add_argument("-f","--input_fasta_file", type=str, required=False,help="Path to the input FASTA file. If provided, it will be used instead of generating from PDB.")
    parser.add_argument("-o","--output_directory", type=str, required=True,help="Path to the output directory.")
    parser.add_argument("-n", "--num_jobs", type=int, required=False, default=1, help="Number of jobs to create. Default is 1.")
    parser.add_argument("-c", "--covalent_docking", action='store_true', default=False,help="Whether ligand must covlanetly interact with protein")
    parser.add_argument("--protein_nmers", type=int, required=False, default=1, help="Number of protein nmers in the complex. Default is 1.")
    args = parser.parse_args()
   

    if args.input_pdb_file is None and args.input_fasta_file is None:
        print("Error: Either PDB file or FASTA file must be provided.")
        sys.exit(1)
    
    if not args.covalent_docking and args.input_csv_file is None:
        parser.error("--input_csv_file is required when --covalent_docking is False for non-covalent docking")

    # Ensure environment variables are set
    ensure_environment_variables()

    
    create_boltz_job(args.input_csv_file, args.input_pdb_file, args.input_fasta_file, args.output_directory, args.num_jobs, args.covalent_docking, protein_nmers=args.protein_nmers)
    # Create SLURM submit script
    if args.input_pdb_file is None:
        pdb_name = os.path.splitext(os.path.basename(args.input_fasta_file))[0]
    else:
        pdb_name= os.path.splitext(os.path.basename(args.input_pdb_file))[0]

    if args.num_jobs > 1:
        for i in range(args.num_jobs):
            create_slurm_submit_script(
                work_dir=f"{args.output_directory}_{i+1}",
                receptor=pdb_name,
                output_dir=f"../slurm_scripts/{pdb_name}_slurm_submit_{i+1}.sh",
                job_name=f"boltz_screen"
            )
    else:
        create_slurm_submit_script(
                work_dir=f"{args.output_directory}",
                receptor=pdb_name,
                output_dir=f"../slurm_scripts/{pdb_name}_slurm_submit.sh",
                job_name=f"boltz_screen"
            )
if __name__ == "__main__":
    main()
