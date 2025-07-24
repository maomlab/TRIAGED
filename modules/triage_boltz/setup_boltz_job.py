import os
import sys
import subprocess

import argparse
from covalent_utils import get_link_atoms
from typing import Union
from triage_biomolecule.triage_biomolecule import TriageBiomolecule
from misc_utils import *



def create_boltz_job(csv_file: str, output_dir: str, num_jobs: int, covalent_docking: bool=False):
    """
    Creates directories with .yaml files based on the input CSV and PDB files.
    :param csv_file: Path to the input CSV file. Not required for covalent docking. 
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

    #creating triage_biomolecule objects from csv files
    triage_list = TriageBiomolecule.from_csv(csv_file)
    
    #creating fasta if missing...
    for triage_row in triage_list:
        for triage_biomolecule in triage_row:
            if triage_biomolecule.entity_type == "protein":
                if triage_biomolecule.sequence is None:
                    if triage_biomolecule.pdb_path is not None:
                        triage_biomolecule.generate_fasta_from_pdb_path()
                    elif triage_biomolecule.pdb_id is not None:
                        triage_biomolecule.fetch_fasta_from_pdb_id()
    #now computing msa's
    for triage_row in triage_list:
        for triage_biomolecule in triage_row:
            if triage_biomolecule.entity_type == "protein":
                if triage_biomolecule.sequence is not None:
                    triage_biomolecule.generate_single_msa() #generate single MSA
                    # Check if there are multiple protein entity types in the triage_row
                    protein_entities = [entity for entity in triage_row if entity.entity_type == "protein"]
                    if len(protein_entities) > 1:
                        print(f"Generating paired MSA for multiple protein entities in triage_row...")
                        proteins = [protein for protein in protein_entities if protein is not None and protein != triage_biomolecule]
                        if proteins:
                            triage_biomolecule.generate_pair_msa(proteins)  # Generate paired MSA
                            triage_biomolecule.create_combined_msa()

    if not covalent_docking:
        chunk_size = len(triage_list) // num_jobs
        chunks = [triage_list[i * chunk_size:(i + 1) * chunk_size] for i in range(num_jobs)]
        if len(triage_list) % num_jobs != 0:
            chunks[-1].extend(triage_list[num_jobs * chunk_size:])
        for job_dir, chunk in zip(job_dirs, chunks):
            for row in chunk:
                
                #finding the first protein and first ligand in the row
                first_ligand = None
                first_protein = None
                for entry in row:
                    if first_ligand is None and entry.entity_type == "ligand":
                        first_ligand = entry
                    if first_protein is None and entry.entity_type == "protein":
                        first_protein = entry
                    if first_ligand and first_protein:
                        break
                yaml_file = os.path.join(job_dir, f"{first_ligand.entity_id}.yaml")
                
                with open(yaml_file, 'w') as yaml:
                    yaml.write(f"version: 1\n")
                    yaml.write("sequences:\n")
                    
                    total_compound_num = 0
                    for entry in row:   
                        if entry.entity_type == "ligand":
                            catalog_id = entry.entity_id
                            smiles = entry.smiles
                            num = entry.num  # Default to 1 if not provided
                            # Sanitize and validate SMILES
                            catalog_id = sanitize_compound_id(catalog_id)
                            smiles = check_smiles(smiles, verbose=True)
                            if smiles is None:
                                print(f"Invalid SMILES for compound ID {catalog_id}, skipping...")
                                continue
                            for compound_num in range(1, num + 1):
                                ligand_chain_id = chr(ord('A') + (total_compound_num))
                                yaml.write("  - ligand:\n")
                                yaml.write(f"      id: {ligand_chain_id}\n")
                                yaml.write(f"      smiles: '{smiles}'\n")
                                total_compound_num += 1
                        elif entry.entity_type == "protein":
                            protein_sequence = entry.sequence
                            protein_num = entry.num  # Default to 1 if not provided
                            protein_msa = entry.msa_path
                            protein_pair_msa = entry.pair_msa_path if hasattr(entry, 'pair_msa_path') else None
                            # Create .yaml entries for each protein 
                            for num in range(1, int(protein_num) + 1):
                                protein_chain_id = chr(ord('A') + (total_compound_num))
                                yaml.write("  - protein:\n")
                                yaml.write(f"      id: {protein_chain_id}\n")
                                yaml.write(f"      sequence: {protein_sequence}\n")
                                if protein_pair_msa:
                                    yaml.write(f"      msa: {entry.combined_msa_path}\n")
                                else:
                                    yaml.write(f"      msa: {protein_msa}\n")
                                total_compound_num += 1
                        elif entry.entity_type == "nucleic_acid":
                            nucleic_acid_sequence = entry.sequence
                            nucleic_acid_num = entry.num  # Default to 1 if not provided
                            # Create .yaml entries for each nucleic acid
                            for num in range(1, int(nucleic_acid_num) + 1):
                                nucleic_acid_chain_id = chr(ord('A') + (total_compound_num))
                                yaml.write("  - dna:\n")
                                yaml.write(f"      id: {nucleic_acid_chain_id}\n")
                                yaml.write(f"      sequence: {nucleic_acid_sequence}\n")
                                total_compound_num += 1
                        elif entry.entity_type == "rna":
                            rna_sequence = entry.sequence
                            rna_num = entry.num
                            # Create .yaml entries for each RNA
                            for num in range(1, int(rna_num) + 1):
                                rna_chain_id = chr(ord('A') + (total_compound_num))
                                yaml.write("  - rna:\n")
                                yaml.write(f"      id: {rna_chain_id}\n")
                                yaml.write(f"      sequence: {rna_sequence}\n")
                                total_compound_num += 1
                    
                    yaml.write("properties:\n")
                    yaml.write("  - affinity:\n")
                    yaml.write(f"      binder: A\n")

    else:
        chunk_size = len(triage_list) // num_jobs
        chunks = [triage_list[i * chunk_size:(i + 1) * chunk_size] for i in range(num_jobs)]
        if len(triage_list) % num_jobs != 0:
            chunks[-1].extend(triage_list[num_jobs * chunk_size:])
        for job_dir, chunk in zip(job_dirs, chunks):
            for row in chunk:
                #finding the first protein in the row...
                for entry in row:
                    if first_protein is None and entry.entity_type == "protein":
                        first_protein = entry
                        if first_protein.pdb_path is not None:
                            pdb_file = first_protein.pdb_path
                        else: 
                            KeyError(f"Protein {first_protein.entity_id} does not have a PDB path, required for covelant docking.")
                    if first_protein:
                        break
                pdb_name = os.path.splitext(os.path.basename(pdb_file))[0]
                prot_atom, res_name, res_idx, _, ccd, lig_atom, _ = get_link_atoms(pdb_file)

                yaml_file = os.path.join(output_dir, f"{pdb_name}_{ccd}.yaml")
                sequence = first_protein.sequence
                
                with open(yaml_file, 'w') as yaml:
                    yaml.write("version: 1\n")
                    yaml.write("sequences:\n")
                    yaml.write("  - protein:\n")
                    yaml.write("      id: A\n")
                    yaml.write(f"      sequence: {sequence}\n")
                        # yaml.write(f"      msa: {msa_file}\n")
                    yaml.write(f"      modification:\n")
                    yaml.write(f"        - position: {res_idx}\n") # fix res idx 
                    yaml.write(f"          ccd: {res_name}\n")

                    yaml.write("  - ligand:\n")
                    yaml.write(f"      id: LIG\n")
                    yaml.write(f"      ccd: '{ccd}'\n")
                        
                    yaml.write("constraints:\n")
                    yaml.write("    - bond:\n")
                    yaml.write(f"        atom1: [A, {res_idx}, {prot_atom}]\n") # fix res idx 
                    yaml.write(f"        atom2: [LIG, 1, {lig_atom}]\n") # lig atm idx needs to be 1 regardless of PDB idx 

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
boltz predict {work_dir} --out_dir ../outputs/{receptor} --num_workers 8
""")
    print(f"SLURM submit script created at {slurm_script_path}")

def main():
    parser = argparse.ArgumentParser(
    description="Setup Boltz job directories and YAML files.",
    epilog="Additional Information:\n"
            "- check_smiles(smiles: str, verbose: bool = True): Validates and canonicalizes a SMILES string.\n"
            "- ensure_environment_variables(): Ensures necessary environment variables are set.\n"
            "- sanitize_compound_id(compound_ID: str): Sanitizes compound IDs for file-system compatibility.\n"
            "- parse_input_csv(csv_file: str): Parses a CSV file with compound and protein data.\n"
            "- generate_msa(project_dir: str, pdb_name: str, sequence: Union[str, list[str]]): Generates a Multiple Sequence Alignment (MSA).\n"
            "- create_boltz_job(csv_file: str, output_dir: str, num_jobs: int, covalent_docking: bool=False, protein_nmers: int=1): Creates job directories and YAML files.\n"
            "- create_slurm_submit_script(work_dir, receptor, output_dir, job_name='boltz_screen'): Creates a SLURM submit script for Boltz jobs.",
    formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("-i", "--input_csv_file", type=str, required=False, help="Path to the input CSV file. Required for non-covalent docking.")
    parser.add_argument("-o", "--output_directory", type=str, required=True, help="Path to the output directory.")
    parser.add_argument("-n", "--num_jobs", type=int, required=False, default=1, help="Number of jobs to create. Default is 1.")
    parser.add_argument("-c", "--covalent_docking", action='store_true', default=False, help="Whether ligand must covalently interact with protein.")
    parser.add_argument("--name", type=str, required=True, help="Name of the receptor (used for SLURM script naming).")
    args = parser.parse_args()

    #if args.input_pdb_file is None and args.input_fasta_file is None:
    #    print("Error: Either PDB file or FASTA file must be provided.")
    #    sys.exit(1)
    
    if not args.covalent_docking and args.input_csv_file is None:
        parser.error("--input_csv_file is required")

    # Ensure environment variables are set
    ensure_environment_variables()
    project_dir = os.getenv("PROJECT_DIR")

    
    create_boltz_job(args.input_csv_file, args.output_directory, args.num_jobs, args.covalent_docking)
    # Create SLURM submit script
    
    #need to create slurm scripts, need name of receptor.
    pdb_name = args.name
    if args.num_jobs > 1:
        for i in range(args.num_jobs):
            create_slurm_submit_script(
                work_dir=f"{args.output_directory}_{i+1}",
                receptor=pdb_name,
                output_dir=f"{project_dir}/slurm_scripts/{pdb_name}_slurm_submit_{i+1}.sh",
                job_name=f"boltz_screen"
            )
    else:
        create_slurm_submit_script(
                work_dir=f"{args.output_directory}",
                receptor=pdb_name,
                output_dir=f"{project_dir}/slurm_scripts/{pdb_name}_slurm_submit.sh",
                job_name=f"boltz_screen"
            )
if __name__ == "__main__":
    main()
