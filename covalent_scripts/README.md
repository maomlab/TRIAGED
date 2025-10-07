# Covalent Boltz2 Inference
Environment required: ccd_pkl and boltz2 (need to upload yaml of the environemnt after testing)
## Running Inference

The following are scripts to run inference with Boltz2 for covalent ligands on covalent proteins. 
The main script, `submit_jobs.py`, is written to support SLURM job management on GreatLakes clusters. 

### `covalent_scripts/submit_jobs.py`

Run predictions in command-line from the project directory (eg. `TRIAGED/`):
    `python -m covalent_scripts.submit_jobs [OPTIONS]`

#### Input Format
| Argument               | Type    | Required | Description                                                                                                    |
|------------------------|---------|----------|---------------------------------------------------------------------------------------------------------------|
| `--name`          | `str`   | Yes      | Name of protein. Used for naming output files. Can be PDB ID, Uniprot ID, anything unique. |
| `--prot_file`     | `str`   | Yes      | Path to either a PDB file or a TXT file with a single chain sequence. Do not include any headers.                                                                                                       |
| `--res_idx`       | `int`   | Yes      | Index of the residue to be covalently targeted by a covalent ligand. Starting at 1. Please confirm index matches expected residue in sequence/PDB provided.                                                       |         
| `--lig_chain`.    | `str`   | Yes      | Chain interacting with ligand in PDB file. Single character.        |                              
| `--lig_csv`       | `bool`  | Yes      | Path to CSV with Ligand info. Do not include any headers with rows `Compound ID,SMILES`. |
| `--outdir`       | `str`    | Yes      | Main output directory for all jobs. |
| `--ccd_db`       | `str`    | No      | Path to output CSV. Will be formatted to work with setup_cov_job.py. Path to directory with covalent compound pkl files. Default: `/home/$USER/.boltz/mols`.        |
| `--slurm`         | `str`    | No        | Path to SLURM template file. Default: `covalent_scripts/covalent_slurm.sh` |

## Input Preprocessing for Covalent Inference

### `preprocessing/make_input_csv.py`

This script sets up a CSV file with necessary information regarding the ligand and protein for writing `.yaml` files using `submit_jobs.py` for covalent docking with Boltz2. Also creates `.pkl` files necessary for docking in default location: `/home/$USER/.boltz/mols/`.

#### Command-Line Arguments 

| Argument               | Type    | Required | Description                                                                                                    |
|------------------------|---------|----------|---------------------------------------------------------------------------------------------------------------|
| `--name`          | `str`   | Yes      | Name of protein. Used for naming output files. Can be PDB ID, Uniprot ID, anything unique. |
| `--prot_file`     | `str`   | Yes      | Path to either a PDB file or a TXT file with a single chain sequence. Do not include any headers.                                                                                                       |
| `--res_idx`       | `int`   | Yes      | Index of the residue to be covalently targeted by a covalent ligand. Starting at 1. Please confirm index matches expected residue in sequence/PDB provided.                                                       |         
| `--lig_chain`.    | `str`   | Yes      | Chain interacting with ligand in PDB file. Single character.        |                              
| `--lig_csv`       | `bool`  | No       | Path to CSV with Ligand info. Do not include any headers with rows `Compound ID,SMILES`. |
| `--out_csv`       | `str`   | Yes      | Path to output CSV. Will be formatted to work with setup_cov_job.py. |
| `--ccd_db`       | `str`    | No      | Path to output CSV. Will be formatted to work with setup_cov_job.py. Path to directory with covalent compound pkl files. Default: `/home/$USER/.boltz/mols`.                                                  |

#### Example Usage
For a single protein and several ligands:
```bash
python preprocessing/make_input_csv.py --name 5MAJ --prot_file test/5MAJ.pdb --res_idx 25 --lig_chain A --lig_csv test/test_ligs.csv --out_csv test/test_out.csv --ccd_db /home/ymanasa/.boltz/mols
```
Find example input & output files in `preprocessing/test/*`


### `submit_jobs.py`

Makes .yamls needed for covalent docking with Boltz2 and submits an array of jobs with Slurm scripts. 

#### Command-Line Arguments




#### Example Usage
For a single protein and several ligands:
```bash
python setup_cov_job.py --input_csv_file /path/to/input.csv --output_directory /path/to/output_directory
```



Scripts to populate your columns with appropriate information can be found in `preprocessing/`. Without the correct covalent information, Boltz will not run! 

### `setup_cov_yamls.py`

This script generates `.yaml` files based on the input CSV, suitable for covalent docking with Boltz2. This does not need to be run seperately if `submit_jobs.py` is run directly with the appropriate inputs. The specific format required as input is listed below.

#### Command-Line Arguments 

| Argument               | Type    | Required | Description                                                                                                    |
|------------------------|---------|----------|---------------------------------------------------------------------------------------------------------------|
| `--input_csv_file`     | `str`   | Yes      | Path to the input CSV file. MUST have these columns: Compound_ID, SMILES, CCD, WH_Type, Lig_Atom, Prot_ID, Prot_Seq, Res_Idx, Res_Name, Res_Atom. Look at `Input CSV Requirements` for more info.                                                                                                          |
| `--output_directory`   | `int`   | Yes      | Path to the output directory for generated yamls.                                                                                                         |

#### Input CSV Requirements
For coavlent docking, the input CSV file must contain the following columns: 
- `Compound_ID`: Your unique identifier for the compound. 
- `SMILES`: SMILES string representing the compound structure. 
- `CCD`: Unique 5 alphanumeric for the pkl file of the ligand. Automatically generated when `make_input_csv.py` is run.
- `WH_Type`: Only these warhead types are supported currently (nitrile, alkylhalide, vinyl-sulfone, acrylamide, nitrile2) 
- `Lig_Atom`: Name of the ligand atom involved in the covalent bond. 
- `Prot_ID`: PDB or UNIPROT protein identifier code.
- `Res_Name`: Name of residue involved in the covalent bond.
- `Res_Idx`: Index of residue involved in the covalent bond, starting at 1.
- `Res_Atom`: Name of the atom in residue involved in the covalent bond. 
- `Prot_Seq`: Sequence of protein that matches indexing for res_idx.

#### Example Usage
```bash
python python setup_cov_job.py --input_csv_file preprocessing/test/test_out.csv --output_directory test_yamls/
```
Find example of input file in `preprocessing/test/*`
Fine example of output files in `test_yamls/`
