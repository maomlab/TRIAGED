## Input Preprocessing for Covalent Inference
Environment required: ccd_pkl (need to upload yaml of the environemnt after testing)
Following scripts are in `covalent_scripts/`
### `preprocessing/make_input_csv.py`

This script sets up a CSV file with necessary information regarding the ligand and protein for writing `.yaml` files using setup_cov_jobs.py for covalent docking with Boltz2. Also creates `.pkl` files necessary for docking in default location: `/home/$USER/.boltz/mols/`.

#### Command-Line Arguments 

| Argument               | Type    | Required | Description                                                                                                    |
|------------------------|---------|----------|---------------------------------------------------------------------------------------------------------------|
| `--name`          | `str`   | Yes      | Name of protein. Used for naming output files. Can be PDB ID, Uniprot ID, anything unique. |
| `--prot_file`     | `str`   | Yes      | Path to either a PDB file or a TXT file with a single chain sequence. Do not include any headers.                                                                                                       |
| `--res_idx`       | `int`   | Yes      | Index of the residue to be covalently targeted by a covalent ligand. Starting at 1. Please confirm index matches expected residue in sequence/PDB provided.                                                                                       |                                       
| `--lig_csv`       | `bool`  | No       | Path to CSV with Ligand info. Do not include any headers with rows `Compound ID,SMILES`. |
| `--out_csv`       | `str`   | Yes      | Path to output CSV. Will be formatted to work with setup_cov_job.py. |
| `--ccd_db`       | `str`    | No      | Path to output CSV. Will be formatted to work with setup_cov_job.py. Path to directory with covalent compound pkl files. Default: `/home/$USER/.boltz/mols`.                                                  |

#### Example Usage
For a single protein and several ligands:
```bash
python preprocessing/make_input_csv.py --prot_file test/5MAJ.pdb --res_idx 25 --lig_csv test/test_ligs.csv --out_csv test/test_out.csv
```
Find example input & output files in `preprocessing/test/*`

### `setup_cov_job.py`

This script generates `.yaml` files based on the input CSV, suitable for covalent docking with Boltz2.

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


### `submit_jobs.py`



#### Command-Line Arguments




#### Example Usage
For a single protein and several ligands:
```bash
python setup_cov_job.py --input_csv_file /path/to/input.csv --output_directory /path/to/output_directory
```



Scripts to populate your columns with appropriate information can be found in `preprocessing/`. Without the correct covalent information, Boltz will not run! 
#### Generating Input CSV
Run `make_input_csv.py`

#### Output
- `.yaml` files are created in the specified output directory
---
## SLURM job scripts 
### slurm_run_boltz.sh

The script will automatically generate SLURM batch scripts for running Boltz jobs on a high-performance computing cluster with the pathing and inputs for boltz already set 
to the generated input directories. By default these scripts will be generated in {Project_dir}/slurm_scripts

#### Example Usage

```bash
sbatch slurm_run_boltz.sh
```

### covalent_jobs.py

This script helps setup ...

### covalent_slurm.sh 

...
