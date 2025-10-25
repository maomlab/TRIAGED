# Boltz2 Inference for Covalent Docking

Environment required: ccd_pkl and boltz2 (need to upload yaml of the environemnt after testing)
Environment requied to run `single_inference.py` and `submit_jobs.py` is ccd_pkl 
Set environmental variables in setup_environment.sh 

## Single Ligand Inference 

The following script can be used to perform docking of a single covalent ligand to a protein with Boltz2.
This script requires a GPU to run. Can be run with and without SLURM. If using SLURM, please set SLURM arguments appropriately in slurm template (example: slurm_args.sh).

### `covalent_module/single_inference.py`
Run predictions in command-line from the project directory (eg. `TRIAGED/`):
    `python -m covalent_module.single_inference [OPTIONS]`

#### Input Format
| Argument   | Type | Required | Description               |
|-------------|------|-----------|---------------------------|
| `--prot_file`     | `str`   | Yes      | Path to a PDB file.                  |
| `--res_idx`       | `int`   | Yes      | Index of the residue to be covalently targeted by a covalent ligand. Starting at 1. Please confirm index matches expected residue in sequence/PDB provided. |         
| `--lig_chain`    | `str`   | No      | Chain interacting with ligand in PDB file. Default is 'A'. |                          
| `--msa_path` | `str` | No | Path to MSA file in csv format. If not provided, msa server will be used. |
| `--smiles` | `str` | Yes | SMILES of ligand. Will be processed to remove leaving group before docking. |
| `--id` | `str` | Yes | Compound ID for the ligand. | 
| `--outdir`       | `str`    | Yes      | Main output directory for Boltz. |
| `--slurm` | `str` | No | Path to SLURM template file. Default: None. |

#### Example Usage
``` bash
python -m covalent_module.single_inference -p covalent_module/test/single/5MAJ.pdb -r 25 -c A -m /home/ymanasa/turbo/ymanasa/opt/maom_boltz/covalent_module/test/single/5MAJ_msa.csv -i Y7FCI -s "O=C(N[C@@H](CC1=CC(C)=CC=C1)C(NCC#N)=O)C2=CC(C(C)(C)C)=NN2C" -o /home/ymanasa/turbo/ymanasa/opt/maom_boltz/covalent_module/test/single
```

## Running Virtual Screening

The following are scripts to run inference with Boltz2 for docking covalent ligands into covalent proteins. 
The main script, `submit_jobs.py`, is written to support SLURM job management on GreatLakes clusters in a virtual screening scenario. 

### `covalent_module/submit_jobs.py`
Run predictions in command-line from the project directory (eg. `TRIAGED/`):
    `python -m covalent_module.submit_jobs [OPTIONS]`

#### Input Format
| Argument   | Type | Required | Description               |
|-------------|------|-----------|---------------------------|
| `--name`          | `str`   | Yes      | Unique name of protein. |
| `--prot_file`     | `str`   | Yes      | Path to either a PDB file or a TXT file with a single sequence.   |
| `--res_idx`       | `int`   | Yes      | Index of the residue to be covalently targeted by a covalent ligand. Starting at 1.   |         
| `--lig_chain`    | `str`   | Yes      | Chain interacting with ligand in PDB file. Single character.        |                              
| `--lig_csv`       | `bool`  | Yes      | Path to CSV with Ligand info. Do not include any headers. Rows must have `Compound ID,SMILES` info. |
| `--msa_path` | `str` | No | Path to MSA file in csv format. If not provided, msa server will be used. |
| `--ccd_db`       | `str`    | No      | Path to output CSV. Will be formatted to work with setup_cov_job.py. Path to directory with covalent compound pkl files. Default: `/home/$USER/.boltz/mols`.        |
| `--slurm`         | `str`    | No        | Path to SLURM template file. Default: `covalent_module/covalent_slurm.sh` |
| `--debug` | `str` | No | If set, files will not be rewritten, just submitting slurm boltz job. Default=False. | 
| `--outdir`       | `str`    | Yes      | Main output directory for all jobs. |

#### Example Usage
``` bash
python -m covalent_module.submit_jobs -n 5MAJ -p covalent_module/example_input/5MAJ.pdb -r 25 -c A -l covalent_module/example_input/test_ligs.csv -o covalent_module/test/
```

Examples of run input can be found in `covalent_module/example_input`
Examples of run output can be found in `covalent_module/test`

## Input Preprocessing for Covalent Inference
### `preprocessing/make_input_csv.py`

This script sets up a CSV file with necessary information regarding the ligand and protein for writing `.yaml` files using `submit_jobs.py` for covalent docking with Boltz2. Also creates `.pkl` files necessary for docking in default location: `/home/$USER/.boltz/mols/`. `submit_jobs.py` will not work withouth the necessary `.pkl` files and accurately formated CSV for inference. The CCD ID generated when making the `.pkl` files MUST match CCD IDs assigned within the CSV file for inference with `covalent_module/submit_jobs.py`.

#### Command-Line Arguments 

| Argument   | Type | Required | Description               |
|-------------|------|-----------|---------------------------|
| `--name`          | `str`   | Yes      | Name of protein. Used for naming output files. Can be PDB ID, Uniprot ID, anything unique. |
| `--prot_file`     | `str`   | Yes      | Path to either a PDB file or a TXT file with a single chain sequence. Do not include any headers.                                                                                                       |
| `--res_idx`       | `int`   | Yes      | Index of the residue to be covalently targeted by a covalent ligand. Starting at 1. Please confirm index matches expected residue in sequence/PDB provided.                                                       |         
| `--lig_chain`    | `str`   | Yes      | Chain interacting with ligand in PDB file. Single character.        |                              
| `--lig_csv`       | `bool`  | No       | Path to CSV with Ligand info. Do not include any headers with rows `Compound ID,SMILES`. |
| `--out_csv`       | `str`   | Yes      | Path to output CSV. Will be formatted to work with setup_cov_job.py | 

#### Example Usage
For a single protein and several ligands:
```bash
python preprocessing/make_input_csv.py --name 5MAJ --prot_file test/5MAJ.pdb --res_idx 25 --lig_chain A --lig_csv test/test_ligs.csv --out_csv test/test_out.csv --ccd_db /home/ymanasa/.boltz/mols
```

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

Example inputs can be found in `covalent_module/example_input/` and output CSV example is: `covalent_module/test/generated_test_ligs_5MAJ.csv`.