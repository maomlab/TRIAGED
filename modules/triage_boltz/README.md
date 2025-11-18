
### `setup_boltz_job.py`
This script sets up Boltz job directories and generates `.yaml` files based on the input CSV and biomolecular data. It supports both covalent and non-covalent docking workflows and handles complex systems with multiple biomolecules.

#### Command-Line Arguments

| Argument               | Type    | Required | Description                                                                                                   |
|------------------------|---------|----------|---------------------------------------------------------------------------------------------------------------|
| `--input_csv_file`     | `str`   | Yes      | Path to the input CSV file containing compound information.              |
| `--input_pdb_file`     | `str`   | No       | Path to the input PDB file containing protein structure.                                                      |
| `--output_directory`   | `str`   | Yes      | Path to the output directory where `.yaml` files will be created.                                             |
| `--input_fasta_file`   | `str`   | No       | Path to the input fasta file containing sequence, if inputed will be used over the pdb file                   |
| `--num_jobs`           | `int`   | No       | How may parallel submission directories and slurm submit scripts to generate depending on avaliable gpus      |
| `--protein_nmers`      | `int`   | No       | How many subunits to model the receptor if it is a multimer                                                   |

#### Example Usage

For a simple virtual screen:
```bash
<<<<<<< HEAD:scripts/README.md
python setup_boltz_job.py --input_csv_file /path/to/input.csv --input_pdb_file /path/to/input.pdb --output_directory /path/to/output_directory
```
optionally you can forgo the pdb file and just give the fasta file when using non-covelant docking:
```bash
python setup_boltz_job.py --input_csv_file /path/to/input.csv --input_fasta_file /path/to/input.fasta --output_directory /path/to/output_directory
```

For a simple virtual screen on 4 gpus:
```bash
python setup_boltz_job.py --input_csv_file /path/to/input.csv --input_pdb_file /path/to/input.pdb --output_directory /path/to/output_directory --num_jobs 4
```
For a virtual screen on a tetramer protein on 4 gpus:
```bash
python setup_boltz_job.py --input_csv_file /path/to/input.csv --input_pdb_file /path/to/input.pdb --output_directory /path/to/output_directory --num_jobs 4 --protein_nmers
```
An example real life command:
```bash
python setup_boltz_job.py -i ../input_files/activity_data/kcnq_compounds.csv -f ../input_files/fastas/KCNQ2.fasta -o ../boltz_inputs/KCNQ2 -n 4 --protein_nmers 4
```
=======
python setup_boltz_job.py --input_csv_file /path/to/input.csv --output_directory /path/to/output_directory --name receptor_name
```

For a virtual screen on 4 GPUs:
```bash
python setup_boltz_job.py --input_csv_file /path/to/input.csv --output_directory /path/to/output_directory --num_jobs 4 --name receptor_name
```

For covalent docking:
```bash
python setup_boltz_job.py --input_csv_file --output_directory /path/to/output_directory --covalent_docking --name receptor_name
```
Note for covalent docking the input will be a csv with 
>>>>>>> origin:modules/triage_boltz/README.md

#### Input CSV Requirements

The input CSV file must contain the following columns:
<<<<<<< HEAD:scripts/README.md
- `compound_id`: Unique identifier for the compound.
- `SMILES`: SMILES string representing the compound structure.

#### Output
- `.a3m` MSA files are pregenerated using the mmseqs2 server, you do not need to have --use-msa-server flag in boltz
- `.yaml` files are created in the specified output directory, one for each compound in the CSV file.
- The protein sequence is extracted from the PDB file and included in the `.yaml` files.

---
## SLURM job scripts 
### slurm_run_boltz.sh

The script will automatically generate SLURM batch scripts for running Boltz jobs on a high-performance computing cluster with the pathing and inputs for boltz already set 
to the generated input directories. By default these scripts will be generated in {Project_dir}/slurm_scripts

#### Example Usage
=======
- `compound_ID`: Unique identifier for the compound. the first compound_ID is what will be the binder.
- `SMILES`: SMILES string representing the compound structure.
- `protein_ID`, `dna_ID`, `rna_ID`: Optional columns for specifying proteins, DNA, or RNA.
- `num`: Optional column specifying the number of entities to model (default is 1).

Example CSV Format for Complex Systems

Below is an example of how to format a CSV file for setting up a virtual screen with additional biomolecules (e.g., cholesterol, ATP, and peptides):

```csv
compound_ID_1,compound_num_1,SMILES_1,compound_ID_2,compound_num_2,SMILES_2,protein_ID_1,protein_sequence_1,protein_ID_2,protein_sequence_2,
ZINC0001,1,CC(=O)Nc1ccc(C)cc1,cholesterol,4,C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C,test_receptor,APQQINDIVHRTITPLIEQQKIPGMAVAVIYQGKPYYFTWGYADIAKKQPVTQQTLFELGSVSKTFTGVLGGDAIARGEIKLSDPTTKYWPELTAKQWNGITLLHLATYTAGGLPLQVPDEVKSSSDLLRFYQNWQPAWAPGTQRLYANSSIGLFGALAVKPSGLSFEQAMQTRVFQPLKLNHTWINVPPAEEKNYAWGYREGKAVHVSPGALDAEAYGVKSTIEDMARWVQSNLKPLDINEKTLQQGIQLAQSRYWQTGDMYQGLGWEMLDWPVNPDSIINGSDNLAARPVKAITPPTPAVRASWVHKTGATGGFGSYVAFIPEKELGIVMLANKNYPNPARVDAAWQILNALQAPQQINDIVHRTITPLIEQQKIPGMAVAVIYQGKPYYFTWGYADIAKKQPVTQQTLFELGSVSKTFTGVLGGDAIARGEIKLSDPTTKYWPELTAKQWNGITLLHLATYTAGGLPLQVPDEVKSSSDLLRFYQNWQPAWAPGTQRLYANSSIGLFGALAVKPSGLSFEQAMQTRVFQPLKLNHTWINVPPAEEKNYAWGYREGKAVHVSPGALDAEAYGVKSTIEDMARWVQSNLKPLDINEKTLQQQGIQLAQSRYWQTGDMYQGLGWEMLDWPVNPDSIINGSDNKIALAARPVKAITPPTPAVRASWVHKTGATGGFGSYVAFIPEKELGIVMLANKNYPNPARVDAAWQILNALQ,ubiq,MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG
ZINC0002,1,CC(=O)Nc1ccc(C)cc1,cholesterol,4,C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C,test_receptor,APQQINDIVHRTITPLIEQQKIPGMAVAVIYQGKPYYFTWGYADIAKKQPVTQQTLFELGSVSKTFTGVLGGDAIARGEIKLSDPTTKYWPELTAKQWNGITLLHLATYTAGGLPLQVPDEVKSSSDLLRFYQNWQPAWAPGTQRLYANSSIGLFGALAVKPSGLSFEQAMQTRVFQPLKLNHTWINVPPAEEKNYAWGYREGKAVHVSPGALDAEAYGVKSTIEDMARWVQSNLKPLDINEKTLQQGIQLAQSRYWQTGDMYQGLGWEMLDWPVNPDSIINGSDNLAARPVKAITPPTPAVRASWVHKTGATGGFGSYVAFIPEKELGIVMLANKNYPNPARVDAAWQILNALQAPQQINDIVHRTITPLIEQQKIPGMAVAVIYQGKPYYFTWGYADIAKKQPVTQQTLFELGSVSKTFTGVLGGDAIARGEIKLSDPTTKYWPELTAKQWNGITLLHLATYTAGGLPLQVPDEVKSSSDLLRFYQNWQPAWAPGTQRLYANSSIGLFGALAVKPSGLSFEQAMQTRVFQPLKLNHTWINVPPAEEKNYAWGYREGKAVHVSPGALDAEAYGVKSTIEDMARWVQSNLKPLDINEKTLQQQGIQLAQSRYWQTGDMYQGLGWEMLDWPVNPDSIINGSDNKIALAARPVKAITPPTPAVRASWVHKTGATGGFGSYVAFIPEKELGIVMLANKNYPNPARVDAAWQILNALQ,ubiq,MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG
ZINC0003,1,CC(=O)Nc1ccc(C)cc1,cholesterol,4,C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C,test_receptor,APQQINDIVHRTITPLIEQQKIPGMAVAVIYQGKPYYFTWGYADIAKKQPVTQQTLFELGSVSKTFTGVLGGDAIARGEIKLSDPTTKYWPELTAKQWNGITLLHLATYTAGGLPLQVPDEVKSSSDLLRFYQNWQPAWAPGTQRLYANSSIGLFGALAVKPSGLSFEQAMQTRVFQPLKLNHTWINVPPAEEKNYAWGYREGKAVHVSPGALDAEAYGVKSTIEDMARWVQSNLKPLDINEKTLQQGIQLAQSRYWQTGDMYQGLGWEMLDWPVNPDSIINGSDNLAARPVKAITPPTPAVRASWVHKTGATGGFGSYVAFIPEKELGIVMLANKNYPNPARVDAAWQILNALQAPQQINDIVHRTITPLIEQQKIPGMAVAVIYQGKPYYFTWGYADIAKKQPVTQQTLFELGSVSKTFTGVLGGDAIARGEIKLSDPTTKYWPELTAKQWNGITLLHLATYTAGGLPLQVPDEVKSSSDLLRFYQNWQPAWAPGTQRLYANSSIGLFGALAVKPSGLSFEQAMQTRVFQPLKLNHTWINVPPAEEKNYAWGYREGKAVHVSPGALDAEAYGVKSTIEDMARWVQSNLKPLDINEKTLQQQGIQLAQSRYWQTGDMYQGLGWEMLDWPVNPDSIINGSDNKIALAARPVKAITPPTPAVRASWVHKTGATGGFGSYVAFIPEKELGIVMLANKNYPNPARVDAAWQILNALQ,ubiq,MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG
```

Example Rows:
Row 1: Includes ZINC0001 (compound 1), cholesterol (compound 2), test_receptor (protein 1) ubiquitin (protein 2)
Row 2: Includes ZINC0002 (compound 1), cholesterol (compound 2), test_receptor (protein 1) ubiquitin (protein 2)

This format allows you to specify multiple compounds and proteins to be docked to complex systems in a single CSV file. The script will also generate the .yaml inputs in convient number parallel batch directories with associated slurm submit scripts with the use of the --num_jobs flag if you so desire. As of right now, the affinity will be only computed for the first compound i.e compound_ID_1 given, though by default confidence model scores will be generated for all biomolecules

if you wish to do covalent docking, you must pass a csv file with the protein_ID and pdb_path columns like so:
```csv
protein_ID,pdb_path
aa2a,/some/path
a2ba,/some/path
2web,/some/path
```

#### Features

- **CSV Parsing**: Reads input CSV files to create `TriageBiomolecule` objects for ligands, proteins, DNA, and RNA.
- **FASTA Generation**: Automatically generates FASTA sequences from PDB files or fetches them using PDB IDs if missing.
- **MSA Computation**: Computes single and paired MSAs for proteins and creates combined MSAs.
- **YAML File Generation**: Creates YAML configuration files for Boltz, including sequences, MSAs, and properties.
- **Covalent Docking Support**: Processes PDB files to extract covalent interaction details for ligand-protein docking.
- **Parallel Job Setup**: Splits input data into chunks for parallel job execution.

#### Output

For non-covalent docking:
- `.a3m` MSA files are pregenerated using the mmseqs2 server.
- `.yaml` files are created in the specified output directory, one for each compound in the CSV file.

For covalent docking:
- `.yaml` files are created in the specified output directory.
- The protein sequence is extracted from the PDB file and included in the `.yaml` files.
- Covalent bond information is automatically extracted and included in the `.yaml` files.

#### SLURM Scripts

The script automatically generates SLURM batch scripts for running Boltz jobs on a high-performance computing cluster. These scripts are saved in `{PROJECT_DIR}/slurm_scripts`.

#### Example SLURM Script Usage
>>>>>>> origin:modules/triage_boltz/README.md

```bash
sbatch slurm_run_boltz.sh
```

<<<<<<< HEAD:scripts/README.md
---

## Future Features: aka work in progress

- generate_boltz_constrains: a method to automatically generate constraints to condition diffusion towards specific protein active states based on two input pdbs.


=======
#### Notes

- Ensure the `PROJECT_DIR` environment variable is set before running the script.
- Input CSV files must follow the expected format, with columns for `compound_ID`, `protein_ID`, `dna_ID`, `rna_ID`, and their associated properties.
- For covalent docking, only PDB files are required, and the `--use-msa-server` flag must be used when running Boltz.

---

>>>>>>> origin:modules/triage_boltz/README.md
# Analyze Boltz Predictions

This script processes Boltz predictions, computes metrics, and optionally performs bootstrapping for further analysis. The results are saved as CSV and JSON files for easy access.

## Requirements

- Python 3.x
- Required Python libraries:
  - `pandas`
  - `argparse`
  - `json`
  - `scoring_utils` (custom module)
  - Other dependencies listed in `scoring_utils.py`

- Use the provided `analysis_env.yaml` to set up the conda environment for dependencies.

## Usage

### Command-Line Arguments

| Argument                   | Type    | Required | Default | Description                                                                |
|----------------------------|---------|----------|---------|----------------------------------------------------------------------------|
| `-i`, `--input_directory`  | `str`   | Yes      | None    | Directory containing Boltz predictions.                                    |
| `-o`, `--output_directory` | `str`   | Yes      | None    | Directory to save the processed output files.                              |
| `-m`, `--compute_metrics`  | `bool`  | No       | `False` | Whether to compute metrics or not.                                         |
| `-b`, `--bootstrap`        | `bool`  | No       | `False` | Whether to compute bootstrap metrics or not.                               |
| `-v`, `--in-vitro-data`    | `str`   | No       | `None`  | Path to the in vitro data file, requires a column named `is_binder`.       |

### Example Usage

#### Basic Processing
```bash
python analyize_boltz_predictions.py -i /path/to/input_directory -o /path/to/output_directory
```

#### Compute Metrics
```bash
python analyize_boltz_predictions.py -i /path/to/input_directory -o /path/to/output_directory -m True -v /path/to/in_vitro_data.csv
```

#### Compute Metrics with Bootstrapping
```bash
python analyize_boltz_predictions.py -i /path/to/input_directory -o /path/to/output_directory -m True -b True -v /path/to/in_vitro_data.csv
```

### Output Files

1. **Processed Data**:  
   - File: `processed_boltz_data_full.csv`  
   - Description: Contains the processed Boltz predictions.

2. **Computed Metrics**:  
   - File: `computed_metrics.json`  
   - Description: Contains computed metrics and bootstrap metrics (if enabled) for each score column.
3. **Histograms**:
   - File: `{metric}_bootstrap_histogram_{score}.png`  
   - Description: script will automatically generate bootstraped histograms for each score column,. 
   
### Input Data Requirements

#### Boltz Predictions
The input directory should contain subdirectories with JSON files for each compound. Each subdirectory must include:
- `affinity_<compound_name>.json`
- `confidence_<compound_name>_model_0.json`

#### In Vitro Data
The in vitro data file must be a CSV with a column named `is_binder` indicating whether a compound is a binder (`True`) or not (`False`).

### Notes

- Ensure the input directory and in vitro data file are correctly formatted.
- Bootstrapping may take longer depending on the dataset size and number of iterations.
___ 

# Combine Boltz model outputs for visualization
### `compile_best_model_structures.py`

This script aligns substructures from multiple Boltz2 job directories to a given
 parent structure and saves the results as a multi-model PDB file. It assumes 
 that all run Boltz2 jobs are contained in a single directory with predictions 
 located in subdirectories.

By default, if a parent is specified, it will be the first model in the written PDB.
 If the `-ep` flag is specified, the PARENT structure will be excluded.

#### Command-Line Arguments

| Argument               | Type    | Required | Description                                                                 |
|------------------------|---------|----------|-----------------------------------------------------------------------------|
| `-p`, `--parent`       | `str`   | No       | Path to the parent file (CIF, PDB, or MOL2) to align to. Defaults to `first._0.cif` in predictions directory. |
| `-d`, `--directory`    | `str`   | Yes      | Root directory containing 'predictions' subdirectories with CIF files.     |
| `-o`, `--output`       | `str`   | No       | Name for the output multi-model PDB file. Defaults to `aligned_models.pdb`.|
| `--max-models`         | `int`   | No       | Maximum number of models to include in the output PDB. Defaults to all found.|
| `-ep`, `--exclude_parent`|       | No       | If a PARENT is specified and flag is used, parent will not be saved to final output. |
| `-v`, `--verbose`      |         | No       | Enable verbose output.                                                      |

#### Example Usage

```bash
python compile_best_model_structures.py -p 3JQZ.pdb -d boltz_jobs -o 3JQZ_top_models.pdb
```

### `compile_individual_mol_structures.py`

This script aligns substructures from a single Boltz2 job directory to a 
 specified parent structure and saves the results as a multi-model PDB file. If
 no parent is given, will default to the first modeled structure.

By default, if a parent is specified, it will be the first model in the written PDB.
 If the `-ep` flag is specified, the PARENT structure will be excluded.

#### Command-Line Arguments

| Argument               | Type    | Required | Description                                                                 |
|------------------------|---------|----------|-----------------------------------------------------------------------------|
| `-p`, `--parent`       | `str`   | No       | Path to the parent file (CIF, PDB, or MOL2) to align to. Defaults to the first structure in the directory. |
| `-d`, `--directory`    | `str`   | Yes      | Directory containing CIF substructures from a Boltz2 job.                  |
| `-o`, `--output`       | `str`   | No       | Name for the output multi-model PDB file. Defaults to `aligned_models.pdb`.|
| `--max-models`         | `int`   | No       | Maximum number of models to include in the output PDB. Defaults to all found in the directory. |
| `-ep`, `--exclude_parent`|       | No       | If a PARENT is specified and flag is used, parent will not be saved to final output. |
| `-v`, `--verbose`      |         | No       | Enable verbose output.                                                      |

#### Example Usage

```bash
<<<<<<< HEAD:scripts/README.md
<<<<<<< HEAD:scripts/README
=======
<<<<<<< HEAD:scripts/README.md
>>>>>>> origin:modules/triage_boltz/README.md
python compile_individual_mol_structures.py -p 3JQZ.pdb -d boltz_jobs/mol1/boltz_results_mol1/predictions/mol1 -o 3JQZ_top_models.pdb```
=======
python compile_individual_mol_structures.py -p 3JQZ.pdb -d boltz_jobs/mol1/boltz_results_mol1/predictions/mol1 -o 3JQZ_top_models.pdb
```
<<<<<<< HEAD:scripts/README.md
>>>>>>> master:scripts/README.md
=======
>>>>>>> origin:modules/triage_boltz/README.md
>>>>>>> origin:modules/triage_boltz/README.md



### Contact

<<<<<<< HEAD:scripts/README.md
For issues or questions, please contact limcaoco@umich.edu
=======
For issues or questions, please contact limcaoco@umich.edu
>>>>>>> origin:modules/triage_boltz/README.md
