from dataclasses import dataclass, field
from typing import Optional, List
import os
from rdkit import Chem
from triage_boltz.misc_utils import parse_input_csv, generate_msa, sanitize_compound_id, check_smiles
import torch

from triage_boltz.pdb_to_fasta import pdb_to_fasta
from triage_boltz.fasta_utils import build_fasta_seq
@dataclass
class TriageBiomolecule:
    """
    A data class to store information about intermolecular entities.
    needs to be compatable with all modules in the triage pipeline.
    """
    entity_id: str
    entity_type: str  # e.g., 'ligand', 'protein', 'dna', 'rna'
    smiles: Optional[str] = None
    inchi: Optional[str] = None
    num: Optional[int] = 1  # Number of entities to model, e.g., for ligands or proteins default is 1. 
    sequence: Optional[str] = None
    pdb_path: Optional[str] = None
    pdb_id: Optional[str] = None
    msa_path: Optional[str] = None
    pair_msa_path: Optional[str] = None
    combined_msa_path: Optional[str] = None
    pose_path: Optional[str] = None

    #fingerprints: 
    ecfp: Optional[torch.Tensor] = None # ECFP fingerprint
    prolif_count: Optional[torch.Tensor] = None # ProLIF count vector
    prolif_bit: Optional[torch.Tensor] = None # ProLIF bit vector

    # add class variables for boltz affinity- to be connected in the pipline
    #
    #

    # add class variables for clustering to be connected in the pipline.

    def __post_init__(self):
        """
        Validate the inputs after initialization.
        """
        #if self.coordinates is not None and not isinstance(self.coordinates, torch.Tensor):
        #    raise TypeError("Coordinates must be a torch.Tensor.")
        if self.smiles is None and self.sequence is None and self.pdb_path is None and self.inchi is None:
            raise ValueError("At least one of 'smiles', 'sequence', 'pdb_path', or 'inchi' must be provided.")
        if self.entity_type == "ligand":
            self.sanitize_entity_id()
            if self.smiles:
                self.check_smiles()
                self.generate_inchi_from_smiles()
            elif self.inchi:
                self.generate_smiles_from_inchi()

    def sanitize_entity_id(self) -> None:
        """
        Sanitize the entity ID by removing any non-alphanumeric
        characters and converting it to uppercase.
        """
        self.entity_id = sanitize_compound_id(self.entity_id)


    def check_smiles(self, verbose: bool = True) -> None:
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
            mol = Chem.MolFromSmiles(self.smiles, sanitize=False)
            if mol is None:
                if verbose:
                    print(f"[ERROR] MolFromSmiles failed for: {self.smiles}")
                return None
            # Attempt sanitization (includes valence check, aromaticity, Hs)
            Chem.SanitizeMol(mol)
            # Return the canonical SMILES
            self.smiles = Chem.MolToSmiles(mol, canonical=True)
        except Exception as e:
            if verbose:
                print(f"[ERROR] Sanitization failed for SMILES: {self.smiles}\n{e}")
            return None
    
    def generate_inchi_from_smiles(self) -> None:
        """
        Generate an InChI string from the SMILES string and store it in the class.
        """
        if self.smiles is None:
            raise ValueError("SMILES string is not provided.")
        mol = Chem.MolFromSmiles(self.smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string, cannot generate InChI.")
        self.inchi = Chem.MolToInchi(mol)

    def generate_smiles_from_inchi(self) -> None:
        """
        Generate a SMILES string from the InChI string and store it in the class.
        """
        if self.inchi is None:
            raise ValueError("InChI string is not provided.")
        mol = Chem.MolFromInchi(self.inchi)
        if mol is None:
            raise ValueError("Invalid InChI string, cannot generate SMILES.")
        self.smiles = Chem.MolToSmiles(mol)

    def generate_single_msa(self) -> None:
        """
        Generate a single Multiple Sequence Alignment (MSA) for the biomolecule and store the path in msa_path.
        """
        if self.entity_type != "protein":
            raise ValueError("MSA generation is only supported for proteins.")
        
        if self.sequence is not None:
            msa_file = generate_msa(self.entity_id, self.sequence)
            self.msa_path = msa_file
        else:
            raise ValueError("Sequence must be provided to generate MSA.")

    def generate_pair_msa(self, other_biomolecules: 'TriageBiomolecule' | List['TriageBiomolecule']) -> None:
        """
        Generate a pair Multiple Sequence Alignment (MSA) for the biomolecule with another entity or entities and store the path in pair_msa_path.

        Parameters:
        - other_biomolecules (TriageBiomolecule | List[TriageBiomolecule]): Another TriageBiomolecule instance or a list of TriageBiomolecule instances to align with.
        """
        if self.entity_type != "protein":
            raise ValueError("Pair MSA generation is only supported for proteins.")
        
        if isinstance(other_biomolecules, TriageBiomolecule):
            other_biomolecules = [other_biomolecules]
        elif not isinstance(other_biomolecules, list) or not all(isinstance(e, TriageBiomolecule) for e in other_biomolecules):
            raise TypeError("other_biomolecules must be a TriageBiomolecule instance or a list of TriageBiomolecule instances.")
        
        sequences = [self.sequence] + [biomolecule.sequence for biomolecule in other_biomolecules if biomolecule.sequence is not None]
        if not sequences:
            raise ValueError("At least one sequence must be provided for pair MSA generation.")
        #print("these are the sequences to be aligned:")
        #print(sequences)
        
        pair_msa_file = generate_msa(self.entity_id, sequences)
        self.pair_msa_path = pair_msa_file
    

    def create_combined_msa(self) -> None:
        """
        creates a combined MSA from the single and pair MSAs and stores the path of the resulting MSA in combined_msa_path.
        """
        print("self.msa_path:", self.msa_path)
        print("self.pair_msa_path:", self.pair_msa_path)
        if self.msa_path is None or self.pair_msa_path is None:
            raise ValueError("Single MSA and paired MSA must be generated before creating a combined MSA.")
        if not os.path.exists(self.msa_path):
            raise FileNotFoundError(f"Single MSA file not found at {self.msa_path}")
        if not os.path.exists(self.pair_msa_path):  
            raise FileNotFoundError(f"Pair MSA file not found at {self.pair_msa_path}")
        
        with open(self.msa_path, 'r') as msa_file:
            unpaired = [line.strip() for line in msa_file if line.strip() and not line.startswith(">")]
        with open(self.pair_msa_path, 'r') as pair_msa_file:
            paired = [line.strip() for line in pair_msa_file if line.strip() and not line.startswith(">")]

        seqs = paired + unpaired
        keys = list(range(len(paired))) + [-1] * len(unpaired)
        
        project_dir = os.getenv("PROJECT_DIR")
        msa_dir = os.path.join(project_dir, "input_files/msa/")
        combined_msa_path = os.path.join(msa_dir, f"{self.entity_id}_combined_msa.csv")

        with open(combined_msa_path,"w") as f:
            f.write("key,sequence\n")
            for key, seq in zip(keys, seqs):
                f.write(f"{key},{seq}\n")

        print(f"Combined MSA written to {combined_msa_path}")
        self.combined_msa_path = combined_msa_path  # Store the path in the instance


    def generate_fasta_from_pdb_path(self) -> None:
        """
        generate a fasta sequence from a PDB file and returns the fasta sequence as a string.
        """
        if self.entity_type != "protein":
            raise ValueError("FASTA generation is only supported for proteins.")

        if self.pdb_path is None:
            raise ValueError("PDB path is not provided.")
        
        if not os.path.exists(self.pdb_path):
            raise FileNotFoundError(f"PDB file not found at {self.pdb_path}")
        
        fasta_sequence = pdb_to_fasta(self.pdb_path)
        if not fasta_sequence:
            raise ValueError("Failed to generate FASTA sequence from PDB file.")
        self.sequence = fasta_sequence  # Store the sequence in the instance
    
    def fetch_fasta_from_pdb_id(self) -> None:
        """
        Fetch FASTA sequence from PDB ID and store it in the instance.
        
        Returns:
        - str: The FASTA sequence.
        """
        if self.entity_type != "protein":
            raise ValueError("FASTA fetching is only supported for proteins.")
        
        if not self.entity_id:
            raise ValueError("Entity ID is not provided.")
        
        fasta_sequence = build_fasta_seq(self.entity_id)
        if not fasta_sequence:
            raise ValueError(f"Failed to fetch FASTA sequence for PDB ID {self.entity_id}.")
        self.sequence = fasta_sequence

    @staticmethod
    def from_csv(csv_file: str) -> list[list['TriageBiomolecule']]:
        """
        Create a list of lists of TriageBiomolecule instances from a CSV file using the parse_input_csv function.

        Parameters:
        - csv_file (str): Path to the CSV file.

        Returns:
        - List[List[TriageBiomolecule]]: A list of lists of TriageBiomolecule instances.
        """
        parsed_data = parse_input_csv(csv_file)
        biomolecules = []

        for row in parsed_data:
            row_biomolecules = []
            for entry in row:
                #print(f"Processing entry: {entry}")
                try:
                    if "compound_ID" in entry:
                        row_biomolecules.append(
                            TriageBiomolecule(
                                entity_id=entry["compound_ID"],
                                entity_type="ligand",
                                smiles=entry.get("SMILES"),
                                inchi=entry.get("InChI"),
                                num=int(entry.get("compound_num")),
                                pose_path=entry.get("pose_path"),
                            )
                        )
                    elif "protein_ID" in entry:
                        row_biomolecules.append(
                            TriageBiomolecule(
                                entity_id=entry["protein_ID"],
                                entity_type="protein",
                                sequence=entry.get("protein_sequence"),
                                pdb_path=entry.get("pdb_path"),
                                pdb_id=entry.get("pdb_id"),
                                num=int(entry.get("protein_num"))
                            )
                        )
                    elif "dna_ID" in entry:
                        row_biomolecules.append(
                            TriageBiomolecule(
                                entity_id=entry["dna_ID"],
                                entity_type="dna",
                                sequence=entry.get("dna_sequence"),
                                num=int(entry.get("dna_num"))
                            )
                        )
                    elif "rna_ID" in entry:
                        row_biomolecules.append(
                            TriageBiomolecule(
                                entity_id=entry["rna_ID"],
                                entity_type="rna",
                                sequence=entry.get("rna_sequence"),
                                num=int(entry.get("rna_num"))
                            )
                        )
                
                except Exception as e:
                    print(f"[ERROR] Failed to create TriageBiomolecule from entry {entry}: {e}")
                    # Skip the entire row if any entry is invalid
                    row_biomolecules = []
                    break
            if row_biomolecules:
                biomolecules.append(row_biomolecules)

        return biomolecules
