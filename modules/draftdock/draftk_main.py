import argparse
import os
import sys
import yaml
import pandas as pd
import torch
from torch.utils.data import DataLoader, Dataset

# Dynamically set the import path based on the BOLTZ_PATH environment variable
boltz_path = os.getenv("BOLTZ_PATH", "/home/limcaoco/opt/boltz/src")
sys.path.append(boltz_path)

from draftk_train import draftk_train
from boltz.model.modules.diffusionv2 import AtomDiffusion

class CustomDataset(Dataset):
    """
    Custom dataset to handle the specified CSV format.
    """
    def __init__(self, csv_path):
        data = pd.read_csv(csv_path)
        self.names = data["name"].values
        self.fastas = data["FASTA_of_protein"].values
        self.smiles = data["SMILES"].values
        self.pdb_paths = data["path_to_pdb"].values

    def __len__(self):
        return len(self.names)

    def __getitem__(self, idx):
        # Load features and atom mask from the PDB file
        pdb_path = self.pdb_paths[idx]
        feats, atom_mask = self.process_pdb(pdb_path)
        return {"feats": feats, "atom_mask": atom_mask}

    def process_pdb(self, pdb_path):
        # Placeholder for processing PDB files
        # Replace this with actual logic to extract features and atom masks
        feats = torch.randn(100, 3)  # Example: random coordinates
        atom_mask = torch.ones(100)  # Example: all atoms included
        return feats, atom_mask

def load_training_data(csv_path):
    """
    Load training data from a CSV file.

    Parameters:
        csv_path (str): Path to the CSV file.

    Returns:
        DataLoader: A DataLoader containing the training data.
    """
    dataset = CustomDataset(csv_path)
    return DataLoader(dataset, batch_size=32, shuffle=True)

def main():
    parser = argparse.ArgumentParser(description="Fine-tune Boltz diffusion model with DRaFT-K.")
    parser.add_argument("--config", type=str, required=True, help="Path to the YAML configuration file.")
    args = parser.parse_args()

    # Load configuration from YAML file
    with open(args.config, "r") as file:
        config = yaml.safe_load(file)

    # Load training data
    train_dataloader = load_training_data(config["csv"])

    # Initialize base model
    base_model = AtomDiffusion(score_model_args=config["score_model_args"])

    # Call the training function
    draftk_train(
        base_model=base_model,
        train_dataloader=train_dataloader,
        num_epochs=config["num_epochs"],
        lora_rank=config["lora_rank"],
        T=config["T"],
        K=config["K"],
        learning_rate=config["learning_rate"],
        weight_decay=config["weight_decay"],
        device=config["device"],
        wandb_project=config["wandb_project"],
        loss_function=config["loss_function"],
        loss_function_kwargs=config.get("loss_function_kwargs", {}),
    )

if __name__ == "__main__":
    main()
