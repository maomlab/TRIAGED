import torch
import os 
import Vina
import sys
import numpy as np
boltz_path = os.getenv("BOLTZ_PATH")
if boltz_path is None:
    sys.error.write("BOLTZ_PATH environment variable is not set. Please set it to the path of the Boltz source code.\n")
    exit(1)
sys.path.append(boltz_path)
from boltz.data.write.pdb import to_pdb
from boltz.data.types import Coords, Interface, Record, Structure, StructureV2
from boltz.data.write.mmcif import to_mmcif
from dataclasses import asdict, replace
from os import run
def center_of_mass_loss(x_0: torch.Tensor, atom_mask: torch.Tensor, target_point: torch.Tensor) -> torch.Tensor:
    """
    Compute a loss to ensure the ligand's center of mass is near a specified target point.

    Parameters:
        x_0 (torch.Tensor): Tensor of shape (B, N, 3) representing the denoised atomic coordinates.
        atom_mask (torch.Tensor): Tensor of shape (B, N) where 1 indicates ligand atoms and 0 otherwise.
        target_point (torch.Tensor): Tensor of shape (3,) representing the target point.

    Returns:
        torch.Tensor: The computed loss (scalar).
    """
    # Extract ligand coordinates using the atom mask
    ligand_coords = x_0[atom_mask.bool()]

    # Compute the center of mass of the ligand
    center_of_mass = ligand_coords.mean(dim=0)

    # Compute the L2 distance between the center of mass and the target point
    loss = torch.norm(center_of_mass - target_point, p=2)

    return loss

def coord_tensor_to_pdbqt(denoised_coords, feats, path, multiplicity=1) -> None:
    """
    convert boltz diffusion tensor into pdbqt
    """
    pad_token_mask = feats["token_pad_mask"].repeat_interleave(multiplicity, 0)
    
    '''
    rec_mask = (feats["mol_type"] == 0).repeat_interleave(multiplicity, 0)
    rec_mask = rec_mask * pad_token_mask
    lig_mask = (
        feats["affinity_token_mask"]
        .repeat_interleave(multiplicity, 0)
        .to(torch.bool)
    )
    lig_mask = lig_mask * pad_token_mask
    cross_pair_mask = (
        lig_mask[:, :, None] * rec_mask[:, None, :]
        + rec_mask[:, :, None] * lig_mask[:, None, :]
        + lig_mask[:, :, None] * lig_mask[:, None, :]
    )
    '''
    structure: StructureV2 = StructureV2.load(path)
    structure = structure.remove_invalid_chains()

    chain_map = {}

    for i, mask in enumerate(structure.mask):
        if mask:
            chain_map[len(chain_map)] = i

    denoised_coords = denoised_coords.unsqueeze(0)  
    coord_unpad = denoised_coords[pad_token_mask.bool()]
    coord_unpad = coord_unpad.cpu().numpy()

    atoms = structure.atoms
    atoms["coords"] = coord_unpad
    atoms["is_present"] = True

    structure: StructureV2
    coord_unpad = [(x,) for x in coord_unpad]
    coord_unpad = np.array(coord_unpad, dtype=Coords)

    residues = structure.residues
    residues["is_present"] = True
    interfaces = np.array([], dtype=Interface)

    new_structure: StructureV2 = replace(
    structure,
    atoms=atoms,
    residues=residues,
    interfaces=interfaces,
    coords=coord_unpad,
    )

    chain_info = []
    for chain in new_structure.chains:
        old_chain_idx = chain_map[chain["asym_id"]]
        old_chain_info = record.chains[old_chain_idx]
        new_chain_info = replace(
            old_chain_info,
            chain_id=int(chain["asym_id"]),
            valid=True,
        )
        chain_info.append(new_chain_info)

    # Save the structure
    project_path = os.getenv("PROJECT_PATH")
    if project_path is None:
        sys.stderr.write("PROJECT_PATH environment variable is not set. Please set it to the project directory.\n")
        exit(1)

    temp_dir = project_path / "temp"
    temp_dir.mkdir(exist_ok=True)
    
    temp_path = temp_dir / "temp.pdb"
    with temp_path.open("w") as f:
        f.write(
            to_pdb(new_structure, boltz2=True)
        )
    
    # Convert the PDB file to PDBQT using Open Babel
    temp_pdbqt_path = os.path.join(temp_dir, "temp.pdbqt")
    command = ["obabel", temp_pdb_path, "-O", temp_pdbqt_path, "--partialcharge", "gasteiger"]
    result = run(command, capture_output=True, text=True)
    if result == None:
        sys.stderr.write(f"Error converting PDB to PDBQT: {result.stderr}\n")
        exit(1)
    # Save the PDBQT file
    #temp_path = temp_dir / "temp.pdbqt"
    #with open(temp_path, "w") as f:
    #    f.write(result.stdout)


    return result.stdout  

def vina_scoring_function(denoised_coords: torch.Tensor, feats) -> torch.Tensor:
    """
    Compute binding affinity using AutoDock Vina.

    Parameters:
        ligand_coords (torch.Tensor): Tensor of shape (N, 3) representing ligand atomic coordinates.
        receptor_path (str): Path to the receptor PDBQT file.
        box_center (list): Center of the docking box [x, y, z].
        box_size (list): Size of the docking box [x, y, z].

    Returns:
        torch.Tensor: Binding affinity (negative value, lower is better).
    """
    # Convert ligand coordinates to PDBQT format
    complex_pdbqt = coord_tensor_to_pdbqt(denoised_coords, feats)

    # Initialize Vina
    v = Vina(sf_name='vina')

    # Set receptor and docking box
    v.set_receptor(receptor_path)
    v.set_box(center=box_center, size=box_size)

    # Set ligand
    v.set_ligand_from_file(ligand_pdbqt)
    # Split the complex PDBQT into receptor and ligand PDBQT files
    receptor_pdbqt = "receptor_temp.pdbqt"
    ligand_pdbqt = "ligand_temp.pdbqt"

    # Extract receptor and ligand from the complex PDBQT
    with open(complex_pdbqt, "r") as complex_file:
        lines = complex_file.readlines()

    with open(receptor_pdbqt, "w") as receptor_file, open(ligand_pdbqt, "w") as ligand_file:
        is_receptor = True
        for line in lines:
            if line.startswith("MODEL") or line.startswith("ENDMDL"):
                continue
            if line.startswith("REMARK VINA RESULT"):
                is_receptor = False
            if is_receptor:
                receptor_file.write(line)
            else:
                ligand_file.write(line)

    # Use the separated ligand PDBQT file
    ligand_pdbqt

    # Perform docking
    v.dock()
    affinity = v.score()  # Binding affinity in kcal/mol

    # Convert to PyTorch tensor
    return torch.tensor(affinity, requires_grad=True)

