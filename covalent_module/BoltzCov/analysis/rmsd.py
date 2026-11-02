
def get_distance(conf, rdkit_idx1, rdkit_idx2):
    """
    Measures the Euclidean distance between two atoms in a combined protein-ligand RDKit conformer.
    
    :param conf (rdkit.Chem.rdchem.Conformer): RDKit combined conformer object containing 3D coordinates.
    :param rdkit_prot_idx (int): RDKit atom index of the protein atom.
    :param rdkit_lig_idx (int): RDKit atom index of the ligand atom.

    :return float: Distance between the two atoms in Ångströms.
    
    Notes:
        - Use find_atom_index() to get rdkit idxs. 
        - Rdkit idxs are different from atom idxs listed in PDB files. 
    """
    return rdMolTransforms.GetBondLength(conf, rdkit_idx1, rdkit_idx2)

def get_ca_atoms_by_residues(structure):
    """Return a dictionary {(chain_id, resid): atom} for CA atoms"""
    ca_atoms = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    resid = residue.get_id()[1]  
                    chain_id = chain.id.strip()
                    key = (chain_id, resid)
                    ca_atoms[key] = residue['CA']
    return ca_atoms

def compute_rmsd(ref_path, pred_path):
    """
    Computes the Root-Mean-Square Deviation (RMSD) between two PDB structures 
    using C-alpha atoms of common residues.

    Args:
        ref_path (str): File path to the reference PDB structure.
        pred_path (str): File path to the predicted PDB structure.

    Returns:
        float: The RMSD value in Ångstroms if common residues are found.
        None: If no common residues are found between the structures.

    Notes:
        - Only residues with matching chain ID and residue number are compared.
        - Superimposes predicted structure onto the reference before computing RMSD.
    """
    parser = PDBParser(QUIET=True)
    ref_struct = parser.get_structure("ref", ref_path)
    pred_struct = parser.get_structure("pred", pred_path)

    ref_atoms_map = get_ca_atoms_by_residues(ref_struct)
    pred_atoms_map = get_ca_atoms_by_residues(pred_struct)

    # Get common residues
    common_keys = sorted(set(ref_atoms_map) & set(pred_atoms_map))
    print(f"Found {len(common_keys)} common residues for RMSD.")

    if len(common_keys) == 0:
        print("No overlapping residues found.")
        return

    ref_atoms = [ref_atoms_map[key] for key in common_keys]
    pred_atoms = [pred_atoms_map[key] for key in common_keys]

    # Superimpose and calculate RMSD
    sup = Superimposer()
    sup.set_atoms(ref_atoms, pred_atoms)
    print(f"RMSD: {sup.rms:.3f} Å")

    return sup.rms

def get_neighbor(atom, exclude_idx):
    """
    Returns the index of a neighboring atom, prioritizing atoms with 'C' in their PDB name,
    while excluding a specified atom index.

    Args:
        atom (rdkit.Chem.Atom): RDKit Atom object whose neighbors are to be examined.
        exclude_idx (int): Atom index to be excluded (e.g., the atom it's bonded to for a dihedral).

    Returns:
        int or None: Index of the selected neighboring atom, or None if no suitable neighbor is found.

    Notes:
        - Only neighbors with 'C' in their PDB atom name are considered (e.g., CB, CG).
        - This is a heuristic to prioritize carbon atoms in dihedral definitions.
    """
    for neighbor in atom.GetNeighbors():
        pdb_info = neighbor.GetPDBResidueInfo()
        name = pdb_info.GetName().strip()
        
        if 'C' in name:
            print(f"Neighbor name: {name}")
            if neighbor.GetIdx() != exclude_idx:
                return neighbor.GetIdx()
    return None

def get_dihedral(combined, conf, atomB_idx, atomC_idx):
    """
    Measures the dihedral (torsion) angle defined by four atoms: A-B-C-D, where B and C are specified,
    and A and D are their respective bonded neighbors.

    Args:
        combined (rdkit.Chem.Mol): RDKit molecule containing both protein and ligand.
        conf (rdkit.Chem.rdchem.Conformer): Conformer object of the combined molecule with 3D coordinates.
        atomB_idx (int): Index of atom B (typically the protein atom in the covalent bond).
        atomC_idx (int): Index of atom C (typically the ligand atom in the covalent bond).

    Returns:
        float: Dihedral angle in degrees.

    Raises:
        ValueError: If a neighboring atom (A or D) cannot be found to define the dihedral angle.
    """
    atom_B = combined.GetAtomWithIdx(atomB_idx)
    atom_C = combined.GetAtomWithIdx(atomC_idx)
    # get neighboring atoms
    atom_A_idx = get_neighbor(atom_B, atom_C.GetIdx())
    atom_D_idx = get_neighbor(atom_C, atom_B.GetIdx())

    # check
    if atom_A_idx is None or atom_D_idx is None:
        raise ValueError("Couldn't find bonded neighbor to define dihedral.")
    
    #dihedral
    return rdMolTransforms.GetDihedralDeg(conf, atom_A_idx, atomB_idx, atomC_idx, atom_D_idx)