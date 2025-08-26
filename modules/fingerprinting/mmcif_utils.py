from typing import Iterable, Optional, Set
from Bio.PDB import MMCIFParser, MMCIFIO, PDBIO, Select, Polypeptide
from rdkit.Chem import rdDetermineBonds
from rdkit.Chem import SDWriter
from rdkit.Chem import RWMol
from rdkit.Chem import Conformer
from rdkit.Chem import Atom
from rdkit.Geometry import Point3D
from rdkit.Chem.rdmolfiles import SDWriter
from rdkit.Chem import rdDetermineBonds
class ProteinChainSelect(Select):
    def __init__(self,
                 keep_chains: Optional[Set[str]] = None,
                 keep_altloc: str = "",
                 keep_h: bool = True):
        """
        keep_chains: if None -> keep all protein chains; else subset like {"A","B"}
        keep_altloc: which altloc to keep. "" keeps blank; "A" keeps altloc A, etc.
        keep_h: keep hydrogens (True) or drop them (False)
        """
        super().__init__()
        self.keep_chains = keep_chains
        self.keep_altloc = keep_altloc
        self.keep_h = keep_h

    def accept_chain(self, chain):
        if self.keep_chains is None:
            return 1
        return 1 if chain.id in self.keep_chains else 0

    def accept_residue(self, residue):
        # keep only amino-acid residues (drop waters/ligands/ions)
        return 1 if Polypeptide.is_aa(residue, standard=False) else 0

    def accept_atom(self, atom):
        # filter altlocs
        alt = atom.get_altloc() or ""
        if self.keep_altloc != "" and alt not in ("", self.keep_altloc):
            return 0
        # optional: drop hydrogens
        if not self.keep_h and (atom.element == "H" or atom.get_name().startswith("H")):
            return 0
        return 1

def extract_protein_receptor(
    cif_path: str,
    out_path: str,
    chain_ids: Optional[Iterable[str]] = None,
    model_index: int = 0,
    keep_altloc: str = "",          # "", or "A" to prefer altloc A, etc.
    drop_hydrogens: bool = False,
    out_format: str = "cif",        # "cif" or "pdb"
) -> int:
    """
    Read an mmCIF and write only protein chain(s) to a new file.

    Returns
    -------
    int
        Number of protein chains written.
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("X", cif_path)
    models = list(structure.get_models())
    if not models:
        raise ValueError("No models in structure.")
    if model_index >= len(models):
        raise IndexError(f"model_index {model_index} out of range (n_models={len(models)}).")
    model = models[model_index]

    # Determine which chains to keep
    protein_chain_ids = []
    for chain in model:
        any_protein = any(Polypeptide.is_aa(res, standard=False) for res in chain)
        if any_protein:
            protein_chain_ids.append(chain.id)

    if not protein_chain_ids:
        # nothing to write; create empty file
        open(out_path, "w").close()
        return 0

    if chain_ids is None:
        # Default: keep the first protein chain only (as “the receptor chain”)
        selected = {protein_chain_ids[0]}
    else:
        wanted = set(chain_ids)
        selected = {cid for cid in protein_chain_ids if cid in wanted}
        if not selected:
            raise ValueError(f"None of the requested chains {wanted} are protein chains. "
                             f"Available protein chains: {protein_chain_ids}")

    selector = ProteinChainSelect(
        keep_chains=selected,
        keep_altloc=keep_altloc,
        keep_h=not drop_hydrogens,
    )

    if out_format.lower() == "cif":
        io = MMCIFIO()
        io.set_structure(structure)
        io.save(out_path, select=selector)
    elif out_format.lower() == "pdb":
        io = PDBIO()
        io.set_structure(structure)
        io.save(out_path, select=selector)
    else:
        raise ValueError("out_format must be 'cif' or 'pdb'.")

    return len(selected)

def extract_ligands_to_sdf(
    cif_path: str,
    out_sdf: str,
    model_index: int = 0,
    exclude_resnames: Iterable[str] = ("HOH", "WAT", "DOD", "H2O"),
    exclude_common_ions: bool = True,
    min_atoms: int = 3,
) -> int:
    """
    Extract hetero small-molecule ligands from an mmCIF and save them to an SDF.

    Parameters
    ----------
    cif_path : str
        Path to input .cif (or .cif.gz).
    out_sdf : str
        Path to output .sdf (multi-molecule).
    model_index : int
        Which model to use (0 = first).
    exclude_resnames : Iterable[str]
        Residue names to ignore (waters by default).
    exclude_common_ions : bool
        If True, skip common monoatomic/diatomic ions (NA, CL, MG, CA, ZN, ...).
    min_atoms : int
        Skip residues with fewer than this many atoms (avoids single ions).

    Returns
    -------
    int
        Number of ligand molecules written.
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("X", cif_path)
    models = list(structure.get_models())
    if not models:
        raise ValueError("No models found in structure.")
    if model_index >= len(models):
        raise IndexError(f"Requested model_index {model_index} but only {len(models)} model(s) present.")
    model = models[model_index]

    exclude_resnames = set(exclude_resnames)
    common_ions: Set[str] = {
        "NA", "K", "CL", "CA", "MG", "MN", "ZN", "FE", "CU", "CO", "NI", "CD", "SR",
        "CS", "BR", "F", "I", "HG", "PB", "BA", "AL", "RB", "YB", "PT", "AG", "AU"
    }

    def is_hetero_residue(res) -> bool:
        hetflag = res.id[0]
        return hetflag.strip() != ""  # True for HETATM

    def is_water_or_excluded(res) -> bool:
        rn = res.get_resname().strip().upper()
        if rn in exclude_resnames:
            return True
        if exclude_common_ions and (rn in common_ions or len(rn) <= 2):
            return True
        return False

    def choose_altloc(atoms):
        chosen = {}
        for a in atoms:
            alt = a.get_altloc() or ""
            key = (a.get_id(), a.get_name())
            if key not in chosen:
                chosen[key] = a
            else:
                curr = chosen[key]
                curr_alt = curr.get_altloc() or ""
                if curr_alt != "" and (alt == "" or (curr_alt != "A" and alt == "A")):
                    chosen[key] = a
        return list(chosen.values())

    lig_residues = []
    for chain in model.get_chains():
        for res in chain.get_residues():
            if not is_hetero_residue(res):
                continue
            if is_water_or_excluded(res):
                continue
            if sum(1 for _ in res.get_atoms()) < min_atoms:
                continue
            lig_residues.append((chain.id, res))

    writer = SDWriter(out_sdf)
    n_written = 0

    for chain_id, res in lig_residues:
        atoms = list(res.get_atoms())
        atoms = choose_altloc(atoms)
        if len(atoms) < min_atoms:
            continue

        rw = RWMol()
        conf = Conformer(len(atoms))
        idx_map = {}
        for i, a in enumerate(atoms):
            elem = (a.element or a.get_name().strip()[0]).capitalize()
            if not elem.isalpha() or len(elem) > 2:
                elem = "C"
            try:
                rd_atom = Atom(elem)
            except Exception:
                rd_atom = Atom("C")
            rd_idx = rw.AddAtom(rd_atom)
            x, y, z = a.coord
            conf.SetAtomPosition(rd_idx, Point3D(float(x), float(y), float(z)))
            idx_map[i] = rd_idx

        rw.AddConformer(conf, assignId=True)
        mol = rw.GetMol()
        try:
            rdDetermineBonds.DetermineBonds(mol)
            Chem.SanitizeMol(mol)
        except Exception:
            try:
                rdDetermineBonds.DetermineConnectivity(mol)
                Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ADJUSTHS)
            except Exception:
                continue

        resname = res.get_resname().strip()
        resseq = res.id[1]
        icode = res.id[2].strip() if res.id[2] else ""
        title = f"{resname}_{chain_id}_{resseq}{icode}"
        mol.SetProp("_Name", title)
        mol.SetProp("resname", resname)
        mol.SetProp("chain_id", chain_id)
        mol.SetProp("resseq", str(resseq))
        mol.SetProp("icode", icode)

        writer.write(mol)
        n_written += 1

    writer.close()
    return n_written

def extract_ligand(
    cif_path: str,
    out_cif: str,
    model_index: int = 0,
    exclude_resnames: Iterable[str] = ("HOH", "WAT", "DOD", "H2O"),
    exclude_common_ions: bool = True,
    min_atoms: int = 3,
) -> int:
    """
    Extract small-molecule ligands from an mmCIF and write them to a new mmCIF file.

    Parameters
    ----------
    cif_path : str
        Path to input .cif (or .cif.gz).
    out_cif : str
        Path to output .cif (ligands only).
    model_index : int
        Which model to use (0 = first).
    exclude_resnames : Iterable[str]
        Residue names to ignore (waters by default).
    exclude_common_ions : bool
        If True, skip common monoatomic/diatomic ions (NA, CL, MG, CA, ZN, ...).
    min_atoms : int
        Skip residues with fewer than this many atoms (avoids single ions).

    Returns
    -------
    int
        Number of ligand residues written.
    """
    from Bio.PDB import MMCIFParser, MMCIFIO
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("X", cif_path)
    models = list(structure.get_models())
    if not models:
        raise ValueError("No models found in structure.")
    if model_index >= len(models):
        raise IndexError(f"Requested model_index {model_index} but only {len(models)} model(s) present.")
    model = models[model_index]

    exclude_resnames = set(exclude_resnames)
    common_ions: Set[str] = {
        "NA", "K", "CL", "CA", "MG", "MN", "ZN", "FE", "CU", "CO", "NI", "CD", "SR",
        "CS", "BR", "F", "I", "HG", "PB", "BA", "AL", "RB", "YB", "PT", "AG", "AU"
    }

    def is_hetero_residue(res) -> bool:
        hetflag = res.id[0]
        return hetflag.strip() != ""  # True for HETATM

    def is_water_or_excluded(res) -> bool:
        rn = res.get_resname().strip().upper()
        if rn in exclude_resnames:
            return True
        if exclude_common_ions and (rn in common_ions or len(rn) <= 2):
            return True
        return False

    ligand_chains = []
    n_ligands = 0
    for chain in model.get_chains():
        ligand_residues = []
        for res in chain.get_residues():
            if not is_hetero_residue(res):
                continue
            if is_water_or_excluded(res):
                continue
            if sum(1 for _ in res.get_atoms()) < min_atoms:
                continue
            ligand_residues.append(res)
        if ligand_residues:
            # Create a new chain with only ligand residues
            from Bio.PDB.Chain import Chain
            new_chain = Chain(chain.id)
            for res in ligand_residues:
                new_chain.add(res.copy())
            ligand_chains.append(new_chain)
            n_ligands += len(ligand_residues)

    if not ligand_chains:
        # Write an empty CIF file
        open(out_cif, "w").close()
        return 0

    # Build a new structure with only ligand chains
    from Bio.PDB.Model import Model
    from Bio.PDB.Structure import Structure
    new_structure = Structure("LIGANDS")
    new_model = Model(0)
    for ch in ligand_chains:
        new_model.add(ch)
    new_structure.add(new_model)

    io = MMCIFIO()
    io.set_structure(new_structure)
    io.save(out_cif)
    return n_ligands
