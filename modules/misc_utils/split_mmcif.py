#!/usr/bin/env python3


from rdkit.Chem import rdDetermineBonds

from pathlib import Path
import argparse
from typing import List, Tuple, Optional, Set, Dict

from Bio.PDB import MMCIFParser, MMCIFIO, StructureBuilder
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom

from rdkit import Chem
from rdkit.Chem import AllChem

AA3 = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
    "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
    "SEC","PYL"
}
NA3 = {"DA","DT","DG","DC","A","U","G","C","DI","DU","I"}  # simple list
WATER = {"HOH","WAT","H2O"}
IONS = {
    "NA","K","CL","CA","MG","ZN","MN","FE","CU","CO","NI","CD","SR","CS",
    "YB","RB","HG","PB","AG","AU","AL","GA","IN","LI","F","BR","I","BA"
}

def is_polymer_res(res: Residue) -> bool:
    # In mmCIF via Biopython, standard residues have hetflag " "
    return res.id[0] == " "

def chain_type(chain: Chain) -> str:
    aa = any((r.get_resname().strip().upper() in AA3) for r in chain if is_polymer_res(r))
    na = any((r.get_resname().strip().upper() in NA3) for r in chain if is_polymer_res(r))
    if aa and not na:
        return "protein"
    if na and not aa:
        return "nucleic"
    if aa and na:
        return "mixed"
    return "nonpoly"

def is_water(res: Residue) -> bool:
    return res.get_resname().strip().upper() in WATER

def is_ion(res: Residue) -> bool:
    rn = res.get_resname().strip().upper()
    return rn in IONS and len(res) == 1

def is_small_molecule(res: Residue) -> bool:
    if is_polymer_res(res):
        return False
    if is_water(res) or is_ion(res):
        return False
    return len(res) > 0

# --- RDKit helpers (distance-based bonding)
PT = Chem.GetPeriodicTable()
COV_RADII = {1:0.31,6:0.76,7:0.71,8:0.66,9:0.57,15:1.07,16:1.05,17:1.02,35:1.20,53:1.39}
METALS = {"Na","K","Mg","Ca","Mn","Fe","Co","Ni","Cu","Zn","Cd","Ag","Sr","Cs","Ba","Al","Ga","In","Au","Hg","Pb"}

def infer_element(atom_name: str, element_field: Optional[str]) -> str:
    if element_field and element_field.strip():
        return element_field.strip().capitalize()
    # PDB-style: element often right-justified in columns; fall back to letters in name
    name = atom_name.strip()
    # common cases: 'C1', 'CA', 'N', 'OXT', 'Cl', 'BR'
    letters = ''.join([c for c in name if c.isalpha()])
    if not letters:
        return "C"
    if letters[:2].capitalize() in ("Cl","Br"):
        return letters[:2].capitalize()
    return letters[0].upper()

def residue_to_rdkit(res):
    em = Chem.RWMol()
    conf = Chem.Conformer()
    for at in res:
        sym = infer_element(at.get_name(), getattr(at, "element", None))
        idx = em.AddAtom(Chem.Atom(sym))
        x, y, z = at.get_coord()
        conf.SetAtomPosition(idx, (float(x), float(y), float(z)))
    if em.GetNumAtoms() < min_atoms:
        return None
    em.AddConformer(conf)
    rdDetermineBonds.DetermineBonds(em)   # <- automatic bond inference
    try:
        Chem.SanitizeMol(em)
    except Exception:
        pass
    return em

def write_residue_sdf(res: Residue, out_path: Path) -> bool:
    mol = residue_to_rdkit(res)
    if mol is None:
        return False
    rn = res.get_resname().strip().upper()
    resseq = res.get_id()[1]
    icode = res.get_id()[2] if isinstance(res.get_id()[2], str) else ""
    title = f"{rn}_{resseq}{icode}"
    mol.SetProp("_Name", title)
    mol.SetProp("resname", rn)
    mol.SetProp("resnum", str(resseq))
    if icode:
        mol.SetProp("icode", icode)
    w = Chem.SDWriter(str(out_path))
    w.write(mol)
    w.close()
    return True

# --- Build a new Biopython Structure with a single chain (optionally keep solvent/ions)
def copy_chain_to_structure(src_struct: Structure, src_chain: Chain, keep_solvent=False, keep_ions=False) -> Structure:
    sb = StructureBuilder.StructureBuilder()
    sb.init_structure(src_struct.get_id() or "split")
    sb.init_model(0)
    sb.init_chain(src_chain.id)
    for res in src_chain:
        rn = res.get_resname().strip().upper()
        if not is_polymer_res(res):
            if is_water(res) and not keep_solvent:
                continue
            if is_ion(res) and not keep_ions:
                continue
            # by default skip non-polymers entirely
            if not (keep_solvent or keep_ions):
                continue
        het, resseq, icode = res.id
        sb.init_seg("    ")
        sb.init_residue(res.get_resname(), het, resseq, icode if isinstance(icode, str) else " ")
        for at in res:
            name = at.get_fullname().strip()
            coord = at.get_coord()
            bfac = float(at.get_bfactor())
            occ = float(at.get_occupancy() if at.get_occupancy() is not None else 1.0)
            altloc = at.get_altloc() if isinstance(at.get_altloc(), str) else " "
            element = infer_element(at.get_name(), getattr(at, "element", None))
            new_atom = Atom(name=name, coord=coord, bfactor=bfac, occupancy=occ,
                            altloc=altloc, fullname=name, serial_number=0, element=element)
            sb.structure[0][src_chain.id].add(new_atom)
    return sb.get_structure()

def split_mmcif(input_cif: str, outdir: str, include_na=False, keep_solvent=False, keep_ions=False) -> None:
    parser = MMCIFParser(QUIET=True)
    struct = parser.get_structure(Path(input_cif).stem, input_cif)
    out = Path(outdir)
    (out / "proteins").mkdir(parents=True, exist_ok=True)
    (out / "ligands_sdf").mkdir(parents=True, exist_ok=True)

    # Proteins: one file per protein chain
    for chain in struct[0]:
        ctype = chain_type(chain)
        if ctype == "protein" or (include_na and ctype == "nucleic"):
            sub = copy_chain_to_structure(struct, chain, keep_solvent=keep_solvent, keep_ions=keep_ions)
            io = MMCIFIO()
            io.set_structure(sub)
            io.save(str((out / "proteins" / f"{Path(input_cif).stem}_chain-{chain.id}.mmcif")))

    # Ligands: one SDF per small-molecule residue
    for chain in struct[0]:
        for res in chain:
            if is_small_molecule(res):
                rn = res.get_resname().strip().upper()
                resseq = res.get_id()[1]
                icode = res.get_id()[2] if isinstance(res.get_id()[2], str) else ""
                fname = f"{Path(input_cif).stem}_{rn}_{chain.id}{resseq}{icode}.sdf"
                write_residue_sdf(res, out / "ligands_sdf" / fname)

# --- CLI
def main():
    ap = argparse.ArgumentParser(description="Split a .mmcif into per-protein-chain mmCIFs and per-ligand SDFs using Biopython + RDKit")
    ap.add_argument("--input", required=True, help="Input .mmcif")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--include-na", action="store_true", help="Also export nucleic-acid chains as 'proteins'")
    ap.add_argument("--keep-solvent", action="store_true", help="Keep waters in protein outputs")
    ap.add_argument("--keep-ions", action="store_true", help="Keep simple ions in protein outputs")
    args = ap.parse_args()
    split_mmcif(args.input, args.outdir, include_na=args.include_na,
                keep_solvent=args.keep_solvent, keep_ions=args.keep_ions)

if __name__ == "__main__":
    main()
