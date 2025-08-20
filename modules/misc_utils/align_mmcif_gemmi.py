#!/usr/bin/env python3
import argparse, csv, math
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import multiprocessing as mp

import gemmi

AtomSel = Tuple[int, str, str]  # (resnum, icode, atom_name)

def is_poly_res(res: gemmi.Residue) -> bool:
    return res.is_polymer() and res.name not in ("HOH",)

def res_key(res: gemmi.Residue) -> Tuple[int, str]:
    return (res.seqid.num, res.seqid.icode or "")

def pick_model(st: gemmi.Structure) -> gemmi.Model:
    # use first model
    return st[0]

def pick_chain(model: gemmi.Model, chain_id: Optional[str]) -> Optional[gemmi.Chain]:
    if chain_id:
        for ch in model:
            if ch.name == chain_id:
                return ch
        return None
    # auto-pick first polymer chain
    for ch in model:
        if any(is_poly_res(res) for res in ch):
            return ch
    return None

def highest_occ_pos(atom: gemmi.Atom) -> gemmi.Position:
    # gemmi collapses altlocs into separate atoms; occupancy is on atom. Use as-is.
    return atom.pos

def get_fit_atoms(chain: gemmi.Chain, mode: str, residues_filter: Optional[set]=None) -> Dict[AtomSel, gemmi.Position]:
    """Collect candidate atom positions for fitting from a protein chain."""
    wanted = {"CA"} if mode == "ca" else {"N", "CA", "C"}
    out: Dict[AtomSel, gemmi.Position] = {}
    for res in chain:
        if not is_poly_res(res):
            continue
        rk = res_key(res)
        if residues_filter is not None and rk not in residues_filter:
            continue
        for at in res:
            if at.name in wanted:
                out[(rk[0], rk[1], at.name)] = highest_occ_pos(at)
    return out

def residues_within_radius_of_ligand(model: gemmi.Model,
                                     ref_chain: gemmi.Chain,
                                     lig_resname: str,
                                     radius: float) -> set:
    """
    Build a set of protein residue keys within radius Å of ANY atom in any residue with name lig_resname (in reference model).
    """
    # collect ligand atoms (filter by residue name, exclude polymer)
    lig_positions: List[gemmi.Position] = []
    for ch in model:
        for res in ch:
            if res.name == lig_resname and not res.is_polymer():
                for at in res:
                    lig_positions.append(at.pos)
    if not lig_positions:
        return set()

    near_keys = set()
    for res in ref_chain:
        if not is_poly_res(res):
            continue
        # check any atom in residue within radius
        take = False
        for at in res:
            p = at.pos
            for lp in lig_positions:
                if p.dist(lp) <= radius:
                    take = True
                    break
            if take:
                break
        if take:
            near_keys.add(res_key(res))
    return near_keys

def superpose(ref_pts: List[gemmi.Position], mob_pts: List[gemmi.Position]) -> Tuple[gemmi.Transform, float]:
    # gemmi.superpose_positions returns (transform, rmsd)
    tf, rmsd = gemmi.superpose_positions(mob_pts, ref_pts)  # move mob onto ref
    return tf, float(rmsd)

def apply_transform_to_model(model: gemmi.Model, tf: gemmi.Transform) -> None:
    for ch in model:
        for res in ch:
            for at in res:
                at.pos = tf.apply(at.pos)

def write_mmCIF(struct: gemmi.Structure, out_path: Path) -> None:
    # Write aligned structure as mmCIF
    doc = gemmi.cif.Document()
    block = gemmi.make_mmcif_block(struct, gemmi.MmcifOutputOptions())
    doc.add_block(block)
    doc.write_file(str(out_path))

def write_pdb(struct: gemmi.Structure, out_path: Path) -> None:
    pdb = struct.make_pdb_string()
    out_path.write_text(pdb)

def build_point_pairs(ref_chain: gemmi.Chain,
                      mob_chain: gemmi.Chain,
                      mode: str,
                      region_keys: Optional[set]) -> Tuple[List[gemmi.Position], List[gemmi.Position]]:
    ref_atoms = get_fit_atoms(ref_chain, mode, region_keys)
    mob_atoms = get_fit_atoms(mob_chain, mode, region_keys)
    # match by (resnum, icode, atom_name)
    keys = sorted(set(ref_atoms.keys()).intersection(mob_atoms.keys()))
    ref_pts = [ref_atoms[k] for k in keys]
    mob_pts = [mob_atoms[k] for k in keys]
    return ref_pts, mob_pts

def process_one(args) -> Tuple[str, str, str, float, int, str]:
    (ref_path, mob_path, out_path, ref_chain_id, mob_chain_id,
     fit_mode, pocket_ligand, pocket_radius, overwrite) = args
    try:
        ref = gemmi.read_structure(str(ref_path))
        mob = gemmi.read_structure(str(mob_path))
        ref_model = pick_model(ref)
        mob_model = pick_model(mob)

        ref_chain = pick_chain(ref_model, ref_chain_id)
        if ref_chain is None:
            raise ValueError(f"Reference chain '{ref_chain_id}' not found (or no polymer chain).")
        mob_chain = pick_chain(mob_model, mob_chain_id)
        if mob_chain is None:
            # heuristic: pick chain with most matching residues
            candidates = [ch for ch in mob_model if any(is_poly_res(r) for r in ch)]
            best, best_n = None, -1
            for ch in candidates:
                rkeys = set(get_fit_atoms(ref_chain, fit_mode).keys())
                mkeys = set(get_fit_atoms(ch, fit_mode).keys())
                n = len(rkeys & mkeys)
                if n > best_n:
                    best, best_n = ch, n
            mob_chain = best
        if mob_chain is None:
            raise ValueError("No suitable mobile chain found.")

        region_keys = None
        if pocket_ligand:
            region_keys = residues_within_radius_of_ligand(ref_model, ref_chain, pocket_ligand, pocket_radius)
            if not region_keys:
                raise ValueError(f"No residues found within {pocket_radius} Å of ligand '{pocket_ligand}' in reference.")

        ref_pts, mob_pts = build_point_pairs(ref_chain, mob_chain, fit_mode, region_keys)
        if len(ref_pts) < (3 if fit_mode == "ca" else 3):
            raise ValueError(f"Not enough common points for fit (found {len(ref_pts)}).")

        tf, rmsd = superpose(ref_pts, mob_pts)
        apply_transform_to_model(mob_model, tf)

        out_path.parent.mkdir(parents=True, exist_ok=True)
        if out_path.exists() and not overwrite:
            raise FileExistsError(f"{out_path} exists. Use --overwrite.")

        # Prefer mmCIF; fall back to PDB if needed
        try:
            write_mmCIF(mob, out_path)
            fmt = "mmcif"
        except Exception:
            write_pdb(mob, out_path.with_suffix(".pdb"))
            fmt = "pdb"

        return (mob_path.name, ref_chain.name, mob_chain.name, rmsd, len(ref_pts), fmt)
    except Exception as e:
        return ("ERROR", str(ref_path), f"{mob_path.name}: {e}", float("nan"), 0, "")

def pdb_to_mmcif(pdb_path: str, cif_path: str = None) -> str:
    """
    Convert a PDB file to mmCIF using gemmi. Returns the output path.
    Preserves chain IDs, altlocs, occupancies, anisotropic B-factors, etc.
    """
    pdb_path = Path(pdb_path)
    if cif_path is None:
        cif_path = str(pdb_path.with_suffix(".mmcif"))

    # read PDB (gz also supported)
    st = gemmi.read_structure(str(pdb_path))
    st.remove_alternative_conformations()  # optional: keep highest-occ; remove this line to keep all alts

    # write mmCIF
    doc = gemmi.cif.Document()
    block = gemmi.make_mmcif_block(st, gemmi.MmcifOutputOptions())
    doc.add_block(block)
    doc.write_file(cif_path)
    return cif_path
    
def main():
    ap = argparse.ArgumentParser(description="Fast batch alignment of mmCIFs to a reference (ligands move but are not used for fitting).")
    ap.add_argument("--ref", required=True, help="Reference .mmcif file")
    ap.add_argument("--inputs", required=True, help="Directory with .mmcif/.cif files to align")
    ap.add_argument("--outdir", required=True, help="Directory for aligned outputs")
    ap.add_argument("--ref-chain", default=None, help="Reference chain ID (default: auto first polymer)")
    ap.add_argument("--mob-chain", default=None, help="Mobile chain ID (default: auto best match)")
    ap.add_argument("--fit", choices=["ca", "backbone"], default="ca", help="Fit using Cα or backbone atoms")
    ap.add_argument("--pocket-ligand", default=None, help="Ligand residue name in REFERENCE (e.g., HEM, NAD, LIG) to define pocket for region-limited fit")
    ap.add_argument("--pocket-radius", type=float, default=6.0, help="Radius (Å) around the ligand for region-limited fit")
    ap.add_argument("--nproc", type=int, default=max(1, mp.cpu_count()//2), help="Parallel workers")
    ap.add_argument("--overwrite", action="store_true", help="Overwrite outputs if present")
    ap.add_argument("--logcsv", default="alignment_log.csv", help="Output CSV log")
    args = ap.parse_args()

    ref_path = Path(args.ref)
    # If reference is a PDB (or similar), convert it to mmCIF in the same location
    ref_str = str(ref_path).lower()
    if ref_str.endswith((".pdb", ".ent", ".pdb.gz")) and not ref_str.endswith((".cif", ".mmcif")):
        # generate mmCIF from PDB; pdb_to_mmcif returns the output path (string)
        try:
            cif_out = pdb_to_mmcif(str(ref_path))
            ref_path = Path(cif_out)
        except Exception as e:
            raise SystemExit(f"Failed to convert reference PDB to mmCIF: {e}")

    in_dir = Path(args.inputs)
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    inputs = sorted([*in_dir.glob("*.mmcif"), *in_dir.glob("*.cif")])
    if not inputs:
        raise SystemExit("No mmCIF/CIF files found.")

    job_args = []
    for p in inputs:
        out = out_dir / p.with_suffix(".mmcif").name
        job_args.append((
            ref_path, p, out,
            args.ref_chain, args.mob_chain,
            args.fit, args.pocket_ligand, args.pocket_radius,
            args.overwrite
        ))

    with mp.Pool(processes=args.nproc) as pool:
        results = list(pool.map(process_one, job_args))

    # write CSV log
    with open(args.logcsv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["file", "ref_chain", "mob_chain", "rmsd", "n_points", "format"])
        for r in results:
            w.writerow(r)

    ok = [r for r in results if r[0] != "ERROR"]
    fail = [r for r in results if r[0] == "ERROR"]
    print(f"Aligned {len(ok)} / {len(results)} files.")
    if ok:
        mean_rmsd = sum(r[3] for r in ok) / len(ok)
        print(f"Mean RMSD: {mean_rmsd:.3f} Å using {ok[0][4]}± (varies) points per structure.")
    if fail:
        print(f"{len(fail)} failures. See {args.logcsv} for details.")
        
if __name__ == "__main__":
    main()
