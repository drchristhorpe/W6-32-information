"""compare_structures -- fine-scale difference between two structures.

Self-contained tool (own helpers + CLI) so it can be lifted into a standalone
library or skill file unchanged.

Reports global and per-residue deviation between a reference and a mobile
structure over a residue selection. Residues are paired by **local sequence
alignment** per chain, so it compares structures with different numbering or
chain layouts (e.g. a separate-chain experimental complex vs a single-chain
trimer prediction).

By default it measures deviation **in the input coordinate frame, without
re-superposing** -- the right behaviour when the structures have already been
brought into a common frame (by `align_structure`, or already co-framed). Pass
`--superpose` to first Kabsch-fit the selection (pure internal conformational
difference).

Coordinates and RMSD use BioPandas; sequence pairing uses BioPython.

Contract: takes two input filepaths (reference, mobile) and an output folder;
writes `<refstem>__vs__<mobstem>_compare.json` into the output folder.

Example (4a, footprint-restricted, same frame):
    python -m compare_structures \\
        structures/pdb/hla_b_27_05__RRFSRSPIRR.pdb \\
        structures/pdb/hla_b_27_05__FRYNGLIHR__w632__1.pdb out/ \\
        --ref-chains A,B --mobile-chains A,B \\
        --footprint interface_description/hla_b_27_05__FRYNGLIHR__w632__1_contacts.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.SeqUtils import seq1
from biopandas.pdb import PandasPdb


def _residue_table(atom_df, chain_id):
    """Ordered residues for one chain. Returns a list of dicts with the residue
    key fields and an atom_name -> (x,y,z) map, ordered as they appear."""
    sub = atom_df[atom_df.chain_id == chain_id]
    residues = []
    seen = {}
    for row in sub.itertuples(index=False):
        rid = (row.residue_number, getattr(row, "insertion", "") or "")
        entry = seen.get(rid)
        if entry is None:
            entry = {
                "chain": chain_id,
                "resseq": int(row.residue_number),
                "icode": (getattr(row, "insertion", "") or "").strip(),
                "resname": row.residue_name,
                "atoms": {},
            }
            seen[rid] = entry
            residues.append(entry)
        entry["atoms"][row.atom_name] = (row.x_coord, row.y_coord, row.z_coord)
    return residues


def _sequence(residues):
    return "".join(seq1(r["resname"], undef_code="X") for r in residues)


def _pair_residues(ref_res, mob_res):
    """Pair residues of two chains by local sequence alignment."""
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.mode = "local"
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    alignment = aligner.align(_sequence(ref_res), _sequence(mob_res))[0]
    pairs = []
    for (r0, r1), (m0, m1) in zip(*alignment.aligned):
        for ri, mi in zip(range(r0, r1), range(m0, m1)):
            pairs.append((ref_res[ri], mob_res[mi]))
    return pairs


def _kabsch(mobile_xyz, ref_xyz):
    """Rotation+translation mapping mobile onto ref (least-squares)."""
    mc, rc = mobile_xyz.mean(0), ref_xyz.mean(0)
    h = (mobile_xyz - mc).T @ (ref_xyz - rc)
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return rot, rc - rot @ mc


def compare_structures(
    reference_filepath: str | Path,
    input_filepath: str | Path,
    output_folder: str | Path,
    ref_chains: list[str] | None = None,
    mobile_chains: list[str] | None = None,
    footprint: str | Path | None = None,
    atoms: str = "ca",
    superpose: bool = False,
) -> dict:
    reference_filepath = Path(reference_filepath)
    input_filepath = Path(input_filepath)
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    ref_df = PandasPdb().read_pdb(str(reference_filepath)).df["ATOM"]
    mob_df = PandasPdb().read_pdb(str(input_filepath)).df["ATOM"]

    if ref_chains is None or mobile_chains is None:
        common = [c for c in ref_df.chain_id.unique() if c in set(mob_df.chain_id.unique())]
        ref_chains = mobile_chains = list(common)
    if len(ref_chains) != len(mobile_chains):
        raise ValueError("ref-chains and mobile-chains must have equal length.")

    keep = _footprint_filter(footprint)  # set of (chain, resseq) or None

    per_residue = []
    ref_pts, mob_pts = [], []  # paired atom coordinates for global RMSD / fit
    for rc, mc in zip(ref_chains, mobile_chains):
        ref_res = _residue_table(ref_df, rc)
        mob_res = _residue_table(mob_df, mc)
        if not ref_res or not mob_res:
            continue
        for ref_r, mob_r in _pair_residues(ref_res, mob_res):
            if keep is not None and (rc, ref_r["resseq"]) not in keep:
                continue
            names = ["CA"] if atoms == "ca" else sorted(
                set(ref_r["atoms"]) & set(mob_r["atoms"])
            )
            paired = [(ref_r["atoms"][n], mob_r["atoms"][n])
                      for n in names if n in ref_r["atoms"] and n in mob_r["atoms"]]
            if not paired:
                continue
            r_xyz = np.array([p[0] for p in paired], float)
            m_xyz = np.array([p[1] for p in paired], float)
            ref_pts.append(r_xyz)
            mob_pts.append(m_xyz)
            per_residue.append({
                "ref_chain": rc, "resseq": ref_r["resseq"], "resname": ref_r["resname"],
                "_r": r_xyz, "_m": m_xyz,
            })

    if not per_residue:
        raise ValueError("No residues matched the selection.")

    ref_all = np.vstack(ref_pts)
    mob_all = np.vstack(mob_pts)
    if superpose:
        rot, trans = _kabsch(mob_all, ref_all)
        mob_all = mob_all @ rot.T + trans
        # re-split transformed coords back per residue
        i = 0
        for entry in per_residue:
            k = len(entry["_m"])
            entry["_m"] = mob_all[i:i + k]
            i += k

    # finalise per-residue deviations
    table = []
    for entry in per_residue:
        d = np.linalg.norm(entry["_r"] - entry["_m"], axis=1)
        table.append({
            "ref_chain": entry["ref_chain"],
            "resseq": entry["resseq"],
            "resname": entry["resname"],
            "rmsd": round(float(np.sqrt((d ** 2).mean())), 3),
            "max_atom_dev": round(float(d.max()), 3),
            "n_atoms": int(len(d)),
        })

    global_rmsd = float(np.sqrt(((ref_all - mob_all) ** 2).sum(axis=1).mean()))
    table.sort(key=lambda e: (e["ref_chain"], e["resseq"]))
    result = {
        "reference": str(reference_filepath),
        "mobile": str(input_filepath),
        "atoms": atoms,
        "superposed": superpose,
        "footprint": str(footprint) if footprint else None,
        "n_residues": len(table),
        "n_atoms": int(len(ref_all)),
        "global_rmsd": round(global_rmsd, 3),
        "max_residue_rmsd": max(table, key=lambda e: e["rmsd"]),
        "per_residue": table,
    }
    out = output_folder / f"{reference_filepath.stem}__vs__{input_filepath.stem}_compare.json"
    with open(out, "w") as fh:
        json.dump(result, fh, indent=2)
    result["_output"] = str(out)
    return result


def _footprint_filter(footprint):
    if footprint is None:
        return None
    data = json.load(open(footprint))
    return {(r["chain"], r["resseq"]) for r in data["epitope_footprint"]}


def _chain_list(spec):
    return [c.strip() for c in spec.split(",") if c.strip()] if spec else None


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("reference_filepath")
    ap.add_argument("input_filepath", help="mobile structure to compare to the reference")
    ap.add_argument("output_folder")
    ap.add_argument("--ref-chains", default=None, help="comma list, parallel to --mobile-chains")
    ap.add_argument("--mobile-chains", default=None, help="comma list, parallel to --ref-chains")
    ap.add_argument("--footprint", default=None, help="contacts.json; restrict to its epitope_footprint residues")
    ap.add_argument("--atoms", choices=("ca", "heavy"), default="ca")
    ap.add_argument("--superpose", action="store_true", help="Kabsch-fit the selection before measuring")
    args = ap.parse_args(argv)

    result = compare_structures(
        args.reference_filepath, args.input_filepath, args.output_folder,
        ref_chains=_chain_list(args.ref_chains),
        mobile_chains=_chain_list(args.mobile_chains),
        footprint=args.footprint, atoms=args.atoms, superpose=args.superpose,
    )
    mr = result["max_residue_rmsd"]
    print(
        f"{Path(args.reference_filepath).stem} vs {Path(args.input_filepath).stem}: "
        f"{result['n_residues']} residues, global RMSD {result['global_rmsd']:.3f} Å "
        f"({'superposed' if args.superpose else 'in-frame'}); "
        f"max {mr['resname']}{mr['resseq']}/{mr['ref_chain']} {mr['rmsd']:.3f} Å "
        f"-> {result['_output']}"
    )


if __name__ == "__main__":
    main()
