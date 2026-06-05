"""align_structure -- superpose a structure onto a canonical reference frame.

Self-contained tool (own helpers + CLI) so it can be lifted into a standalone
library or skill file unchanged.

The fitting selection is matched by *local sequence alignment* of the reference
selection against the mobile chain, rather than by assuming residue numbers line
up. This lets the same tool align both the separate-chain experimental
structures and the single-chain single-chain-trimer (SCT) predictions (whose
HLA alpha chain is fused mid-sequence with its own numbering).

Contract: takes an input filepath and an output folder; writes the aligned
structure plus an `align.json` describing the transform into the output folder.

Example:
    python -m w632_tools.align_structure \\
        structures/pdb/hla_a_02_01__SLLQHLIGL.pdb out/ \\
        --reference structures/1hhk_1_aligned.pdb \\
        --ref-chain A --ref-residues 1-180 --mobile-chain A
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Superimposer
from Bio.SeqUtils import seq1

CA = "CA"


def _load_structure(path: str | Path, name: str):
    """Parse a structure, auto-selecting the parser from the file extension so
    the tool ingests both .pdb (e.g. ESMFold) and .cif (e.g. Boltz/AF3)."""
    path = Path(path)
    if path.suffix.lower() in (".cif", ".mmcif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure(name, str(path))


def _standard_ca_residues(chain, residue_range=None):
    """Residues in `chain` that have a CA atom, optionally within an inclusive
    (lo, hi) residue-number range. Returns a list of Residue objects in order."""
    out = []
    for res in chain:
        if res.id[0] != " ":  # skip hetero/water
            continue
        if CA not in res:
            continue
        if residue_range is not None:
            lo, hi = residue_range
            if not (lo <= res.id[1] <= hi):
                continue
        out.append(res)
    return out


def _sequence(residues):
    return "".join(seq1(r.resname, custom_map={}, undef_code="X") for r in residues)


def _match_by_local_alignment(ref_residues, mob_residues):
    """Pair reference residues with mobile residues via local (Smith-Waterman)
    sequence alignment. Returns matched (ref_residue, mob_residue) pairs."""
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.mode = "local"
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1

    alignment = aligner.align(_sequence(ref_residues), _sequence(mob_residues))[0]
    ref_blocks, mob_blocks = alignment.aligned
    pairs = []
    for (r0, r1), (m0, m1) in zip(ref_blocks, mob_blocks):
        for ri, mi in zip(range(r0, r1), range(m0, m1)):
            pairs.append((ref_residues[ri], mob_residues[mi]))
    return pairs, alignment.score


def align_structure(
    input_filepath: str | Path,
    output_folder: str | Path,
    reference: str | Path,
    ref_chain: str = "A",
    ref_residues: tuple[int, int] | None = (1, 180),
    mobile_chain: str = "A",
    mobile_residues: tuple[int, int] | None = None,
) -> dict:
    """Superpose `input_filepath` onto `reference` over the matched CA atoms and
    write the transformed structure + transform metadata into `output_folder`."""
    input_filepath = Path(input_filepath)
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    ref_struct = _load_structure(reference, "ref")
    mob_struct = _load_structure(input_filepath, "mob")

    ref_res = _standard_ca_residues(ref_struct[0][ref_chain], ref_residues)
    mob_res = _standard_ca_residues(mob_struct[0][mobile_chain], mobile_residues)
    if not ref_res or not mob_res:
        raise ValueError("No CA residues found for the requested selection.")

    pairs, score = _match_by_local_alignment(ref_res, mob_res)
    if len(pairs) < 3:
        raise ValueError(f"Too few matched residues ({len(pairs)}) to superpose.")

    ref_atoms = [r[CA] for r, _ in pairs]
    mob_atoms = [m[CA] for _, m in pairs]

    sup = Superimposer()
    sup.set_atoms(ref_atoms, mob_atoms)
    sup.apply(mob_struct.get_atoms())  # transform the whole mobile structure

    io = PDBIO()
    io.set_structure(mob_struct)
    aligned_path = output_folder / f"{input_filepath.stem}_aligned.pdb"
    io.save(str(aligned_path))

    rot, trans = sup.rotran
    result = {
        "input": str(input_filepath),
        "reference": str(reference),
        "aligned_output": str(aligned_path),
        "ref_selection": {"chain": ref_chain, "residues": ref_residues},
        "mobile_selection": {"chain": mobile_chain, "residues": mobile_residues},
        "n_matched_residues": len(pairs),
        "alignment_score": float(score),
        "rmsd": float(sup.rms),
        "rotation": rot.tolist(),
        "translation": trans.tolist(),
    }
    with open(output_folder / f"{input_filepath.stem}_align.json", "w") as fh:
        json.dump(result, fh, indent=2)
    return result


def _residue_range(spec: str | None) -> tuple[int, int] | None:
    if spec is None:
        return None
    lo, hi = spec.split("-")
    return int(lo), int(hi)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("input_filepath", help="PDB to align (mobile)")
    ap.add_argument("output_folder", help="folder for aligned PDB + align.json")
    ap.add_argument(
        "--reference",
        default="structures/1hhk_1_aligned.pdb",
        help="canonical reference PDB (default: %(default)s)",
    )
    ap.add_argument("--ref-chain", default="A")
    ap.add_argument(
        "--ref-residues",
        default="1-180",
        help="inclusive residue range for the fit, e.g. 1-180 (default: %(default)s)",
    )
    ap.add_argument("--mobile-chain", default="A")
    ap.add_argument(
        "--mobile-residues",
        default=None,
        help="optional residue range to restrict the mobile selection",
    )
    args = ap.parse_args(argv)

    result = align_structure(
        args.input_filepath,
        args.output_folder,
        reference=args.reference,
        ref_chain=args.ref_chain,
        ref_residues=_residue_range(args.ref_residues),
        mobile_chain=args.mobile_chain,
        mobile_residues=_residue_range(args.mobile_residues),
    )
    print(
        f"aligned {Path(args.input_filepath).name}: "
        f"{result['n_matched_residues']} residues, RMSD {result['rmsd']:.3f} Å "
        f"-> {result['aligned_output']}"
    )


if __name__ == "__main__":
    main()
