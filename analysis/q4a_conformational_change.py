"""Design question 4a -- how much does W6/32 change the B*27:05 conformation?

Compares the W6/32-bound B*27:05 complex to the apo B*27:05/RRFSRSPIRR structure
at the W6/32 interface footprint only (peptide ignored). Both experimental
structures are already in the same coordinate frame, so we measure in-frame
deviation (total change relative to the frame) and superposed deviation (internal
distortion only). Both crystallographic copies of the complex are used.

Run from the experiment root:  uv run python analysis/q4a_conformational_change.py
"""

from __future__ import annotations

import json
from pathlib import Path

from compare_structures import compare_structures

ROOT = Path(__file__).resolve().parent.parent
APO = ROOT / "structures" / "pdb" / "hla_b_27_05__RRFSRSPIRR.pdb"
OUT = ROOT / "analysis" / "results" / "q4a"

# (complex copy pdb, epitope chains, footprint->apo chain map)
# The apo reference uses chains A (HLA heavy) and B (b2m); copy 2 of the complex
# is chains D,E, so its footprint is remapped onto A,B.
COPIES = [
    ("hla_b_27_05__FRYNGLIHR__w632__1", ["A", "B"], {"A": "A", "B": "B"}),
    ("hla_b_27_05__FRYNGLIHR__w632__2", ["D", "E"], {"D": "A", "E": "B"}),
]


def _control() -> None:
    """Specificity control: is the change footprint-specific, or is the whole
    molecule offset? Compare apo vs bound copy 1 over the whole HLA heavy chain
    and split α1/α2 (the framing region) from α3 (the footprint region)."""
    import math

    bound = ROOT / "structures" / "pdb" / "hla_b_27_05__FRYNGLIHR__w632__1.pdb"
    res = compare_structures(APO, bound, OUT / "control",
                             ref_chains=["A"], mobile_chains=["A"], atoms="ca")

    def rms(rows):
        return math.sqrt(sum(e["rmsd"] ** 2 for e in rows) / len(rows)) if rows else 0.0

    a12 = [e for e in res["per_residue"] if e["resseq"] <= 180]
    a3 = [e for e in res["per_residue"] if e["resseq"] >= 182]
    print(
        f"CONTROL chain A in-frame: α1/α2 (1-180) RMSD={rms(a12):.3f} Å  "
        f"vs  α3 (≥182) RMSD={rms(a3):.3f} Å  "
        f"-> change is {'α3-specific' if rms(a3) > 2 * rms(a12) else 'not localised'}"
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for stem, epi, chain_map in COPIES:
        bound = ROOT / "structures" / "pdb" / f"{stem}.pdb"
        footprint = ROOT / "interface_description" / f"{stem}_contacts.json"
        # apo chains A (HLA) and B (b2m) map to the complex epitope chains
        for mode in ("inframe", "superposed"):
            res = compare_structures(
                APO, bound, OUT / mode,
                ref_chains=["A", "B"], mobile_chains=epi,
                footprint=footprint, footprint_chain_map=chain_map,
                atoms="ca", superpose=(mode == "superposed"),
            )
            top = sorted(res["per_residue"], key=lambda e: e["rmsd"], reverse=True)[:5]
            print(
                f"{stem} [{mode:10s}] n={res['n_residues']} "
                f"globalRMSD={res['global_rmsd']:.3f} Å  top: "
                + ", ".join(f"{t['resname']}{t['resseq']}/{t['ref_chain']}={t['rmsd']:.2f}" for t in top)
            )
    _control()


if __name__ == "__main__":
    main()
