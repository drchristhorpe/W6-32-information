"""Design question 4d (extension) -- which predictor is the best design start?

Ranks the A*02:01/PRAME single-chain-trimer predictions (AlphaFold3, Boltz2,
ESMFold2, ESMFold2-fast) by how close their W6/32 interface region is to the
*bound* structure. No W6/32-bound A*02:01 exists, so the bound reference is the
W6/32-bound B*27:05 complex (chains A,B); all predictors are scored against the
same reference, so the small (conserved, ~1.4 Å) allele offset is constant and the
ranking is fair. In-frame (all SCTs aligned to canonical), restricted to the W6/32
footprint; α3-only and β2m-only sub-RMSDs reported too.

Run from the experiment root:  uv run python analysis/q4d_predictor_interface_to_bound.py
"""

from __future__ import annotations

import math
from pathlib import Path

from compare_structures import compare_structures

ROOT = Path(__file__).resolve().parent.parent
BOUND = ROOT / "structures" / "pdb" / "hla_b_27_05__FRYNGLIHR__w632__1.pdb"
FOOTPRINT = ROOT / "interface_description" / "hla_b_27_05__FRYNGLIHR__w632__1_contacts.json"
OUT = ROOT / "analysis" / "results" / "q4d"
AL = ROOT / "sct_predictions" / "aligned"

PREDICTIONS = {
    "alphafold3":    AL / "alphafold3/fold_hla_a_02_01_single_chain_trimer_sllqhligl/fold_hla_a_02_01_single_chain_trimer_sllqhligl_model_0_aligned.pdb",
    "boltz":         AL / "boltz/hla_a_02_01__single_chain_trimer__SLLQHLIGL/sample_0_predicted_structure_aligned.pdb",
    "esmfold2":      AL / "esmfold2/hla_a_02_01__single_chain_trimer__SLLQHLIGL/hla_a_02_01__single_chain_trimer__SLLQHLIGL_aligned.pdb",
    "esmfold2-fast": AL / "esmfold2-fast/hla_a_02_01__single_chain_trimer__SLLQHLIGL/hla_a_02_01__single_chain_trimer__SLLQHLIGL_aligned.pdb",
}


def _rms(rows):
    return math.sqrt(sum(e["rmsd"] ** 2 for e in rows) / len(rows)) if rows else float("nan")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    table = []
    for name, pdb in PREDICTIONS.items():
        res = compare_structures(
            BOUND, pdb, OUT / name,
            ref_chains=["A", "B"], mobile_chains=["A", "A"],
            footprint=FOOTPRINT, atoms="ca",
        )
        # chain A = HLA heavy (α3 footprint), chain B = β2m
        a3 = [e for e in res["per_residue"] if e["ref_chain"] == "A"]
        b2m = [e for e in res["per_residue"] if e["ref_chain"] == "B"]
        table.append((name, res["n_residues"], res["global_rmsd"], _rms(a3), _rms(b2m)))

    table.sort(key=lambda r: r[2])
    print(f"{'predictor':14s} {'n':>3s} {'footprint':>10s} {'α3':>7s} {'β2m':>7s}   (RMSD to W6/32-bound B*27:05, Å)")
    for name, n, g, a3, b2m in table:
        print(f"{name:14s} {n:>3d} {g:>10.3f} {a3:>7.3f} {b2m:>7.3f}")
    print(f"\nBest interface match to bound: {table[0][0]} ({table[0][2]:.3f} Å footprint RMSD)")


if __name__ == "__main__":
    main()
