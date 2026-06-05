"""Design question 4c (extension) -- best site for an α3 ↔ β2m-linker disulphide.

Scans residue pairs between the β2m→α1 linker (SCT res 124-143) and the α3 domain
(SCT res 324-419) of the A*02:01/PRAME single-chain trimer for geometry compatible
with engineering a disulphide (i.e. "if both were mutated to Cys"). Because the
linker is poly-Gly/Ser, a pseudo-Cβ is built from each residue's backbone so all
positions are scored consistently.

Engineering criteria (Disulfide-by-Design style):
  * Cβ-Cβ distance 3.0-4.5 Å (ideal ~3.85),
  * Cα-Cα distance 4.0-7.5 Å.
Candidates are ranked by closeness of Cβ-Cβ to the ideal.

Run from the experiment root:  uv run python analysis/q4c_disulphide_scan.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser

ROOT = Path(__file__).resolve().parent.parent
SCT = ROOT / "sct_predictions" / "aligned" / "boltz" / "hla_a_02_01__single_chain_trimer__SLLQHLIGL" / "sample_0_predicted_structure_aligned.pdb"
OUT = ROOT / "analysis" / "results" / "q4c"

# β2m → α1 linker. The native α1 domain starts G-S-HSMRYF (Gly1, Ser2, His3 at
# SCT res 142, 143, 144), so the Gly at 142 belongs to α1, not the linker: the
# true linker ends at 141.
LINKER = range(124, 142)   # β2m → α1 linker (124-141)
ALPHA3 = range(324, 420)   # α3 domain
IDEAL_CB = 3.85
CB_MIN, CB_MAX = 3.0, 4.5
CA_MIN, CA_MAX = 4.0, 7.5


def pseudo_cb(residue):
    """Idealised Cβ position from backbone N, CA, C (works for Gly too)."""
    try:
        n, ca, c = (residue["N"].coord, residue["CA"].coord, residue["C"].coord)
    except KeyError:
        return None, None
    b, cvec = ca - n, c - ca
    a = np.cross(b, cvec)
    cb = -0.58273431 * a + 0.56802827 * b - 0.54067466 * cvec + ca
    return ca, cb


def mean_plddt(residue):
    """Predicted local confidence (Boltz stores pLDDT in the B-factor column)."""
    vals = [a.bfactor for a in residue]
    return round(float(np.mean(vals)), 1) if vals else None


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    chain = PDBParser(QUIET=True).get_structure("sct", str(SCT))[0]["A"]
    res = {r.id[1]: r for r in chain if r.id[0] == " "}

    geom = {}
    for num, r in res.items():
        ca, cb = pseudo_cb(r)
        if cb is not None:
            geom[num] = (r.resname, ca, cb, mean_plddt(r))

    candidates = []
    for i in LINKER:
        if i not in geom:
            continue
        for j in ALPHA3:
            if j not in geom:
                continue
            d_cb = float(np.linalg.norm(geom[i][2] - geom[j][2]))
            d_ca = float(np.linalg.norm(geom[i][1] - geom[j][1]))
            if CB_MIN <= d_cb <= CB_MAX and CA_MIN <= d_ca <= CA_MAX:
                candidates.append({
                    "linker_res": i, "linker_resname": geom[i][0], "linker_plddt": geom[i][3],
                    "alpha3_res": j, "alpha3_resname": geom[j][0], "alpha3_plddt": geom[j][3],
                    "cb_cb": round(d_cb, 2), "ca_ca": round(d_ca, 2),
                    "score": round(abs(d_cb - IDEAL_CB), 3),
                })
    candidates.sort(key=lambda c: c["score"])

    linker_plddt = [geom[i][3] for i in LINKER if i in geom]
    mean_linker_plddt = round(sum(linker_plddt) / len(linker_plddt), 1) if linker_plddt else None

    result = {
        "structure": str(SCT),
        "linker_range": [LINKER.start, LINKER.stop - 1],
        "alpha3_range": [ALPHA3.start, ALPHA3.stop - 1],
        "criteria": {"cb_cb": [CB_MIN, CB_MAX], "ca_ca": [CA_MIN, CA_MAX], "ideal_cb": IDEAL_CB},
        "mean_linker_plddt": mean_linker_plddt,
        "n_candidates": len(candidates),
        "candidates": candidates,
    }
    with open(OUT / "disulphide_candidates.json", "w") as fh:
        json.dump(result, fh, indent=2)

    print(f"linker mean pLDDT = {mean_linker_plddt} (low = flexible; exact linker partner is unreliable)")
    print(f"{len(candidates)} candidate pair(s) (Cβ-Cβ {CB_MIN}-{CB_MAX} Å, Cα-Cα {CA_MIN}-{CA_MAX} Å):")
    for c in candidates[:10]:
        print(
            f"  linker {c['linker_resname']}{c['linker_res']} (pLDDT {c['linker_plddt']}) ↔ "
            f"α3 {c['alpha3_resname']}{c['alpha3_res']} (pLDDT {c['alpha3_plddt']})  "
            f"Cβ-Cβ={c['cb_cb']} Å  Cα-Cα={c['ca_ca']} Å"
        )


if __name__ == "__main__":
    main()
