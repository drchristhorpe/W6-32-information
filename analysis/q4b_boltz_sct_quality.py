"""Design question 4b -- is the Boltz A*02:01/PRAME SCT a sound design input?

Compares the aligned Boltz single-chain-trimer prediction to the experimental
HLA-A*02:01/PRAME structure at three scopes:
  * overall  -- HLA heavy + β2m (general quality),
  * peptide  -- the SLLQHLIGL peptide conformation,
  * interface-- the W6/32 footprint (α3 + β2m), reusing the B*27:05 footprint
                (α3 numbering is conserved and β2m is identical across class I).
Both structures are in the canonical frame (experimental pre-framed; Boltz
aligned by `align_structure`), so deviation is measured in-frame. Boltz's own
confidence (pLDDT/pTM) is read for context.

Run from the experiment root:  uv run python analysis/q4b_boltz_sct_quality.py
"""

from __future__ import annotations

import json
from pathlib import Path

from compare_structures import compare_structures

ROOT = Path(__file__).resolve().parent.parent
EXPERIMENTAL = ROOT / "structures" / "pdb" / "hla_a_02_01__SLLQHLIGL.pdb"
BOLTZ_DIR = ROOT / "sct_predictions" / "aligned" / "boltz" / "hla_a_02_01__single_chain_trimer__SLLQHLIGL"
BOLTZ = BOLTZ_DIR / "sample_0_predicted_structure_aligned.pdb"
BOLTZ_METRICS = ROOT / "sct_predictions" / "raw" / "boltz" / "hla_a_02_01__single_chain_trimer__SLLQHLIGL" / "outputs" / "files" / "prediction" / "metrics.json"
FOOTPRINT = ROOT / "interface_description" / "hla_b_27_05__FRYNGLIHR__w632__1_contacts.json"
OUT = ROOT / "analysis" / "results" / "q4b"

# scope -> (ref_chains, mobile_chains, footprint?)
SCOPES = {
    "overall_hla_b2m": (["A", "B"], ["A", "A"], None),
    "peptide":         (["P"],      ["A"],      None),
    "w632_interface":  (["A", "B"], ["A", "A"], FOOTPRINT),
}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for scope, (rc, mc, fp) in SCOPES.items():
        res = compare_structures(
            EXPERIMENTAL, BOLTZ, OUT,
            ref_chains=rc, mobile_chains=mc, footprint=fp, atoms="ca",
        )
        # name the output by scope (compare_structures names by stems, so copy intent here)
        tag = OUT / f"{scope}_compare.json"
        Path(res["_output"]).replace(tag)
        top = sorted(res["per_residue"], key=lambda e: e["rmsd"], reverse=True)[:5]
        print(
            f"{scope:18s} n={res['n_residues']:3d} globalRMSD={res['global_rmsd']:.3f} Å  top: "
            + ", ".join(f"{t['resname']}{t['resseq']}/{t['ref_chain']}={t['rmsd']:.2f}" for t in top)
        )

    m = json.load(open(BOLTZ_METRICS))["best_sample"]["metrics"]
    print(
        f"\nBoltz confidence: complex_plddt={m['complex_plddt']:.3f}  "
        f"ptm={m['ptm']:.3f}  structure_confidence={m['structure_confidence']:.3f}"
    )


if __name__ == "__main__":
    main()
