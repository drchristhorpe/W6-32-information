"""Align every SCT prediction onto the canonical frame.

Orchestration driver (not a tool): it imports the `align_structure` tool as a
library and runs it over all four predictors, writing each result into
`sct_predictions/aligned/<predictor>/<sct_id>/`. Per-SCT subfolders keep the
non-identifying Boltz/ESMFold output filenames from colliding.

Run from the experiment root:  uv run python analysis/align_all_predictions.py
"""

from __future__ import annotations

from pathlib import Path

from align_structure import align_structure

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "sct_predictions" / "raw"
OUT = ROOT / "sct_predictions" / "aligned"
REFERENCE = ROOT / "structures" / "1hhk_1_aligned.pdb"


def _model_paths():
    """Yield (predictor, sct_id, model_file) for every prediction."""
    for d in sorted((RAW / "alphafold3").glob("fold_*")):
        hits = list(d.glob("*_model_0.cif"))
        if hits:
            yield "alphafold3", d.name, hits[0]

    for d in sorted((RAW / "boltz").iterdir()):
        model = d / "outputs" / "files" / "prediction" / "sample_0_predicted_structure.cif"
        if model.exists():
            yield "boltz", d.name, model

    for predictor in ("esmfold2", "esmfold2-fast"):
        for d in sorted((RAW / predictor).iterdir()):
            if d.is_dir():
                model = d / f"{d.name}.pdb"
                if model.exists():
                    yield predictor, d.name, model


def main() -> None:
    n = 0
    for predictor, sct_id, model in _model_paths():
        out_dir = OUT / predictor / sct_id
        result = align_structure(model, out_dir, reference=REFERENCE)
        n += 1
        print(
            f"[{predictor:13s}] {sct_id:48s} "
            f"{result['n_matched_residues']:3d} res  RMSD {result['rmsd']:.3f} Å"
        )
    print(f"\nAligned {n} predictions -> {OUT}")


if __name__ == "__main__":
    main()
