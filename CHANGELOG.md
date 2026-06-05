# Changelog / lab notebook

Reverse-chronological record of work on the W6/32 information experiment. Each
entry: what was done, why, and what it produced. Paired with [PLAN.md](PLAN.md).

## 2026-06-05

### Project scaffold + Tool 1 (`align-structure`)

- **Scaffolded** a `uv` project (`w632-tools`, Python 3.14) — a member of the
  existing home-level uv workspace. Smoke-tested that BioPython 1.87 and
  BioPandas 0.5.1 parse the canonical reference on 3.14.
  - Note: `uv sync` reshapes the *shared* workspace `.venv` to the active
    project's deps, so it temporarily uninstalled sibling members' heavy
    packages (torch/transformers/etc.). These are restored automatically on the
    next `uv run`/`uv sync` in those projects (the workspace lockfile still
    covers them).
- **Tooling layout decision:** each tool lives in its own self-contained folder
  under `tools/<name>/` (own `pyproject.toml` + `README.md` + module), so it can
  be copied out and used elsewhere (`uv run`/`pip install`). Tools are wired into
  this repo's shared env as editable path dependencies of `w632-tools`.
- **Skill per tool:** each tool also ships a matching Claude skill at
  `.claude/skills/<name>/SKILL.md` that wraps the tool's CLI (single source of
  truth; no duplicated logic). Added `align-structure`.
- **Data-in-git decision:** `sct_predictions/` (680 MB, mostly AF3 MSAs) is kept
  out of git (`.gitignore`); only code, docs, and `structures/` are versioned.
- **Built `tools/align_structure/`** — superpose a structure onto the canonical
  frame ([structures/1hhk_1_aligned.pdb](structures/1hhk_1_aligned.pdb)), fitting
  on chain A α1/α2 Cα (residues 1–180).
  - Fit residues are matched by **local sequence alignment** (BLOSUM62,
    Smith-Waterman), not by residue numbering — so one tool handles both the
    separate-chain experimental complexes and the single-chain SCT predictions
    (HLA α fused mid-sequence with independent numbering).
  - Input auto-detects `.pdb`/`.cif` (ESMFold PDB vs Boltz/AF3 CIF); output is
    always PDB + an `<stem>_align.json` (rotation, translation, matched count,
    fit RMSD).
  - **Verified:** experimental A\*02:01/PRAME → 180 residues, 0.682 Å (already
    near the canonical frame, as expected); Boltz A\*02:01/PRAME SCT single chain
    → 180 residues, 0.362 Å. Standalone isolated install reproduces the result.

Next: align all SCT predictions into `sct_predictions/aligned/`, then build
Tool 2 (`interface_contacts`).
