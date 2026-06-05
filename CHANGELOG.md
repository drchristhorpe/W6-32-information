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

### Aligned all SCT predictions

- `analysis/align_all_predictions.py` (driver; imports the tool as a library)
  aligned all **54** predictions (14 AlphaFold3, 14 Boltz, 13 ESMFold2, 13
  ESMFold2-fast) onto the canonical frame → `sct_predictions/aligned/<predictor>/<sct_id>/`
  (gitignored, with the data).
- Every prediction matched the full 180 α1/α2 residues; fit RMSD 0.36–0.71 Å.
  A\*02:01 SCTs fit tightest (same allele as the 1HHK reference); B\*27 / A\*11 /
  A\*30 slightly higher, consistent with cross-allele superposition. The
  A\*02:01/PRAME Boltz SCT (focus of 4b) is the single tightest at 0.362 Å.

### Tool 2 (`interface-contacts`) + W6/32 footprint

- **Built `tools/interface_contacts/`** (+ skill `interface-contacts`) — heavy-atom
  contacts between an epitope and a paratope chain group via BioPython
  `NeighborSearch`; reports per-side residue footprints + atom pairs. Distance
  only (no H-bond/π typing). Cutoff is a configurable module constant
  `DEFAULT_CUTOFF = 5.0` Å (widened from 4.0 for the design phase, per request).
- **Generated W6/32 footprints** for both crystallographic copies of the complex
  → `interface_description/` (tracked). Heavy chain = paratope, HLA heavy + β2m =
  epitope; light chain excluded (SCT steric clash).
- **Validation:** the footprint is the canonical conformational W6/32 epitope —
  HLA α3 domain (~222–243) + β2m (~1–6, 58–61), with a Glu229–Arg69 salt bridge
  (~2.8 Å). Reproducible across both copies (22 vs 23 epitope residues). Recorded
  in [CONCLUSIONS.md](CONCLUSIONS.md) (Aim 1).
- Added **CONCLUSIONS.md** as the reader-facing answers document for the design
  questions (separate from this lab notebook).

Next: build Tool 3 (`compare_structures`), then answer 4a/4b/4c.
