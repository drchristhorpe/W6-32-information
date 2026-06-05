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

### Tool 3 (`compare-structures`)

- **Built `tools/compare_structures/`** (+ skill `compare-structures`) — global +
  per-residue deviation between two structures over a residue selection. Residues
  paired by **local sequence alignment** per chain (handles experimental↔SCT
  numbering/chain-layout differences). Measures **in the input frame by default**
  (inputs assumed pre-aligned), with `--superpose` for a Kabsch fit. Optional
  `--footprint` restricts to an `interface-contacts` footprint. BioPandas for
  coordinates/RMSD, BioPython for pairing.
- **Verified:** self-comparison = 0.000 Å; 4a validation run (apo vs W6/32-bound
  B\*27:05 at the footprint) gave 2.22 Å in-frame / 1.23 Å superposed — sensible.

All three composable tools (+ skills) are now complete.

### Design-question analyses (Aims 1–4; 4a/4b/4c)

Drivers in `analysis/` (import the tools as libraries); result JSONs in
`analysis/results/`. Full write-ups in [CONCLUSIONS.md](CONCLUSIONS.md).

- Added `compare_structures --footprint-chain-map` (remap footprint chains onto
  the reference, e.g. copy-2 D,E → A,B).
- **4a** (`q4a_conformational_change.py`): W6/32 binding reorients the α3 domain
  (footprint 2.2 Å in-frame / 1.2 Å superposed, both copies), while a specificity
  control shows α1/α2 unchanged (0.48 Å) vs α3 (2.78 Å) — change is α3-specific.
- **Aim 2:** W6/32 epitope conserved A\*02:01 vs B\*27:05 (1.36 Å at the footprint).
- **4b** (`q4b_boltz_sct_quality.py`): Boltz SCT is design-grade overall (1.01 Å)
  and for the PRAME peptide (1.00 Å), with high confidence (pLDDT 0.90/pTM 0.93);
  W6/32 epitope mostly good except the plastic α3 222–228 loop (~4 Å) — the same
  loop implicated in 4a.
- **4c** (`q4c_disulphide_scan.py`): geometric Cβ–Cβ scan of β2m→α1 linker (124–141,
  after correcting the α1 G-S-HSMRYF boundary) × α3. Top disulphide sites
  Gly135↔Leu407 and Gly132↔Lys327 (α3 anchors confidently placed, outside the
  W6/32 footprint); linker pLDDT ~36 means the linker partner needs explicit
  modelling.

Experiment complete: three composable tools + skills, all four aims and the
extension answered.

### Tool 4 (`suggest-disulphides`) + rigorous 4c

- **Built `tools/suggest_disulphides/`** (+ skill) — proposes engineerable
  disulphides between residue groups by modelling Sγ over three χ1 rotamers and
  requiring real geometry (Sγ–Sγ ≈ 2.03 Å, χ3 ≈ ±90°, Cβ–Sγ–Sγ ≈ 105°), not just
  Cβ–Cβ distance. `--exclude-footprint`, pLDDT-aware. BioPython + NumPy. Kept the
  standalone Cβ–Cβ scan (`analysis/q4c_*`) alongside, per request.
- **Validated:** scanning the whole SCT chain, the tool recovers all three native
  disulphides (α2 Cys242–Cys305, α3 Cys344–Cys400, β2m Cys49–Cys104) with near-ideal
  geometry and no false positives.
- **Re-ran 4c rigorously:** the 6 Cβ–Cβ proximity hits are all rejected by proper
  Sγ/χ3 geometry → **no viable linker×α3 disulphide in the static Boltz model**
  (linker pLDDT ~36). Updated CONCLUSIONS.md 4c: α3 anchor (Leu407/Lys327) is firm,
  but the disordered linker partner must be sited by explicit modelling/MD.

### 4d — best predictor for design (interface proximity to bound)

- `analysis/q4d_predictor_interface_to_bound.py`: ranked all four predictors'
  A\*02:01/PRAME SCTs by W6/32-footprint RMSD to the W6/32-bound B\*27:05 complex
  (common reference → fair ranking despite the constant allele offset).
- **Result:** ESMFold2-fast and Boltz tie on the overall footprint (1.40 Å);
  Boltz is best at the design-critical α3 loop (1.61 Å). AlphaFold3 and standard
  ESMFold2 are worse at α3 (1.92 / 2.04 Å) despite AF3 modelling β2m best. All
  still 1.6–2.0 Å from the bound α3, so the loop needs refinement regardless.
- **Conclusion:** Boltz-2 / ESMFold2-fast are the best design starting points,
  Boltz with a slight α3 edge. Written up in CONCLUSIONS.md 4d.
