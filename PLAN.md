# PLAN

Working plan for the W6/32 information experiment. Authored by Claude from [EXPERIMENTAL_DESIGN.md](EXPERIMENTAL_DESIGN.md). Execute only after approval; update [CLAUDE.md](CLAUDE.md) and `CHANGELOG.md` (lab notebook) after each piece of work; commit on approval of the commit message.

## Decisions locked in

- **Environment:** Python 3.14, `uv` for package management.
- **Libraries:** BioPython (alignment, parsing, contacts) + BioPandas (RMSD / dataframes). **Arpeggio dropped** per decision — contacts via BioPython `NeighborSearch` (distance-based heavy-atom contacts; no H-bond/π chemical typing).
- **Canonical reference frame:** [structures/1hhk_1_aligned.pdb](structures/1hhk_1_aligned.pdb) — chain A = HLA-A\*02:01 heavy (res 1–275), B = β2m, C = peptide.
- **Superposition selection ("antigen binding domain"):** chain A α1/α2 groove **Cα, residues 1–180**.
- **Tool contract (every tool):** self-contained in its own folder `tools/<name>/` (own `pyproject.toml` + `README.md` + module), so it can be lifted out as a standalone library; takes an **input filepath** and an **output folder** as arguments; writes machine-readable JSON (+ derived structure files) into the output folder; no hidden global state.
- **Skill per tool:** every tool also ships a matching Claude skill at `.claude/skills/<name>/SKILL.md` that wraps the tool's CLI (single source of truth — no duplicated logic).

## Open prerequisites

- None blocking. Canonical reference confirmed present.
- Note: predictions in `sct_predictions/raw/` are **not yet aligned** — they must pass through Tool 1 before any cross-structure comparison.

## Reference data map (verified on disk)

| Structure | Chains | Use |
|---|---|---|
| `structures/1hhk_1_aligned.pdb` | A heavy / B β2m / C peptide | **canonical frame & alignment target** |
| `structures/pdb/hla_a_02_01__SLLQHLIGL.pdb` | A / B / P | A\*02:01 + PRAME experimental — **ground truth (4b)** |
| `structures/pdb/hla_b_27_05__FRYNGLIHR__w632__1.pdb` | A,B,C + H,L | B\*27:05 complex + W6/32 heavy(H)+light(L), copy 1 |
| `structures/pdb/hla_b_27_05__FRYNGLIHR__w632__2.pdb` | D,E,F + G,K | copy 2 |
| `structures/pdb/hla_b_27_05__RRFSRSPIRR.pdb` | A / B / C | apo B\*27:05 + RRFSRSPIRR — **comparator (4a)** |
| `structures/pdb/w632__heavy_chain__{1,2}.pdb` | H / G | isolated W6/32 heavy chains |

SCT construct order (all predictors, identical sequence): `peptide –(GGGGS)3– β2m –(GGGGS)4– HLA_alpha`, single chain. Segment offsets derived from the construct (cross-checked against esmfold `summary.json` `segments`).

## Tools to build

### Tool 1 — `align_structure` (build first)
Superpose a structure onto the canonical reference.
- **In:** `input_filepath` (PDB), `output_folder`; opts: `--reference` (default `1hhk_1_aligned.pdb`), `--ref-chain A`, `--mobile-chain`, `--residues 1-180`, atom `CA`.
- **How:** BioPython `Superimposer` on matched α1/α2 Cα. For SCT predictions, the mobile "heavy-chain" selection is the HLA_alpha segment of the single chain (mapped via construct offsets).
- **Out:** aligned PDB (→ `sct_predictions/aligned/` for predictions), `align.json` (rotation, translation, n_atoms, alignment RMSD).

### Tool 2 — `interface_contacts`
Map an interface as a residue/atom contact set.
- **In:** `input_filepath`, `output_folder`; `--epitope-chains`, `--paratope-chains`, `--cutoff` (configurable; **default 5.0 Å** heavy-atom — wider footprint for the design phase).
- **How:** BioPython `NeighborSearch` over heavy atoms; report contacting residue pairs, per-residue min distance, and the epitope-side footprint residue list.
- **Out:** `contacts.json` (pairs + footprint), written for `interface_description/`.
- **W6/32 usage:** paratope = heavy chain (H/G); epitope = HLA heavy + β2m (A,B / D,E). Light chain excluded (sterically hindered in SCT).

### Tool 3 — `compare_structures`
Fine-scale difference between two structures.
- **In:** `reference_filepath`, `input_filepath`, `output_folder`; `--selection` (chains/residue ranges, optional footprint list from Tool 2), per-residue toggle.
- **How:** BioPandas to build Cα (and optionally all-atom) frames on a residue-equivalence mapping; global RMSD + per-residue RMSD/displacement. Assumes inputs share the canonical frame (via Tool 1).
- **Out:** `compare.json` (global RMSD, per-residue table, max-displacement residues).

## Execution plan for the design questions

**4a — Does W6/32 binding change the B\*27:05 conformation? (interface only)**
1. `interface_contacts` on `hla_b_27_05__FRYNGLIHR__w632__1` → W6/32-heavy footprint on {A,B}.
2. `compare_structures` between the W6/32-bound B\*27:05 and apo `hla_b_27_05__RRFSRSPIRR`, **restricted to the footprint residues** (peptide ignored). Both already in one frame → no Tool 1 needed. Output: per-residue conformational change at the W6/32 interface.

**4b — Is the Boltz A\*02:01/PRAME SCT good enough as a design input?**
1. `align_structure` the Boltz SCT (`sct_predictions/raw/boltz/hla_a_02_01__single_chain_trimer__SLLQHLIGL`) onto canonical.
2. `compare_structures` vs `hla_a_02_01__SLLQHLIGL` experimental at three scopes: (a) global groove, (b) **peptide only** (SLLQHLIGL conformation), (c) **W6/32 interface footprint** (footprint defined on the experimental A\*02:01 via Tool 2 against a W6/32 heavy chain superposed from the B\*27:05 complex).
3. Cross-read against Boltz `metrics.json` confidence. Optional: same for AF3/ESMFold as sanity context.

**4c — Best location for an α3 ↔ β2m-linker disulphide (extension)**
- Geometric scan over residue pairs (α3 domain of HLA_alpha × the β2m→α1 linker) for Cβ–Cβ distance ≈ 3.5–4.5 Å and compatible χ geometry, on the aligned SCT. Reuse Tool 1 output; report ranked candidate pairs. Build a small `disulphide_scan` helper only if 4a/4b tooling doesn't already cover it.

## Order of work (milestones)

1. Scaffold `uv` project (pyproject, src layout, BioPython + BioPandas). Smoke-test parse of one PDB.
2. Tool 1 `align_structure` + align all SCT predictions into `sct_predictions/aligned/`.
3. Tool 2 `interface_contacts` + produce W6/32 footprints into `interface_description/`.
4. Tool 3 `compare_structures`.
5. Apply to 4a, then 4b, then 4c (extension).
6. After each milestone: update CLAUDE.md + CHANGELOG.md; request commit approval.

Each tool ships with a CLI entry point and is importable, so it can later be lifted into a standalone library or skill file unchanged.
