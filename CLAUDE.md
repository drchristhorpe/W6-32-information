# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

A **data repository** (not a software project — no build, lint, or test). It collects structural-biology data about **W6/32** (a pan-HLA-class-I monoclonal antibody) and about **HLA single-chain trimers (SCTs)**, for protein-design tasks. There is no code to run; work here is reading, parsing, comparing, and organizing structure files and their associated confidence metrics.

Key biological context:
- **W6/32** — pan-class-I anti-HLA antibody. Its heavy-chain structures live in `structures/`.
- **PRAME peptide** = `SLLQHLIGL`, presented by HLA-A*02:01 (see `README.md`).
- **Single-chain trimer (SCT)** — peptide + β2-microglobulin + HLA heavy (alpha) chain fused into one polypeptide via GGGGS linkers, in that order: `peptide–(GGGGS)x3–β2m–(GGGGS)x4–HLA_alpha`. All `sct_predictions/` entries are SCTs of this form. Segment boundaries for a given prediction are recorded in the esmfold `summary.json` `segments` array.

## Directory layout

- `structures/` — reference/experimental structures, mirrored in two formats: `cif/` and `pdb/`. Same basenames in both. Naming: `<allele>__<peptide>` for HLA–peptide complexes (e.g. `hla_a_02_01__SLLQHLIGL`), and `w632__heavy_chain__<n>` for W6/32 antibody heavy chains.
- `sct_predictions/raw/` — SCT structure predictions, one subdirectory **per predictor**: `alphafold3/`, `boltz/`, `esmfold2/`, `esmfold2-fast/`. The same set of allele×peptide SCTs is predicted by each (AlphaFold3 has a couple extra). Output file layout differs per predictor — see below.
- `sct_predictions/aligned/` — intended for structurally-aligned predictions; currently empty.
- `interface_description/` — currently empty.
- `PLAN.md` — currently empty (working plan for the experiment).
- `w632_sct_prame.pse` — PyMOL session file (binary).

## Naming conventions (differ by predictor — watch the casing and separators)

The same SCT is named three different ways:
- **AlphaFold3**: `fold_<allele>_single_chain_trimer_<peptide>` — all lowercase, single underscores, `fold_` prefix. e.g. `fold_hla_a_02_01_single_chain_trimer_sllqhligl`.
- **Boltz / esmfold2 / esmfold2-fast**: `<allele>__single_chain_trimer__<PEPTIDE>` — **double** underscores between parts, peptide in **UPPERCASE**. e.g. `hla_a_02_01__single_chain_trimer__SLLQHLIGL`.

When matching predictions across predictors, normalize both the separators and the peptide casing.

## Per-predictor output formats

Each predictor writes a different structure on disk and reports confidence differently. When extracting/comparing metrics, branch on the predictor:

- **alphafold3/`<job>/`** — 5 models per job (`*_model_0..4.cif`), with matching `*_full_data_0..4.json` and `*_summary_confidences_0..4.json`. `*_job_request.json` holds the input sequence (single `proteinChain` with the full fused SCT sequence) and `modelSeeds`. Also `msas/` (`.a3m`), `templates/` (`.cif`), `terms_of_use.md`.
- **boltz/`<job>/`** — `run.json` (input sequence, model `boltz-2.1`, timestamps, and full `output` block with per-sample + `best_sample` metrics) plus `outputs/files/prediction/`: `sample_0_predicted_structure.cif`, `sample_0_pae.npz`, `metrics.json`. Metrics include `structure_confidence`, `ptm`, `iptm`, `complex_plddt`, `complex_pde`, etc. (iptm fields are 0 because an SCT is a single chain).
- **esmfold2/** and **esmfold2-fast/** (`<job>/`) — `<job>.pdb` (structure), `<job>.confidence.json`, `<job>.raw_response.json`, and `summary.json`. `summary.json` is the richest metadata: `complex_type`, `model`, `seed`, `sequence_length`, the **`segments`** array (peptide/b2m/hla_alpha start/end indices), `folding_config`, timing, and top-level `plddt_mean`, `ptm`, `interface_ptm`. `esmfold2-fast` is the faster/lower-setting variant with the same layout.

Note cross-predictor metric scales differ: boltz/esmfold report pLDDT and pTM on a **0–1** scale here.

## Tooling (Python analysis code)

Code lives alongside the data and is managed with **`uv`** (Python 3.14). This
directory is a member of the home-level uv workspace (`~/pyproject.toml`), so it
shares the workspace `.venv`. Consequence: `uv sync` here reshapes that shared
env to this project's deps (temporarily removing sibling projects' packages like
torch); they come back on the next `uv run` in those projects.

**Each analysis tool is a self-contained folder under `tools/<name>/`** with its
own `pyproject.toml`, `README.md`, and module — so it can be copied out and used
elsewhere (`uv run`/`pip install`) unchanged. Tools are wired into this repo's
shared env as **editable path dependencies** of the root `w632-tools` project
(see `[tool.uv.sources]` in `pyproject.toml`). Conventions for every tool: takes
an `input_filepath` + `output_folder`, writes machine-readable JSON, no hidden
global state.

Run a tool via its console entry point in the shared env, e.g.:
```bash
uv run align-structure <input.pdb|.cif> <output_folder> [--reference ... --ref-residues 1-180 --mobile-chain A]
```
Or standalone (isolated env, proves extractability): `uv run --isolated --project tools/align_structure align-structure ...`.

**Each tool also ships a matching skill** at `.claude/skills/<name>/SKILL.md`,
a thin wrapper over the tool's CLI (single source of truth — the logic stays in
`tools/<name>/`, the skill just documents when/how to invoke it).

Tools (see each folder's README for the full contract):
- `tools/align_structure/` (skill: `align-structure`) — superpose onto the
  canonical frame; fit residues are matched by **local sequence alignment**, so it
  aligns both separate-chain experimental complexes and single-chain SCT
  predictions. Auto-detects `.pdb`/`.cif`.
- `tools/interface_contacts/` (skill: `interface-contacts`) — heavy-atom contact
  footprint between epitope/paratope chain groups (BioPython `NeighborSearch`,
  configurable cutoff, default 5.0 Å). Distance-based only.
- `tools/compare_structures/` (skill: `compare-structures`) — global + per-residue
  deviation over a selection; residues paired by **local sequence alignment**;
  in-frame by default (`--superpose` for Kabsch); `--footprint` to restrict to an
  interface. BioPandas + BioPython. Takes PDB (align CIF predictions first).
- `tools/suggest_disulphides/` (skill: `suggest-disulphides`) — propose engineerable
  disulphides between residue groups; models Sγ rotamers and checks Sγ–Sγ/χ3
  geometry (not just Cβ–Cβ); `--exclude-footprint`, pLDDT-aware. BioPython + NumPy.
- `tools/mutate_sequence/` (skill: `mutate-sequence`) — apply validated point
  mutations to a sequence from a JSON spec (`{sequence, mutations:[{position,from,to}]}`,
  1-based; `from` is checked). Stdlib only. Outputs mutated JSON + FASTA.
- `tools/prepare_cif_for_boltz/` (skill: `prepare-cif-for-boltz`) — reconstruct
  mmCIF polymer metadata (`_entity_poly_seq`, `_struct_asym`, `_atom_site.auth_seq_id`)
  stripped by viewer exports, so BoltzGen accepts the file. gemmi. Outputs
  `<stem>_for_boltz.cif`.

Orchestration drivers live in `analysis/` (e.g. `align_all_predictions.py`), and
import the tools as libraries. Experiment outputs: design-question answers go in
[CONCLUSIONS.md](CONCLUSIONS.md); interface footprints in `interface_description/`;
aligned predictions in `sct_predictions/aligned/` (gitignored).

Workflow discipline (from [EXPERIMENTAL_DESIGN.md](EXPERIMENTAL_DESIGN.md)): plan
in [PLAN.md](PLAN.md) → execute → update this file + [CHANGELOG.md](CHANGELOG.md)
(lab notebook) → commit only after the user approves the commit message.

## Working notes

- This is a git repo on `main`. The bulky `.cif`/`.pdb`/`.npz`/`.pse` files and macOS `.DS_Store` files are data, not artifacts to regenerate — do not delete them when "cleaning up".
- Prefer the structured `summary.json` (esmfold) / `metrics.json` / `run.json` (boltz) / `*_summary_confidences_*.json` (AF3) for reading confidence values rather than parsing the coordinate files.
