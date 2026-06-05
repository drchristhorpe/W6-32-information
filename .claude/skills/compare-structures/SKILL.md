---
name: compare-structures
description: Measure fine-scale global and per-residue deviation between two protein structures over a residue selection, with residues paired by local sequence alignment (handles different numbering / chain layouts, e.g. experimental complex vs single-chain trimer). Use when asked how similar two structures or a prediction-vs-experiment pair are, for per-residue RMSD, peptide-conformation comparison, or interface-restricted change. Measures in the input frame by default (assumes pre-aligned); --superpose for a Kabsch fit. Outputs a per-residue RMSD JSON.
---

# compare-structures

Wraps the self-contained `compare_structures` tool (`tools/compare_structures/`).
Reports global + per-residue deviation between a reference and a mobile structure
over a residue selection. Residues are paired by **local sequence alignment** per
chain, so it works across different numbering and chain layouts.

By default it measures deviation **in the input coordinate frame** (no
re-superposition) — correct when the structures are already in a common frame
(via `align-structure`, or already co-framed). Use `--superpose` to Kabsch-fit
the selection first.

## When to use

- 4a: change at the W6/32 interface between bound and apo B\*27:05 (same frame,
  `--footprint`).
- 4b: prediction vs experiment — global groove, peptide-only, or interface-only
  RMSD of an aligned SCT against the experimental structure.
- Any per-residue prediction-quality or conformational-change question.

## How to run

```bash
uv run compare-structures <reference.pdb> <mobile.pdb> <output_folder> \
    [--ref-chains A,B --mobile-chains A,B] \
    [--footprint interface_description/<...>_contacts.json] \
    [--atoms ca|heavy] [--superpose]
```

- `--ref-chains`/`--mobile-chains`: parallel lists; default = chains common to
  both. For an SCT mobile (single fused chain), repeat the chain, e.g.
  `--ref-chains A,B --mobile-chains A,A` — local alignment locates each region.
- `--footprint`: restrict to the `epitope_footprint` of an `interface-contacts`
  result (matched by ref chain + resseq).

## Output contract

`<refstem>__vs__<mobstem>_compare.json`: `global_rmsd`, `n_residues`, `n_atoms`,
`superposed`, `max_residue_rmsd`, and a `per_residue` table (`ref_chain`,
`resseq`, `resname`, `rmsd`, `max_atom_dev`, `n_atoms`).

## Notes

- Default in-frame measurement captures real displacement (including
  rigid-body/hinge motion relative to the alignment frame); `--superpose` isolates
  internal conformational difference. Pick deliberately per question.
- Inputs are PDB (BioPandas). Align CIF predictions to PDB with `align-structure`
  first (it outputs PDB).
