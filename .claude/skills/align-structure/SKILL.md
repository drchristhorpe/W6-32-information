---
name: align-structure
description: Superpose a protein structure (.pdb or .cif) onto a canonical reference frame, fitting on a sequence-matched residue selection (default MHC-I α1/α2 Cα, res 1–180). Use when asked to align/superpose an HLA/MHC structure or a single-chain-trimer (SCT) prediction onto the canonical 1HHK frame, or to move predictions into a common coordinate frame before RMSD/interface comparison. Outputs an aligned PDB + transform JSON.
---

# align-structure

Wraps the self-contained `align_structure` tool (`tools/align_structure/`). It
superposes a mobile structure onto a reference, choosing the fit atoms by **local
sequence alignment** of the reference selection against the mobile chain — so it
aligns both separate-chain experimental complexes and single-chain SCT
predictions (whose HLA α chain is fused mid-sequence with its own numbering).
Input may be `.pdb` or `.cif`; output is always PDB.

## When to use

- Bringing Boltz / AlphaFold3 / ESMFold SCT predictions into the canonical frame
  (`structures/1hhk_1_aligned.pdb`) so they can be compared to experiment.
- Superposing any HLA/MHC structure onto the α1/α2 platform.

## How to run

From the experiment root, via the shared workspace env:

```bash
uv run align-structure <input.pdb|.cif> <output_folder> \
    --reference structures/1hhk_1_aligned.pdb \
    --ref-chain A --ref-residues 1-180 \
    --mobile-chain A [--mobile-residues LO-HI]
```

Defaults: reference `structures/1hhk_1_aligned.pdb`, fit selection chain A
residues `1-180` (MHC-I α1/α2). For an SCT prediction pass `--mobile-chain` as
the single chain id (Boltz/AF3 use `A`); the α1/α2 region is located
automatically by sequence.

Standalone (isolated env, no repo install needed):
`uv run --isolated --project tools/align_structure align-structure ...`.

## Output contract

Into `<output_folder>`:
- `<stem>_aligned.pdb` — the transformed structure.
- `<stem>_align.json` — `rotation`, `translation`, `n_matched_residues`,
  `alignment_score`, fit `rmsd`, and the selections used.

## Notes

- A small fit RMSD on an already-co-framed experimental structure (~0.7 Å)
  reflects genuine structural difference, not a frame error.
- Needs ≥3 matched residues or it errors; check `n_matched_residues` in the JSON.
