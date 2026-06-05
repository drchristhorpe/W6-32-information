# align-structure

Superpose a structure onto a canonical reference frame and write the transformed
structure plus the transform metadata.

The fitting selection is matched by **local sequence alignment** of the reference
selection against the mobile chain, not by assuming residue numbers line up. This
lets one tool align both separate-chain experimental complexes and single-chain
single-chain-trimer (SCT) predictions (whose HLA alpha chain is fused mid-sequence
with its own numbering). Input may be `.pdb` or `.cif` (parser auto-selected);
output is always PDB.

## Contract

- Args: `input_filepath` (mobile structure) and `output_folder`.
- Writes `<stem>_aligned.pdb` and `<stem>_align.json` (rotation, translation,
  matched-residue count, fit RMSD) into the output folder.

## Usage

```bash
uv run align-structure <input> <output_folder> \
    --reference structures/1hhk_1_aligned.pdb \
    --ref-chain A --ref-residues 1-180 \
    --mobile-chain A [--mobile-residues LO-HI]
```

Default reference is `structures/1hhk_1_aligned.pdb` and default fit selection is
chain A residues 1–180 (the MHC-I α1/α2 peptide-binding platform).

## Standalone use elsewhere

This folder is self-contained. Copy it out and either `uv run align-structure ...`
or `pip install .` (entry point `align-structure`). Its only runtime dependency is
BioPython.
