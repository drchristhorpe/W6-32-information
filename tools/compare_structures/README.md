# compare-structures

Fine-scale difference between two structures: global and per-residue deviation
over a residue selection. Residues are paired by **local sequence alignment**
per chain, so it compares structures with different numbering or chain layouts
(e.g. a separate-chain experimental complex vs a single-chain trimer).

By default it measures deviation **in the input coordinate frame, without
re-superposing** — correct when the inputs are already in a common frame (via
`align-structure`, or already co-framed). Pass `--superpose` to Kabsch-fit the
selection first (pure internal conformational difference).

Coordinates/RMSD use BioPandas; sequence pairing uses BioPython.

## Contract

- Args: `reference_filepath`, `input_filepath` (mobile), `output_folder`.
- Writes `<refstem>__vs__<mobstem>_compare.json`: `global_rmsd`, `n_residues`,
  `max_residue_rmsd`, and a `per_residue` table (`ref_chain`, `resseq`,
  `resname`, `rmsd`, `max_atom_dev`, `n_atoms`).

## Usage

```bash
uv run compare-structures <reference.pdb> <mobile.pdb> <output_folder> \
    [--ref-chains A,B --mobile-chains A,B] \
    [--footprint path/to/contacts.json] \
    [--atoms ca|heavy] [--superpose]
```

- `--ref-chains`/`--mobile-chains`: parallel lists; default = chains common to
  both. For an SCT mobile, point multiple ref chains at the single SCT chain
  (e.g. `--ref-chains A,B --mobile-chains A,A`); local alignment finds each region.
- `--footprint`: a `contacts.json` from `interface-contacts`; restricts the
  comparison to its `epitope_footprint` residues (by ref chain + resseq).
- `--atoms ca` (default) or `heavy` (all heavy atoms, matched by name per residue).

## Standalone use elsewhere

Self-contained folder. Copy it out and `uv run compare-structures ...` or
`pip install .`. Runtime deps: BioPandas, BioPython, NumPy.
