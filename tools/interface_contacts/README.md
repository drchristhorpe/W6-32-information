# interface-contacts

Map a protein–protein interface as a heavy-atom contact set. Given a structure
and two groups of chains (an **epitope** side and a **paratope** side), it finds
every heavy-atom contact within a distance cutoff (BioPython `NeighborSearch`)
and reports the per-residue footprint on each side plus the underlying atom
pairs. Input may be `.pdb` or `.cif`.

Contacts are **distance-based only** — there is no H-bond / salt-bridge / π
chemical typing.

## Contract

- Args: `input_filepath` and `output_folder`.
- Writes `<stem>_contacts.json`: metadata, `epitope_footprint`,
  `paratope_footprint` (each a list of residues with `min_distance`,
  `n_contacts`, `partners`), and the sorted `contacts` atom-pair list.

## Usage

```bash
uv run interface-contacts <input.pdb|.cif> <output_folder> \
    --epitope-chains A,B --paratope-chains H [--cutoff 5.0]
```

The cutoff is configurable; the default is **5.0 Å** (module constant
`DEFAULT_CUTOFF`), chosen to give the design phase a wider interface description.

For the W6/32 complex, the paratope is the antibody heavy chain (H/G) and the
epitope is the HLA heavy + β2m chains (A,B / D,E); the light chain is excluded
because it is sterically hindered in the single-chain-trimer format.

## Standalone use elsewhere

Self-contained folder. Copy it out and `uv run interface-contacts ...` or
`pip install .` (entry point `interface-contacts`). Only runtime dependency is
BioPython.
