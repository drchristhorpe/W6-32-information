---
name: interface-contacts
description: Map a protein-protein interface as a heavy-atom contact set between two groups of chains (epitope vs paratope), using BioPython NeighborSearch with a configurable distance cutoff (default 5.0 Å). Use when asked for the contact residues / footprint of an antibody-antigen or other interface, e.g. the W6/32 heavy-chain footprint on HLA + β2m, or to define interface residues to restrict a downstream RMSD comparison. Outputs a contacts JSON with per-side residue footprints. Distance-based only (no H-bond/π typing).
---

# interface-contacts

Wraps the self-contained `interface_contacts` tool (`tools/interface_contacts/`).
Finds every heavy-atom contact within a cutoff between an **epitope** chain group
and a **paratope** chain group, and reports the per-residue footprint on each side
plus the atom-pair list. Input may be `.pdb` or `.cif`. Contacts are distance-based
only — no H-bond / salt-bridge / π classification.

## When to use

- Define the W6/32 footprint on the HLA+β2m surface (the interface residues that
  later restrict 4a/4b comparisons).
- Get the contact residues of any antibody–antigen or chain–chain interface.

## How to run

```bash
uv run interface-contacts <input.pdb|.cif> <output_folder> \
    --epitope-chains A,B --paratope-chains H [--cutoff 5.0]
```

- `--cutoff` defaults to **5.0 Å** (wider footprint for design work).
- For the W6/32 complex: paratope = antibody heavy chain (`H` in copy 1, `G` in
  copy 2); epitope = HLA heavy + β2m (`A,B` in copy 1, `D,E` in copy 2). The light
  chain is excluded (sterically hindered in the single-chain-trimer format).

## Output contract

Into `<output_folder>`: `<stem>_contacts.json` with `epitope_chains`,
`paratope_chains`, `cutoff`, `n_contact_atom_pairs`, `epitope_footprint` and
`paratope_footprint` (each residue: `chain`, `resseq`, `resname`, `min_distance`,
`n_contacts`, `partners`), and the distance-sorted `contacts` atom-pair list.

## Notes

- A correct W6/32 footprint spans the HLA α3 domain (~res 222–243) and β2m
  (~res 1–6, 58–61) — a good sanity check that chains were assigned correctly.
- Use `epitope_footprint` residue keys to drive `compare-structures --selection`.
