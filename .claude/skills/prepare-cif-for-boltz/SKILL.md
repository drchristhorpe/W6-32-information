---
name: prepare-cif-for-boltz
description: Reconstruct mmCIF polymer metadata that BoltzGen structure-template design requires but viewer exports strip (_entity_poly_seq.entity_id/num/mon_id, _struct_asym.id/entity_id, _atom_site.auth_seq_id). Use when Boltz/BoltzGen returns a 400 "missing polymer metadata required for structure-template protein design" error, or to repair a stripped/converted mmCIF or PDB into a clean PDBx/mmCIF. Rebuilds from coordinates with gemmi and validates the required items. Outputs <stem>_for_boltz.cif.
---

# prepare-cif-for-boltz

Wraps the self-contained `prepare_cif_for_boltz` tool (`tools/prepare_cif_for_boltz/`).
BoltzGen needs polymer metadata that viewer exports often drop:
`_entity_poly_seq.{entity_id,num,mon_id}`, `_struct_asym.{id,entity_id}`,
`_atom_site.auth_seq_id`. This rebuilds them from the coordinates with gemmi and
writes a clean PDBx/mmCIF, validating the required items (errors if any remain
missing).

## When to use

- BoltzGen 400: "Binder structure mmCIF is missing polymer metadata required for
  structure-template protein design ... use an original PDBx/mmCIF that preserves
  polymer sequence metadata."
- Repairing any mmCIF/PDB exported through a structure viewer before Boltz.

## How to run

```bash
uv run prepare-cif-for-boltz <input.cif|.pdb> <output_folder>
```

E.g.:
```bash
uv run prepare-cif-for-boltz structures/cif/w632_heavy_chain_variable.cif structures/cif/
# -> structures/cif/w632_heavy_chain_variable_for_boltz.cif
```

## Output contract

Writes `<stem>_for_boltz.cif`. The tool returns / reports the chains, polymer
entities, and a per-item presence check for the six required metadata fields;
it exits non-zero if any are still missing.

## Notes

- Rebuilds the sequence from the *modelled* residues; unmodelled gaps are not
  invented. That is what BoltzGen needs for a structure template.
- Method: gemmi `setup_entities()` + `assign_label_seq_id()`, then each polymer
  entity's `full_sequence` is filled from its subchain to regenerate
  `_entity_poly_seq`.
