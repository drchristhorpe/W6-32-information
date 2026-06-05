# prepare-cif-for-boltz

Reconstruct the polymer metadata that **BoltzGen** structure-template protein
design requires but that viewer exports often strip:

- `_entity_poly_seq.entity_id` / `.num` / `.mon_id`
- `_struct_asym.id` / `.entity_id`
- `_atom_site.auth_seq_id`

It rebuilds them from the coordinates with **gemmi** and writes a clean PDBx/mmCIF,
then validates that all required items are present (exit code 1 if not).

## Contract

- Args: `input_filepath` (`.cif` or `.pdb`) and `output_folder`.
- Writes `<stem>_for_boltz.cif`.

## Usage

```bash
uv run prepare-cif-for-boltz <input.cif|.pdb> <output_folder>
```

E.g. for the W6/32 VH region:
```bash
uv run prepare-cif-for-boltz structures/cif/w632_heavy_chain_variable.cif structures/cif/
# -> structures/cif/w632_heavy_chain_variable_for_boltz.cif
```

## Method

`gemmi.read_structure` → `setup_entities()` (assign subchains + entities) →
`assign_label_seq_id()` → populate each polymer entity's `full_sequence` from its
subchain (this regenerates `_entity_poly_seq`) → `make_mmcif_document().write_file`.

## Standalone use elsewhere

Self-contained folder. Copy it out and `uv run prepare-cif-for-boltz ...` or
`pip install .`. Runtime dependency: gemmi.
