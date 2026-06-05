"""prepare_cif_for_boltz -- reconstruct mmCIF polymer metadata for BoltzGen.

Self-contained tool (own helpers + CLI) so it can be lifted into a standalone
library or skill file unchanged.

BoltzGen structure-template protein design requires polymer metadata that viewer
exports often strip:
    _entity_poly_seq.entity_id / .num / .mon_id
    _struct_asym.id / .entity_id
    _atom_site.auth_seq_id
This tool rebuilds that metadata from the coordinates with gemmi and writes a
clean PDBx/mmCIF, then validates the required items are present.

Method: read structure -> setup_entities() (assign subchains/entities) ->
assign_label_seq_id() -> populate each polymer entity's full_sequence from its
subchain (this is what regenerates _entity_poly_seq) -> write mmCIF.

Contract: takes an input filepath (.cif or .pdb) and an output folder; writes
`<stem>_for_boltz.cif` into the output folder.

Example (job 8):
    python -m prepare_cif_for_boltz \\
        structures/cif/w632_heavy_chain_variable.cif structures/cif/
"""

from __future__ import annotations

import argparse
from pathlib import Path

import gemmi

REQUIRED_ITEMS = [
    "_entity_poly_seq.entity_id",
    "_entity_poly_seq.num",
    "_entity_poly_seq.mon_id",
    "_atom_site.auth_seq_id",
    "_struct_asym.id",
    "_struct_asym.entity_id",
]


def prepare_cif_for_boltz(input_filepath: str | Path, output_folder: str | Path) -> dict:
    input_filepath = Path(input_filepath)
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    structure = gemmi.read_structure(str(input_filepath))
    structure.setup_entities()        # assign subchains + entities
    structure.assign_label_seq_id()   # label_seq_id for _atom_site

    # Regenerate _entity_poly_seq: gemmi only writes it from each polymer entity's
    # full_sequence, which is empty when the input had no sequence metadata.
    model = structure[0]
    for entity in structure.entities:
        if entity.entity_type != gemmi.EntityType.Polymer:
            continue
        for subchain_id in entity.subchains:
            subchain = model.get_subchain(subchain_id)
            if subchain and len(subchain):
                entity.full_sequence = [res.name for res in subchain]
                break

    out_path = output_folder / f"{input_filepath.stem}_for_boltz.cif"
    structure.make_mmcif_document().write_file(str(out_path))

    text = out_path.read_text()
    present = {item: (item in text) for item in REQUIRED_ITEMS}
    missing = [item for item, ok in present.items() if not ok]

    result = {
        "input": str(input_filepath),
        "output": str(out_path),
        "chains": [f"{c.name}:{len(c)}" for c in model],
        "polymer_entities": [e.name for e in structure.entities if e.entity_type == gemmi.EntityType.Polymer],
        "required_items_present": present,
        "missing_items": missing,
        "ok": not missing,
    }
    return result


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("input_filepath", help="structure (.cif or .pdb) to repair")
    ap.add_argument("output_folder", help="folder for <stem>_for_boltz.cif")
    args = ap.parse_args(argv)

    result = prepare_cif_for_boltz(args.input_filepath, args.output_folder)
    status = "OK — all required metadata present" if result["ok"] else f"MISSING: {result['missing_items']}"
    print(
        f"{Path(args.input_filepath).name}: chains {result['chains']}; "
        f"{status} -> {result['output']}"
    )
    if not result["ok"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
