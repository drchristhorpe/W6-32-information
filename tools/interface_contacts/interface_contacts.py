"""interface_contacts -- map a protein-protein interface as a contact set.

Self-contained tool (own helpers + CLI) so it can be lifted into a standalone
library or skill file unchanged.

Given a structure and two groups of chains (an "epitope" side and a "paratope"
side), it finds every heavy-atom contact within a distance cutoff and reports:
  * the per-residue epitope footprint (the residues the paratope touches),
  * the symmetric paratope footprint,
  * the underlying atom-atom contact pairs.

Contacts are distance-based only (BioPython NeighborSearch over heavy atoms);
there is no H-bond / salt-bridge / pi chemical typing.

Contract: takes an input filepath and an output folder; writes a single
`<stem>_contacts.json` into the output folder.

Example (W6/32 heavy-chain footprint on HLA + b2m):
    python -m interface_contacts \\
        structures/pdb/hla_b_27_05__FRYNGLIHR__w632__1.pdb out/ \\
        --epitope-chains A,B --paratope-chains H
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from Bio.PDB import MMCIFParser, NeighborSearch, PDBParser

# Default heavy-atom contact distance. Widened from 4.0 to 5.0 Å so the design
# phase gets a broader interface description. Configurable per call.
DEFAULT_CUTOFF = 5.0


def _load_structure(path: str | Path, name: str):
    path = Path(path)
    parser = MMCIFParser(QUIET=True) if path.suffix.lower() in (".cif", ".mmcif") else PDBParser(QUIET=True)
    return parser.get_structure(name, str(path))


def _is_heavy(atom) -> bool:
    el = (atom.element or "").strip().upper()
    if el:
        return el not in ("H", "D")
    return not atom.get_id().lstrip().startswith(("H", "D"))


def _res_key(residue) -> str:
    chain = residue.get_parent().id
    het, resseq, icode = residue.id
    return f"{chain}:{residue.resname}{resseq}{icode}".strip()


def _res_record(residue) -> dict:
    chain = residue.get_parent().id
    _, resseq, icode = residue.id
    return {
        "chain": chain,
        "resseq": resseq,
        "icode": icode.strip(),
        "resname": residue.resname,
        "key": _res_key(residue),
    }


def interface_contacts(
    input_filepath: str | Path,
    output_folder: str | Path,
    epitope_chains: list[str],
    paratope_chains: list[str],
    cutoff: float = DEFAULT_CUTOFF,
) -> dict:
    """Find heavy-atom contacts between the epitope and paratope chain groups and
    write `<stem>_contacts.json` into `output_folder`."""
    input_filepath = Path(input_filepath)
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    epitope_chains = list(epitope_chains)
    paratope_chains = list(paratope_chains)
    group_of = {c: "epitope" for c in epitope_chains}
    group_of.update({c: "paratope" for c in paratope_chains})

    model = _load_structure(input_filepath, input_filepath.stem)[0]

    atoms = []
    for chain in model:
        if chain.id not in group_of:
            continue
        for residue in chain:
            if residue.id[0] != " ":  # skip water / hetero
                continue
            for atom in residue:
                if _is_heavy(atom):
                    atoms.append(atom)
    if not atoms:
        raise ValueError("No heavy atoms found for the requested chains.")

    ns = NeighborSearch(atoms)
    contacts = []
    # per-residue aggregation: key -> {record, min_distance, partners:set}
    epi_foot: dict[str, dict] = {}
    par_foot: dict[str, dict] = {}

    for a, b in ns.search_all(cutoff, level="A"):
        ca = a.get_parent().get_parent().id
        cb = b.get_parent().get_parent().id
        ga, gb = group_of[ca], group_of[cb]
        if ga == gb:  # same side of the interface
            continue
        epi_atom, par_atom = (a, b) if ga == "epitope" else (b, a)
        dist = epi_atom - par_atom
        epi_res = epi_atom.get_parent()
        par_res = par_atom.get_parent()

        contacts.append({
            "epitope": {**_res_record(epi_res), "atom": epi_atom.get_id()},
            "paratope": {**_res_record(par_res), "atom": par_atom.get_id()},
            "distance": round(float(dist), 3),
        })
        _accumulate(epi_foot, epi_res, par_res, dist)
        _accumulate(par_foot, par_res, epi_res, dist)

    result = {
        "input": str(input_filepath),
        "cutoff": cutoff,
        "epitope_chains": epitope_chains,
        "paratope_chains": paratope_chains,
        "n_contact_atom_pairs": len(contacts),
        "epitope_footprint": _finalise(epi_foot),
        "paratope_footprint": _finalise(par_foot),
        "contacts": sorted(contacts, key=lambda c: c["distance"]),
    }
    out_path = output_folder / f"{input_filepath.stem}_contacts.json"
    with open(out_path, "w") as fh:
        json.dump(result, fh, indent=2)
    result["_output"] = str(out_path)
    return result


def _accumulate(table: dict, residue, partner, dist: float) -> None:
    key = _res_key(residue)
    entry = table.get(key)
    if entry is None:
        entry = {**_res_record(residue), "min_distance": float(dist), "partners": set()}
        table[key] = entry
    entry["min_distance"] = min(entry["min_distance"], float(dist))
    entry["partners"].add(_res_key(partner))


def _finalise(table: dict) -> list[dict]:
    out = []
    for entry in table.values():
        partners = sorted(entry.pop("partners"))
        entry["n_contacts"] = len(partners)
        entry["partners"] = partners
        entry["min_distance"] = round(entry["min_distance"], 3)
        out.append(entry)
    # sort by chain then residue number for a stable, readable footprint
    return sorted(out, key=lambda e: (e["chain"], e["resseq"]))


def _chain_list(spec: str) -> list[str]:
    return [c.strip() for c in spec.split(",") if c.strip()]


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("input_filepath", help="structure (.pdb or .cif)")
    ap.add_argument("output_folder", help="folder for <stem>_contacts.json")
    ap.add_argument("--epitope-chains", required=True, help="comma-separated, e.g. A,B")
    ap.add_argument("--paratope-chains", required=True, help="comma-separated, e.g. H")
    ap.add_argument(
        "--cutoff", type=float, default=DEFAULT_CUTOFF,
        help="heavy-atom contact distance in Å (default: %(default)s)",
    )
    args = ap.parse_args(argv)

    result = interface_contacts(
        args.input_filepath,
        args.output_folder,
        epitope_chains=_chain_list(args.epitope_chains),
        paratope_chains=_chain_list(args.paratope_chains),
        cutoff=args.cutoff,
    )
    print(
        f"{Path(args.input_filepath).name}: {result['n_contact_atom_pairs']} contact "
        f"atom-pairs @ {args.cutoff} Å; epitope footprint "
        f"{len(result['epitope_footprint'])} residues -> {result['_output']}"
    )


if __name__ == "__main__":
    main()
