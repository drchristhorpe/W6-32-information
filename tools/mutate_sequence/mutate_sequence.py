"""mutate_sequence -- apply precise point mutations to a sequence.

Self-contained tool (own helpers + CLI) so it can be lifted into a standalone
library or skill file unchanged. Stdlib only.

Input is a JSON file:

    {
      "sequence": "SLLQ...",
      "mutations": [
        {"position": 136, "from": "G", "to": "C"},
        {"position": 325, "from": "A", "to": "C"}
      ]
    }

Positions are 1-based. Every mutation is validated -- `from` must match the
residue actually present at `position` -- so there are no silent edits; a
mismatch raises an error listing the offending positions.

Contract: takes an input filepath and an output folder; writes
`<stem>_mutated.json` (original + mutated sequence, applied mutations) and
`<stem>_mutated.fasta` into the output folder.

Example:
    python -m mutate_sequence mutations.json out/
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


class MutationError(ValueError):
    """Raised when a requested mutation does not match the input sequence."""


def mutate_sequence(input_filepath: str | Path, output_folder: str | Path) -> dict:
    input_filepath = Path(input_filepath)
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    spec = json.loads(input_filepath.read_text())
    sequence = spec["sequence"]
    mutations = spec.get("mutations", [])
    if not sequence:
        raise MutationError("Input 'sequence' is empty.")

    seq = list(sequence)
    applied, errors = [], []
    for m in mutations:
        pos, frm, to = m["position"], m["from"], m["to"]
        if not (1 <= pos <= len(seq)):
            errors.append(f"position {pos} out of range 1..{len(seq)}")
            continue
        actual = seq[pos - 1]
        if actual != frm:
            errors.append(f"position {pos}: expected '{frm}' but sequence has '{actual}'")
            continue
        seq[pos - 1] = to
        applied.append({"position": pos, "from": frm, "to": to})

    if errors:
        raise MutationError("Mutation validation failed: " + "; ".join(errors))

    mutated = "".join(seq)
    label = ",".join(f"{m['from']}{m['position']}{m['to']}" for m in applied) or "none"
    result = {
        "input": str(input_filepath),
        "length": len(mutated),
        "n_mutations": len(applied),
        "mutations": applied,
        "original_sequence": sequence,
        "mutated_sequence": mutated,
    }

    out_json = output_folder / f"{input_filepath.stem}_mutated.json"
    out_json.write_text(json.dumps(result, indent=2))
    out_fasta = output_folder / f"{input_filepath.stem}_mutated.fasta"
    out_fasta.write_text(f">{input_filepath.stem} mutated {label}\n{mutated}\n")

    result["_output"] = str(out_json)
    result["_fasta"] = str(out_fasta)
    return result


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("input_filepath", help="JSON with 'sequence' and 'mutations'")
    ap.add_argument("output_folder", help="folder for <stem>_mutated.json/.fasta")
    args = ap.parse_args(argv)

    result = mutate_sequence(args.input_filepath, args.output_folder)
    label = ", ".join(f"{m['from']}{m['position']}{m['to']}" for m in result["mutations"])
    print(
        f"{Path(args.input_filepath).name}: applied {result['n_mutations']} mutation(s) "
        f"[{label}] to a {result['length']}-residue sequence -> {result['_output']}"
    )


if __name__ == "__main__":
    main()
