---
name: mutate-sequence
description: Apply precise, validated point mutations to a protein/nucleotide sequence from a JSON spec ({sequence, mutations:[{position, from, to}]}, 1-based). Each mutation's `from` is checked against the actual residue, so there are no silent edits. Use when asked to introduce specific mutations into a sequence (e.g. engineering cysteines for a disulphide, G136C/A325C in the A*02:01/PRAME SCT) and get the mutated sequence out as JSON + FASTA.
---

# mutate-sequence

Wraps the self-contained `mutate_sequence` tool (`tools/mutate_sequence/`). Applies
point mutations from a JSON spec, validating each `from` residue against the
sequence (1-based positions) so mistakes surface as errors rather than silent
edits. No dependencies.

## When to use

- Introduce specific point mutations (e.g. engineer Cys for a disulphide such as
  the 4c candidates) and emit the mutated sequence.
- Any precise sequence edit where you want the `from` residues checked.

## How to run

```bash
uv run mutate-sequence <mutations.json> <output_folder>
```

Input JSON:
```json
{
  "sequence": "SLLQ...",
  "mutations": [
    {"position": 136, "from": "G", "to": "C"},
    {"position": 325, "from": "A", "to": "C"}
  ]
}
```

## Output contract

Into `<output_folder>`: `<stem>_mutated.json` (`original_sequence`,
`mutated_sequence`, applied `mutations`, `length`) and `<stem>_mutated.fasta`.

## Notes

- Positions are 1-based; a `from`/actual mismatch or out-of-range position aborts
  with an error listing the offenders.
- The A*02:01/PRAME SCT sequence (419 aa) can be taken from a prediction's
  `run.json`/`job_request.json`; e.g. positions 136 (G) and 325 (A) are valid C targets.
