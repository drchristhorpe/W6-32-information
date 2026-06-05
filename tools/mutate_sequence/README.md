# mutate-sequence

Apply precise point mutations to a sequence from a JSON spec. Every mutation is
**validated** — the `from` residue must match what is actually at `position`
(1-based) — so there are no silent edits; any mismatch is an error.

No dependencies (stdlib only).

## Contract

- Args: `input_filepath` (JSON) and `output_folder`.
- Input JSON:
  ```json
  {
    "sequence": "SLLQ...",
    "mutations": [
      {"position": 136, "from": "G", "to": "C"},
      {"position": 325, "from": "A", "to": "C"}
    ]
  }
  ```
- Writes `<stem>_mutated.json` (original + mutated sequence, applied mutations)
  and `<stem>_mutated.fasta`.

## Usage

```bash
uv run mutate-sequence <mutations.json> <output_folder>
```

## Standalone use elsewhere

Self-contained folder with no runtime dependencies. Copy it out and
`uv run mutate-sequence ...` or `pip install .` (entry point `mutate-sequence`).
