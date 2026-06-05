# suggest-disulphides

Propose residue pairs that, if both mutated to cysteine, could form a viable
engineered disulphide. Unlike a plain Cβ–Cβ distance screen, it **models the Sγ
atoms** (three χ1 rotamers per residue, ideal Cys geometry) and checks the true
engineering criterion — **Sγ–Sγ ≈ 2.03 Å, χ3 (Cβ–Sγ–Sγ–Cβ) ≈ ±90°,
Cβ–Sγ–Sγ angle ≈ 105°** — keeping the best rotamer pair. Input may be `.pdb` or
`.cif`.

A pseudo-Cβ is built from the backbone where no Cβ exists (Gly), so any position
is scorable. Per-residue pLDDT (B-factor column) is reported, and an interface
footprint can be excluded so suggestions don't disrupt a binding surface.

## Contract

- Args: `input_filepath` and `output_folder`.
- Writes `<stem>_disulphides.json`: the criteria used and a ranked `candidates`
  list (residue identities + pLDDT, `cb_cb`, `ca_ca`, modelled `ss`, `chi3`,
  `cb_sg_sg_angles`, `score`, `low_plddt`).

## Usage

```bash
uv run suggest-disulphides <input.pdb|.cif> <output_folder> \
    [--group-a A:124-141] [--group-b A:324-419] \
    [--cb-min 3.0 --cb-max 4.5 --ca-max 7.5 --min-seq-sep 3] \
    [--exclude-footprint contacts.json] [--min-plddt 70]
```

- `--group-a`/`--group-b`: `CHAIN[:LO-HI][,...]` selections. With both set, only
  inter-group bridges are reported (e.g. a linker × a domain). With one/none, it
  scans all-vs-all within the selection.
- `--exclude-footprint`: an `interface-contacts` JSON; pairs touching its
  `epitope_footprint` are dropped (e.g. spare the W6/32 epitope).
- `--min-plddt`: flags pairs whose anchor confidence is below the threshold
  (useful when scanning predicted structures — a flexible linker has low pLDDT).

## Standalone use elsewhere

Self-contained folder. Copy it out and `uv run suggest-disulphides ...` or
`pip install .`. Runtime deps: BioPython, NumPy.
