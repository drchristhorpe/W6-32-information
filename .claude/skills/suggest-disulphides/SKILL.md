---
name: suggest-disulphides
description: Propose engineerable disulphide bonds in a structure - residue pairs that, if mutated to cysteine, could form a viable S-S bond. Models Sγ rotamers and checks Sγ-Sγ ≈ 2.03 Å and χ3 ≈ ±90° (Disulfide-by-Design criterion), not just Cβ-Cβ distance. Use when asked where to add a stabilising disulphide / staple, e.g. between the β2m→α1 linker and the α3 domain of an MHC single-chain trimer. Supports restricting to two residue groups, excluding a protected interface footprint, and reporting per-residue pLDDT. Outputs a ranked candidates JSON.
---

# suggest-disulphides

Wraps the self-contained `suggest_disulphides` tool (`tools/suggest_disulphides/`).
For residue pairs in the selection it asks: if both were mutated to cysteine,
could they form a viable disulphide? It models Sγ atoms (three χ1 rotamers, ideal
Cys geometry) and keeps the best rotamer pair satisfying **Sγ-Sγ ≈ 2.03 Å, χ3 ≈
±90°, Cβ-Sγ-Sγ angle ≈ 105°** — the real engineering criterion, stronger than a
Cβ-Cβ distance screen. Input `.pdb` or `.cif`.

## When to use

- Find where to add a stabilising disulphide / staple, e.g. β2m→α1 linker × α3
  in an MHC single-chain trimer (design question 4c).
- Validate a proposed Cys/Cys pair's geometry.

## How to run

```bash
uv run suggest-disulphides <input.pdb|.cif> <output_folder> \
    [--group-a A:124-141] [--group-b A:324-419] \
    [--cb-min 3.0 --cb-max 4.5 --ca-max 7.5 --min-seq-sep 3] \
    [--exclude-footprint interface_description/<...>_contacts.json] \
    [--min-plddt 70]
```

- `--group-a`/`--group-b`: `CHAIN[:LO-HI]` selections; both set → only A×B
  bridges; one/none → all-vs-all within the selection.
- `--exclude-footprint`: drop pairs touching an `interface-contacts` footprint
  (e.g. don't disrupt the W6/32 epitope).
- `--min-plddt`: flag low-confidence anchors (predicted structures).

## Output contract

`<stem>_disulphides.json`: `criteria`, `n_candidates`, `n_excluded_by_footprint`,
and ranked `candidates` (each: `res_a`/`res_b` with pLDDT, `cb_cb`, `ca_ca`,
modelled `ss`, `chi3`, `cb_sg_sg_angles`, `score`, `low_plddt`).

## Notes

- On a predicted structure, a candidate whose partner sits in a low-pLDDT region
  (e.g. a flexible GGGGS linker) is geometrically suggestive but unreliable —
  treat the high-pLDDT anchor as the firm target and confirm the partner by
  explicit side-chain repacking / re-folding / MD.
- The ad-hoc `analysis/q4c_disulphide_scan.py` (Cβ-Cβ only) is kept for the
  record; this tool is the rigorous, reusable version.
