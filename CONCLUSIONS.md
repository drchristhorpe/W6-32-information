# Conclusions

Findings for the design questions in [EXPERIMENTAL_DESIGN.md](EXPERIMENTAL_DESIGN.md).
This is the reader-facing answers document; the process/lab notebook is in
[CHANGELOG.md](CHANGELOG.md). Each section records the question, the method
(which tools/skills), the result, and the design implication.

Status: 🟡 in progress — tooling complete through `interface_contacts`;
`compare_structures` and the 4a/4b/4c analyses pending.

---

## The W6/32 interface (Aim 1)

**Question:** What is the interface between the W6/32 heavy chain and the
HLA-B\*27:05 + β2m surface in the experimental structure?

**Method:** `interface-contacts` on `hla_b_27_05__FRYNGLIHR__w632` (heavy chain =
paratope, HLA heavy + β2m = epitope), 5.0 Å heavy-atom cutoff. Both
crystallographic copies analysed → `interface_description/`.

**Result (reproducible across both copies, 22–23 epitope residues):**
- **HLA heavy chain (α3 domain):** a dominant cluster at residues ~222–243, plus
  a small α2 patch at ~120–122.
- **β2m:** residues ~1–6 and ~58–61.
- Anchored by polar contacts including a Glu229–Arg69 salt bridge (~2.8 Å) and
  Asp122/Gln226 H-bonds to the heavy chain.

This is the classic conformational, β2m-dependent W6/32 epitope (α3 + β2m), which
validates the pipeline and defines the footprint to be preserved/targeted in
design.

**Implication:** design that must retain W6/32 binding should preserve the α3
(222–243) and β2m (1–6, 58–61) surface conformation.

---

## 4a — Does W6/32 binding change the B\*27:05 conformation? (interface only)

_Pending `compare_structures`._ Plan: compare W6/32-bound B\*27:05 vs apo
`hla_b_27_05__RRFSRSPIRR`, restricted to the footprint above (peptide ignored).

## 4b — Is the Boltz A\*02:01/PRAME SCT good enough as a design input?

_Pending `compare_structures`._ Plan: aligned Boltz SCT vs A\*02:01/PRAME
experimental at global / peptide-only / W6/32-interface scopes, cross-read with
Boltz confidence.

## 4c — Best location for an α3 ↔ β2m-linker disulphide (extension)

_Pending._ Plan: geometric Cβ–Cβ scan over α3 × β2m-linker residue pairs on the
aligned SCT.
