# Conclusions

Findings for the design questions in [EXPERIMENTAL_DESIGN.md](EXPERIMENTAL_DESIGN.md).
This is the reader-facing answers document; the process/lab notebook is in
[CHANGELOG.md](CHANGELOG.md). Each section records the question, the method
(which tools/skills), the result, and the design implication.

Status: 🟢 complete — all three tools built; Aims 1–4 and design questions 4a, 4b,
and 4c (extension) answered.

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

## The W6/32 epitope across alleles (Aim 2)

**Question:** How does the W6/32 interface differ between HLA-B\*27:05 and
HLA-A\*02:01 (the PRAME structure)?

**Method:** `compare-structures`, A\*02:01/PRAME experimental vs B\*27:05 W6/32
complex, Cα, in-frame, restricted to the W6/32 footprint. Result in
`analysis/results/aim2/`.

**Result:** footprint backbone differs by only **1.36 Å** (19 residues); the
largest differences are again the α3 222–228 loop (Asp227, Gln226) and the β2m
N-terminus (Ile1, Arg3).

**Conclusion:** the W6/32 epitope surface is **well conserved between the two
alleles** — consistent with W6/32 being pan-class-I — with the modest differences
localised to the same plastic α3 loop / β2m N-terminus seen in 4a. (This pair is
A\*02:01-free vs B\*27:05-bound, so allele and binding-state effects are conflated;
the conservation is the robust takeaway.)

---

## 4a — Does W6/32 binding change the B\*27:05 conformation? (interface only)

**Question:** How much does W6/32 change the B\*27:05 conformation on binding,
looking only at the W6/32 interface?

**Method:** `compare-structures` (skill), apo `hla_b_27_05__RRFSRSPIRR` vs the
W6/32-bound complex, Cα, restricted to the W6/32 footprint (peptide excluded).
Both experimental structures are pre-framed, so deviation is measured **in-frame**
(total change vs the frame) and **superposed** (internal distortion). Both
crystallographic copies + a specificity control. Driver:
`analysis/q4a_conformational_change.py`; results in `analysis/results/q4a/`.

**Result (consistent across both copies):**
| measure | copy 1 | copy 2 |
|---|---|---|
| footprint RMSD, in-frame | 2.22 Å | 2.13 Å |
| footprint RMSD, superposed | 1.23 Å | 1.16 Å |

- Largest movers are the **α3 contact loop (Glu222–Glu229, peaking at Gln226,
  3.9–4.4 Å)** — exactly the residues that hydrogen-bond/salt-bridge to the W6/32
  heavy chain — plus the **β2m N-terminus (Ile1, ~3.5 Å)**.
- **Specificity control (decisive):** over the whole HLA heavy chain, the α1/α2
  peptide-binding platform (the framing region) differs by only **0.48 Å**, while
  α3 differs by **2.78 Å**. The change is α3-specific, not a global crystal/peptide
  offset — and confirms the two structures are genuinely co-framed on the groove.

**Conclusion:** W6/32 engagement is associated with a substantial, **localised
reorientation of the α3 domain** (~2.8 Å in-frame; ~1.2 Å residual internal
distortion after footprint fit) and a β2m N-terminal shift, while the α1/α2
peptide-binding groove is essentially unchanged (0.48 Å).

**Caveat:** single apo/bound pair with different peptides (RRFSRSPIRR vs
FRYNGLIHR) and crystal forms; the α1/α2 control removes the global-offset confound
but the precise magnitude remains indicative rather than definitive.

**Implication:** designs must let the α3 domain adopt its bound conformation
(α3 is mobile relative to the rigid α1/α2 platform); preserving the α3 222–229
loop and β2m N-terminus is key to W6/32 recognition.

## 4b — Is the Boltz A\*02:01/PRAME SCT good enough as a design input?

**Question:** Can we trust the Boltz-2 prediction of the HLA-A\*02:01/PRAME
single-chain trimer as a protein-design input? Is it close to experiment, good in
general, is the peptide conformation right, and is the W6/32 interface well modelled?

**Method:** `align-structure` then `compare-structures` (skills): aligned Boltz SCT
vs experimental `hla_a_02_01__SLLQHLIGL`, Cα, in-frame, at three scopes. The
W6/32-interface scope reuses the B\*27:05 footprint (α3 numbering conserved, β2m
identical). Driver `analysis/q4b_boltz_sct_quality.py`; results in
`analysis/results/q4b/`.

**Result:**
| scope | n (Cα) | global RMSD | notes |
|---|---|---|---|
| overall HLA + β2m | 339 | **1.01 Å** | design-grade overall fold |
| **peptide (SLLQHLIGL)** | 9 | **1.00 Å** | correct presented pose; largest at His5/Leu9 |
| W6/32 interface (footprint) | 19/22 | 1.56 Å | β2m + α2 patch ~1 Å; **α3 loop off** |

- Boltz self-confidence is high: complex pLDDT **0.90**, pTM **0.93**.
- The interface error is **localised to the α3 222–228 loop** (Asp227 4.4 Å,
  Gln226 3.9 Å); β2m and the α2 patch (Gly120/Lys121) are modelled to ~1 Å.
- **Convergence with 4a:** the α3 222–229 loop that is mispredicted here is exactly
  the loop that (i) contacts W6/32 and (ii) reorients ~2.8 Å on W6/32 binding (4a).
  A plastic, antibody-engaged loop is both hard to predict and conformation-dependent.

**Conclusion:** The Boltz SCT is a **sound design input for the overall fold, the
peptide-binding groove, and the PRAME peptide conformation** (all ~1 Å, with high
model confidence). The W6/32 epitope is **mostly reliable (β2m, α2) but not at the
α3 222–228 contact loop** (~4 Å off). Because that loop is plastic and the
prediction is of the free SCT, its predicted α3 conformation should not be trusted
for W6/32-interface design.

**Caveat:** the experimental A\*02:01/PRAME is also antibody-free, so the α3
discrepancy is free-prediction vs free-experiment in a flexible loop; for
W6/32-bound design geometry, take the α3 conformation from the B\*27:05 W6/32
complex rather than from either free A\*02:01 structure.

**Implication:** safe to design on the groove/peptide; treat the α3 loop as
flexible / graft the bound-state α3 conformation when W6/32 binding must be retained.

## 4c — Best location for an α3 ↔ β2m-linker disulphide (extension)

**Question:** If we add a disulphide between the MHC α3 domain and the β2m→α1
linker, where is the best location?

**Method:** geometric scan (`analysis/q4c_disulphide_scan.py`) on the Boltz
A\*02:01/PRAME SCT (the only format with the linker). For each linker (res
**124–141**) × α3 (res 324–419) pair, build a pseudo-Cβ from the backbone (the
linker is poly-Gly/Ser; we ask "if mutated to Cys") and keep pairs with Cβ–Cβ
3.0–4.5 Å and Cα–Cα 4.0–7.5 Å, ranked toward the ideal Cβ–Cβ ≈ 3.85 Å. pLDDT
recorded per residue. Results in `analysis/results/q4c/`.

> Linker boundary: native α1 begins G-S-HSMRYF, so SCT residue 142 (Gly) is α1,
> not linker — the linker ends at 141.

**Result — 6 candidate pairs (best first):**
| linker (Cys) | α3 (Cys) | Cβ–Cβ | Cα–Cα | α3 pLDDT |
|---|---|---|---|---|
| Gly135 | **Leu407** | 4.04 Å | 4.85 Å | 91 |
| Gly132 | **Lys327** | 3.34 Å | 4.69 Å | 91 |
| Gly125 | Trp345 | 3.64 Å | 6.42 Å | 93 |
| Gly134 | Ala325 | 4.41 Å | 5.70 Å | 91 |
| Ser133 | Thr328 | 4.45 Å | 4.98 Å | 93 |
| Gly126 | Leu347 | 4.49 Å | 6.72 Å | 92 |

- **Top picks:** Gly135↔Leu407 (best-balanced geometry, Cα–Cα 4.85 Å) and
  Gly132↔Lys327 (tightest Cβ–Cβ; Lys→Cys also removes a surface charge).
- All candidate α3 anchors are **outside the W6/32 footprint** (which maps to SCT
  ~363–384), so the staple should not disrupt antibody binding.

**Conclusion:** the most promising sites pair an α3 anchor (**Leu407** or
**Lys327**, both confidently placed and clear of the W6/32 epitope) with a nearby
linker glycine (Gly135 or Gly132 respectively) mutated to cysteine.

**Major caveat:** the linker is intrinsically disordered — its predicted pLDDT is
**~36** vs ~92 for α3 — so the *exact* linker partner and the bond geometry are
unreliable from a single static prediction. Treat the **α3 anchor as the reliable
design target**; then choose/optimise the linker Cys by explicit SG-rotamer
modelling and re-folding/MD, since a flexible linker can sample many
conformations to satisfy the bond.
