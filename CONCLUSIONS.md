# Conclusions

Findings for the design questions in [EXPERIMENTAL_DESIGN.md](EXPERIMENTAL_DESIGN.md).
This is the reader-facing answers document; the process/lab notebook is in
[CHANGELOG.md](CHANGELOG.md). Each section records the question, the method
(which tools/skills), the result, and the design implication.

Status: 🟢 complete — four composable tools built; Aims 1–4 and extension tasks
4a, 4b, 4c, and 4d answered.

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

**Method:** two scans on the Boltz A\*02:01/PRAME SCT (the only format with the
linker), over linker (res **124–141**) × α3 (res 324–419):
1. **Proximity scan** (`analysis/q4c_disulphide_scan.py`) — pseudo-Cβ from
   backbone, keep pairs with Cβ–Cβ 3.0–4.5 Å (crude distance proxy).
2. **Rigorous scan** (`suggest-disulphides` tool/skill) — models Sγ over three χ1
   rotamers and requires real disulphide geometry: Sγ–Sγ ≈ 2.03 Å, χ3 ≈ ±90°,
   Cβ–Sγ–Sγ ≈ 105°; W6/32 footprint excluded; per-residue pLDDT.

Results in `analysis/results/q4c/`.

> Linker boundary: native α1 begins G-S-HSMRYF, so SCT residue 142 (Gly) is α1,
> not linker — the linker ends at 141.
>
> Tool validation: scanning the whole chain, `suggest-disulphides` recovers all
> three native disulphides (α2 Cys242–Cys305, α3 Cys344–Cys400, β2m Cys49–Cys104)
> at Sγ–Sγ 2.0–2.1 Å, χ3 −72…−88°, with no false positives.

**Result:**
- The proximity scan flags **6 linker×α3 pairs** by Cβ–Cβ alone — α3 anchors
  Leu407, Lys327, Trp345, Ala325, Thr328, Leu347 (all pLDDT ~91–93, all outside
  the W6/32 footprint at SCT ~363–384).
- The rigorous scan accepts **0** of them: under proper Sγ/χ3 geometry, **no
  linker×α3 pair forms a viable disulphide in this static model**. The 6 Cβ–Cβ
  hits are false positives.
- Cause: the linker is intrinsically disordered (pLDDT **~36** vs ~92 for α3), so
  its single predicted conformation places no residue in bonding geometry.

**Conclusion:** there is **no reliable disulphide site identifiable from the static
prediction** — the proper-geometry test rejects every proximity hit. The α3 face
nearest the linker is confidently placed (**Leu407 / Lys327**, clear of the W6/32
epitope) and is the sensible **anchor**, but the disordered linker means the
partner Cys cannot be sited from this model.

**Recommended next step (out of scope here):** fix an α3 anchor, then explicitly
sample/repack the linker — or re-fold/MD with a Cys pair restrained — to find a
linker position that reaches true Sγ–Sγ/χ3 geometry. This is exactly why the crude
Cβ–Cβ proxy is insufficient and the Sγ/χ3 tool matters.

## 4d — Which predictor is the best design starting point? (extension)

**Question:** Is the AlphaFold3, ESMFold2(-fast) or Boltz-2 A\*02:01/PRAME SCT
closer to the **bound** structure in the W6/32 interface region?

**Method:** `compare-structures`, each predictor's aligned PRAME SCT vs the
W6/32-bound B\*27:05 complex (the only bound reference — no bound A\*02:01 exists),
Cα, in-frame, restricted to the W6/32 footprint; α3 and β2m reported separately.
All scored against the same reference, so the conserved allele offset (Aim 2,
1.36 Å) is constant and the ranking is fair. Driver
`analysis/q4d_predictor_interface_to_bound.py`; results in `analysis/results/q4d/`.

**Result — RMSD to W6/32-bound B\*27:05 (Å):**
| predictor | footprint | α3 | β2m |
|---|---|---|---|
| **esmfold2-fast** | **1.40** | 1.66 | 0.73 |
| **boltz** | **1.40** | **1.61** | 0.93 |
| alphafold3 | 1.58 | 1.92 | 0.61 |
| esmfold2 | 1.68 | 2.04 | 0.66 |

- All four fall in a narrow band; **β2m is modelled well by everyone (0.6–0.9 Å)**,
  so the discriminator is the **plastic α3 contact loop**.
- **ESMFold2-fast and Boltz tie on the overall footprint (1.40 Å)** and both clearly
  beat AlphaFold3 and standard ESMFold2.
- On the most design-critical α3 region, **Boltz is best (1.61 Å)**, narrowly ahead
  of ESMFold2-fast (1.66). AlphaFold3 models β2m best but its α3 is among the worst.

**Conclusion:** **Boltz-2 and ESMFold2-fast are the best design starting points**
for retaining W6/32 binding; they are tied at the interface overall, with **Boltz
holding a slight edge at the α3 contact loop** that matters most. AlphaFold3 — and
plain ESMFold2 — are *not* better here despite AF3's general reputation. Note all
predictions still sit 1.6–2.0 Å from the bound α3 conformation (they are free-state
predictions of a plastic loop), so the α3 loop should still be refined toward the
bound state regardless of predictor.

**Caveat:** the spread is modest (~0.3 Å at the footprint) and includes the
constant allele offset; the ranking is robust but the absolute closeness should
not be over-read.

## Job 6 — W6/32 interface mapped onto the SCT

**Task:** map the W6/32 contact positions onto the A\*02:01/PRAME SCT; list the
positions **in the interface** and those **in/near the interface that move on
binding**.

**Method:** the footprint (interface-contacts) and apo→bound movements (4a) are in
B\*27:05/mature numbering; mapped onto SCT numbering by local sequence alignment of
the B\*27:05 HLA and β2m chains against the SCT chain (`analysis/q6_map_interface_to_sct.py`;
result `analysis/results/q6/interface_on_sct.json`). The mapping reproduces the
construct offsets (β2m +24, HLA +141) and the SCT residue identities match at every
mapped position — confirming both the mapping and that the epitope is conserved
B\*27→A\*02.

**List 1 — in the W6/32 interface (22 positions, SCT numbering):**
- β2m: Ile25, Arg27, Thr28, Lys30, Lys82, Asp83, Trp84, Ser85
- HLA α2 patch: Gly261, Lys262, Asp263
- HLA α3 loop: Glu363, Asp364, Thr366, Gln367, Asp368, Thr369, Glu370, Leu371,
  Glu373, Thr374, Lys384

**List 2 — in/near the interface and moving on binding (≥1.0 Å, 12 positions):**
- HLA α3: **Gln367 (4.35 Å)**, Asp368 (3.73), Glu363 (3.53), Asp364 (3.46),
  Thr366 (3.42), Thr369 (3.17), Glu370 (3.00), Leu371 (1.98), Glu373 (1.64),
  Lys384 (1.55), Thr374 (1.34)
- β2m: Ile25 (2.54 Å)

**Conclusion / design use:** in SCT numbering, the W6/32 epitope to preserve is the
β2m patch (25–30, 82–85), the α2 patch (261–263) and the α3 loop (363–374, 384).
The α3 loop **363–371** (plus β2m Ile25) is the conformationally responsive core —
the part that moves on binding and was hardest to predict (4b) — so it is both the
key recognition determinant and the region to handle as flexible / refine toward
the bound state in design.

---

*All aims (1–4), extension tasks (4a–4d) and build jobs (1–3, 5–6) are addressed
above; the five composable tools are listed in CLAUDE.md.*
