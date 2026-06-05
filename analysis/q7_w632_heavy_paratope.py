"""Job 7 -- the antibody-side view of the W6/32 heavy-chain variable region (VH).

Two patches on VH (heavy chain residues <= 120 of the Fab):
  (a) ANTIGEN contact patch -- VH residues contacting HLA + β2m (the paratope;
      keep these in a nanobody).
  (b) FORMER LIGHT-CHAIN patch now solvent-exposed -- VH residues that contact the
      light chain in the Fab but become solvent-exposed once the light chain is
      removed (the VH→VHH interface to re-engineer / camouflage).

Uses the `interface_contacts` tool for both contact sets, and BioPython
ShrakeRupley SASA on the isolated heavy chain to find which former-VL contacts are
exposed. Relative SASA (RSA) uses Tien et al. (2013) theoretical max-ASA.

Run from the experiment root:  uv run python analysis/q7_w632_heavy_paratope.py
"""

from __future__ import annotations

import json
from pathlib import Path

from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.SeqUtils import seq1

from interface_contacts import interface_contacts

ROOT = Path(__file__).resolve().parent.parent
COMPLEX = ROOT / "structures" / "pdb" / "hla_b_27_05__FRYNGLIHR__w632__1.pdb"
OUT = ROOT / "analysis" / "results" / "q7"
VH_MAX = 120        # heavy-chain variable domain ends ~120; CH1 follows
RSA_EXPOSED = 50.0  # % relative SASA to call a residue "totally solvent exposed"

# Tien et al. 2013 theoretical max ASA (Å²) for relative SASA.
MAX_ASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLU": 223.0, "GLN": 225.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}


def heavy_chain_rsa() -> dict:
    """Per-residue relative SASA (%) of the heavy chain in isolation."""
    model = PDBParser(QUIET=True).get_structure("cx", str(COMPLEX))[0]
    for cid in [c.id for c in list(model)]:
        if cid != "H":
            model.detach_child(cid)
    ShrakeRupley().compute(model, level="R")
    rsa = {}
    for r in model["H"]:
        if r.id[0] == " " and r.resname in MAX_ASA:
            rsa[r.id[1]] = round(r.sasa / MAX_ASA[r.resname] * 100, 1)
    return rsa


def vh_sequence() -> tuple[str, int, int]:
    """One-letter VH sequence (chain H residues <= VH_MAX) and its resseq span."""
    model = PDBParser(QUIET=True).get_structure("cx", str(COMPLEX))[0]
    res = [r for r in model["H"] if r.id[0] == " " and r.id[1] <= VH_MAX and "CA" in r]
    seq = "".join(seq1(r.resname, undef_code="X") for r in res)
    return seq, res[0].id[1], res[-1].id[1]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    # VH sequence for predictions (fold the isolated VH / design a nanobody)
    vh_seq, vh_lo, vh_hi = vh_sequence()
    (OUT / "w632_VH.fasta").write_text(
        f">w632_heavy_VH residues_{vh_lo}-{vh_hi}\n{vh_seq}\n"
    )

    # (a) VH residues contacting the antigen (HLA + β2m)
    ag = interface_contacts(COMPLEX, OUT / "antigen", epitope_chains=["A", "B"], paratope_chains=["H"])
    vh_antigen = sorted(
        (r for r in ag["paratope_footprint"] if r["resseq"] <= VH_MAX),
        key=lambda r: r["resseq"],
    )

    # (b) VH residues contacting the light chain
    vl = interface_contacts(COMPLEX, OUT / "lightchain", epitope_chains=["H"], paratope_chains=["L"])
    vh_vl = [r for r in vl["epitope_footprint"] if r["resseq"] <= VH_MAX]

    rsa = heavy_chain_rsa()
    for r in vh_vl:
        r["rsa_heavy_alone"] = rsa.get(r["resseq"])
    vh_vl.sort(key=lambda r: (r["rsa_heavy_alone"] is None, -(r["rsa_heavy_alone"] or 0)))
    exposed = [r for r in vh_vl if (r["rsa_heavy_alone"] or 0) >= RSA_EXPOSED]

    result = {
        "complex": str(COMPLEX),
        "vh_definition": f"chain H residues <= {VH_MAX}",
        "vh_sequence": vh_seq,
        "vh_span": [vh_lo, vh_hi],
        "rsa_exposed_threshold": RSA_EXPOSED,
        "antigen_contact_patch": vh_antigen,
        "antigen_contact_positions": [r["resseq"] for r in vh_antigen],
        "former_vl_contacts": vh_vl,
        "former_vl_exposed_positions": [r["resseq"] for r in exposed],
    }
    (OUT / "heavy_vh_patches.json").write_text(json.dumps(result, indent=2))

    # Simple companion file: just the residue numbers.
    summary = {
        "antigen_paratope": [r["resseq"] for r in vh_antigen],
        "light_chain_contacts": sorted(r["resseq"] for r in vh_vl),
    }
    (OUT / "vh_patches_simple.json").write_text(json.dumps(summary, indent=2))

    print(f"W6/32 heavy-chain VH analysis  ({OUT / 'heavy_vh_patches.json'})\n")
    print(f"VH sequence (res {vh_lo}-{vh_hi}, {len(vh_seq)} aa) -> {OUT / 'w632_VH.fasta'}:")
    print(f"    {vh_seq}\n")
    print(f"(a) ANTIGEN contact patch on VH ({len(vh_antigen)} residues, keep in nanobody):")
    print("    " + ", ".join(f"{r['resname']}{r['resseq']}" for r in vh_antigen))
    print(f"\n(b) Former light-chain contacts on VH, by exposure when L removed "
          f"({len(vh_vl)} residues; {len(exposed)} totally exposed at RSA>={RSA_EXPOSED}%):")
    for r in vh_vl:
        flag = "  <- EXPOSED (engineer)" if (r["rsa_heavy_alone"] or 0) >= RSA_EXPOSED else ""
        print(f"    {r['resname']}{r['resseq']}  RSA(H-alone)={r['rsa_heavy_alone']}%{flag}")


if __name__ == "__main__":
    main()
