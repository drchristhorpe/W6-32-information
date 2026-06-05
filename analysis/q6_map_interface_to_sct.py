"""Job 6 -- map the W6/32 interface onto the A*02:01/PRAME single-chain trimer.

Produces two lists, in SCT residue numbering:
  1. positions IN the W6/32 interface (the footprint), and
  2. positions in/near the interface that MOVE on binding (apo -> W6/32-bound).

The footprint and the apo->bound movements were determined on the B*27:05 complex
in mature HLA / β2m numbering. They are mapped onto the SCT (a single fused chain
with its own numbering) by local sequence alignment of the B*27:05 HLA and β2m
chains against the SCT chain -- robust to the numbering offset and the B*27↔A*02
allele differences.

Run from the experiment root:  uv run python analysis/q6_map_interface_to_sct.py
"""

from __future__ import annotations

import json
from pathlib import Path

from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

ROOT = Path(__file__).resolve().parent.parent
COMPLEX = ROOT / "structures" / "pdb" / "hla_b_27_05__FRYNGLIHR__w632__1.pdb"
SCT = ROOT / "sct_predictions" / "aligned" / "boltz" / "hla_a_02_01__single_chain_trimer__SLLQHLIGL" / "sample_0_predicted_structure_aligned.pdb"
FOOTPRINT = ROOT / "interface_description" / "hla_b_27_05__FRYNGLIHR__w632__1_contacts.json"
MOVERS = ROOT / "analysis" / "results" / "q4a" / "inframe" / "hla_b_27_05__RRFSRSPIRR__vs__hla_b_27_05__FRYNGLIHR__w632__1_compare.json"
OUT = ROOT / "analysis" / "results" / "q6"
MOVE_THRESHOLD = 1.0  # Å, in-frame apo->bound displacement to count as "moves on binding"


def _residues(chain):
    return [r for r in chain if r.id[0] == " " and "CA" in r]


def _seq(residues):
    return "".join(seq1(r.resname, undef_code="X") for r in residues)


def _build_map(src_residues, sct_residues):
    """{src_resseq -> sct_resseq} via local sequence alignment."""
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.mode = "local"
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aln = aligner.align(_seq(src_residues), _seq(sct_residues))[0]
    mapping = {}
    for (s0, s1), (t0, t1) in zip(*aln.aligned):
        for si, ti in zip(range(s0, s1), range(t0, t1)):
            mapping[src_residues[si].id[1]] = sct_residues[ti].id[1]
    return mapping


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    parser = PDBParser(QUIET=True)
    cx = parser.get_structure("cx", str(COMPLEX))[0]
    sct_chain = parser.get_structure("sct", str(SCT))[0]["A"]
    sct_res = _residues(sct_chain)
    sct_by_pos = {r.id[1]: r.resname for r in sct_res}

    # chain A = HLA heavy, chain B = β2m in the complex
    map_by_chain = {
        "A": _build_map(_residues(cx["A"]), sct_res),
        "B": _build_map(_residues(cx["B"]), sct_res),
    }

    def to_sct(chain, resseq):
        return map_by_chain.get(chain, {}).get(resseq)

    # ---- List 1: positions in the interface (footprint) ----
    footprint = json.load(open(FOOTPRINT))["epitope_footprint"]
    in_interface = []
    for r in footprint:
        sct_pos = to_sct(r["chain"], r["resseq"])
        if sct_pos is None:
            continue
        in_interface.append({
            "sct_position": sct_pos,
            "sct_resname": sct_by_pos.get(sct_pos),
            "source": f"{r['resname']}{r['resseq']}/{r['chain']}",
            "region": "b2m" if r["chain"] == "B" else "HLA",
            "min_distance": r["min_distance"],
        })
    in_interface.sort(key=lambda e: e["sct_position"])

    # ---- List 2: positions in/near the interface that move on binding ----
    movers = []
    for e in json.load(open(MOVERS))["per_residue"]:
        if e["rmsd"] < MOVE_THRESHOLD:
            continue
        sct_pos = to_sct(e["ref_chain"], e["resseq"])
        if sct_pos is None:
            continue
        movers.append({
            "sct_position": sct_pos,
            "sct_resname": sct_by_pos.get(sct_pos),
            "source": f"{e['resname']}{e['resseq']}/{e['ref_chain']}",
            "region": "b2m" if e["ref_chain"] == "B" else "HLA",
            "displacement": e["rmsd"],
        })
    movers.sort(key=lambda e: e["displacement"], reverse=True)

    result = {
        "sct": str(SCT),
        "footprint_source": str(FOOTPRINT),
        "movers_source": str(MOVERS),
        "move_threshold_A": MOVE_THRESHOLD,
        "interface_positions_sct": [e["sct_position"] for e in in_interface],
        "moving_positions_sct": [e["sct_position"] for e in movers],
        "interface_detail": in_interface,
        "moving_detail": movers,
    }
    out = OUT / "interface_on_sct.json"
    out.write_text(json.dumps(result, indent=2))

    print(f"W6/32 interface mapped onto the A*02:01/PRAME SCT  ({out})\n")
    print(f"List 1 - IN the interface ({len(in_interface)} positions, SCT numbering):")
    print("  " + ", ".join(f"{e['sct_resname']}{e['sct_position']}" for e in in_interface))
    print(f"\nList 2 - move on binding (>= {MOVE_THRESHOLD} Å, {len(movers)} positions):")
    for e in movers:
        print(f"  {e['sct_resname']}{e['sct_position']} (SCT)  <- {e['source']}  {e['displacement']:.2f} Å  [{e['region']}]")


if __name__ == "__main__":
    main()
