"""suggest_disulphides -- propose engineerable disulphide bonds in a structure.

Self-contained tool (own helpers + CLI) so it can be lifted into a standalone
library or skill file unchanged.

For residue pairs in the requested selection, asks: if both were mutated to
cysteine, could they form a viable disulphide? Unlike a plain Cβ-Cβ distance
screen, this models the Sγ atoms (three χ1 rotamers per residue, ideal Cys
geometry) and checks the true engineering criterion -- Sγ-Sγ ≈ 2.03 Å,
χ3 (Cβ-Sγ-Sγ-Cβ) ≈ ±90°, Cβ-Sγ-Sγ angle ≈ 105° -- keeping the best rotamer pair.

A pseudo-Cβ is built from the backbone when no Cβ exists (Gly) so any position is
scorable. Per-residue pLDDT (B-factor column) is reported; an interface footprint
can be excluded so suggestions don't disrupt a binding surface.

Contract: takes an input filepath and an output folder; writes
`<stem>_disulphides.json` into the output folder.

Example (β2m→α1 linker × α3, sparing the W6/32 epitope):
    python -m suggest_disulphides sct.pdb out/ \\
        --group-a A:124-141 --group-b A:324-419 \\
        --exclude-footprint interface_description/..._contacts.json
"""

from __future__ import annotations

import argparse
import json
from math import cos, radians, sin
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser

# Ideal cysteine / disulphide geometry.
CB_SG = 1.808              # Cβ-Sγ bond length (Å)
CA_CB_SG = radians(114.4)  # Cα-Cβ-Sγ angle
CHI1_ROTAMERS = (-60.0, 60.0, 180.0)  # N-Cα-Cβ-Sγ
IDEAL_SS = 2.03            # Sγ-Sγ (Å)
IDEAL_CHI3 = 90.0          # |Cβ-Sγ-Sγ-Cβ| (deg)


def _load(path: str | Path, name: str):
    path = Path(path)
    parser = MMCIFParser(QUIET=True) if path.suffix.lower() in (".cif", ".mmcif") else PDBParser(QUIET=True)
    return parser.get_structure(name, str(path))


def _pseudo_cb(n, ca, c):
    b, cvec = ca - n, c - ca
    a = np.cross(b, cvec)
    return -0.58273431 * a + 0.56802827 * b - 0.54067466 * cvec + ca


def _place(a, b, c, length, angle, dihedral):
    """NeRF: position of D with |C-D|=length, angle(B,C,D)=angle, dihedral(A,B,C,D)."""
    bc = c - b
    bc /= np.linalg.norm(bc)
    nrm = np.cross(b - a, bc)
    nrm /= np.linalg.norm(nrm)
    m = np.cross(nrm, bc)
    d = np.array([-length * cos(angle),
                  length * sin(angle) * cos(dihedral),
                  length * sin(angle) * sin(dihedral)])
    return c + d[0] * bc + d[1] * m + d[2] * nrm


def _dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))


def _angle(p0, p1, p2):
    a, b = p0 - p1, p2 - p1
    cosang = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))
    return np.degrees(np.arccos(np.clip(cosang, -1.0, 1.0)))


def _parse_selection(spec):
    """`A:124-141,B` -> {('A',(124,141)), ('B',None)}. None spec -> None (all)."""
    if not spec:
        return None
    out = []
    for part in spec.split(","):
        part = part.strip()
        if ":" in part:
            ch, rng = part.split(":")
            lo, hi = rng.split("-")
            out.append((ch, (int(lo), int(hi))))
        else:
            out.append((part, None))
    return out


def _in_selection(chain_id, resseq, selection):
    if selection is None:
        return True
    for ch, rng in selection:
        if ch == chain_id and (rng is None or rng[0] <= resseq <= rng[1]):
            return True
    return False


def _collect(model, selection):
    """Residues with backbone N,CA,C in the selection -> list of dicts."""
    res = []
    for chain in model:
        for r in chain:
            if r.id[0] != " " or not all(a in r for a in ("N", "CA", "C")):
                continue
            if not _in_selection(chain.id, r.id[1], selection):
                continue
            n, ca, c = r["N"].coord, r["CA"].coord, r["C"].coord
            cb = r["CB"].coord if "CB" in r else _pseudo_cb(n, ca, c)
            sg = [_place(n, ca, cb, CB_SG, CA_CB_SG, radians(chi)) for chi in CHI1_ROTAMERS]
            res.append({
                "chain": chain.id, "resseq": r.id[1], "resname": r.resname,
                "ca": ca, "cb": cb, "sg": sg,
                "plddt": round(float(np.mean([a.bfactor for a in r])), 1),
            })
    return res


def suggest_disulphides(
    input_filepath, output_folder,
    group_a=None, group_b=None,
    cb_cb=(3.0, 4.5), ca_ca_max=7.5, ss_window=(1.8, 2.3),
    chi3_window=(60.0, 120.0), cbss_window=(90.0, 120.0),
    min_seq_sep=3, exclude_footprint=None, min_plddt=None,
) -> dict:
    input_filepath, output_folder = Path(input_filepath), Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)
    model = _load(input_filepath, input_filepath.stem)[0]

    sel_a = _parse_selection(group_a)
    sel_b = _parse_selection(group_b)
    res_a = _collect(model, sel_a)
    res_b = _collect(model, sel_b) if sel_b is not None else res_a
    inter_group = sel_b is not None

    excluded = _footprint_set(exclude_footprint)
    candidates, n_excluded = [], 0

    for i, ri in enumerate(res_a):
        for jx, rj in enumerate(res_b):
            # avoid double-counting / self when scanning one group all-vs-all
            if not inter_group and jx <= i:
                continue
            if ri["chain"] == rj["chain"] and abs(ri["resseq"] - rj["resseq"]) < min_seq_sep:
                continue
            if np.linalg.norm(ri["cb"] - rj["cb"]) > cb_cb[1] or np.linalg.norm(ri["cb"] - rj["cb"]) < cb_cb[0]:
                continue
            if np.linalg.norm(ri["ca"] - rj["ca"]) > ca_ca_max:
                continue

            best = _best_rotamer_pair(ri, rj, ss_window, chi3_window, cbss_window)
            if best is None:
                continue
            if (ri["chain"], ri["resseq"]) in excluded or (rj["chain"], rj["resseq"]) in excluded:
                n_excluded += 1
                continue
            cand = {
                "res_a": {"chain": ri["chain"], "resseq": ri["resseq"], "resname": ri["resname"], "plddt": ri["plddt"]},
                "res_b": {"chain": rj["chain"], "resseq": rj["resseq"], "resname": rj["resname"], "plddt": rj["plddt"]},
                "cb_cb": round(float(np.linalg.norm(ri["cb"] - rj["cb"])), 2),
                "ca_ca": round(float(np.linalg.norm(ri["ca"] - rj["ca"])), 2),
                **best,
                "low_plddt": bool(min_plddt and (ri["plddt"] < min_plddt or rj["plddt"] < min_plddt)),
            }
            candidates.append(cand)

    candidates.sort(key=lambda c: c["score"])
    result = {
        "input": str(input_filepath),
        "group_a": group_a, "group_b": group_b,
        "criteria": {
            "cb_cb": list(cb_cb), "ca_ca_max": ca_ca_max, "ss_window": list(ss_window),
            "chi3_window": list(chi3_window), "cbss_window": list(cbss_window),
            "min_seq_sep": min_seq_sep, "ideal_ss": IDEAL_SS, "ideal_chi3": IDEAL_CHI3,
        },
        "exclude_footprint": str(exclude_footprint) if exclude_footprint else None,
        "n_excluded_by_footprint": n_excluded,
        "n_candidates": len(candidates),
        "candidates": candidates,
    }
    out = output_folder / f"{input_filepath.stem}_disulphides.json"
    with open(out, "w") as fh:
        json.dump(result, fh, indent=2)
    result["_output"] = str(out)
    return result


def _best_rotamer_pair(ri, rj, ss_window, chi3_window, cbss_window):
    """Best χ1 rotamer combination satisfying the disulphide windows, or None."""
    best = None
    for sgi in ri["sg"]:
        for sgj in rj["sg"]:
            ss = float(np.linalg.norm(sgi - sgj))
            if not (ss_window[0] <= ss <= ss_window[1]):
                continue
            chi3 = float(_dihedral(ri["cb"], sgi, sgj, rj["cb"]))
            if not (chi3_window[0] <= abs(chi3) <= chi3_window[1]):
                continue
            ang_i = float(_angle(ri["cb"], sgi, sgj))
            ang_j = float(_angle(rj["cb"], sgj, sgi))
            if not (cbss_window[0] <= ang_i <= cbss_window[1] and cbss_window[0] <= ang_j <= cbss_window[1]):
                continue
            score = abs(ss - IDEAL_SS) + abs(abs(chi3) - IDEAL_CHI3) / 90.0
            if best is None or score < best["score"]:
                best = {
                    "ss": round(ss, 2), "chi3": round(chi3, 1),
                    "cb_sg_sg_angles": [round(ang_i, 1), round(ang_j, 1)],
                    "score": round(score, 3),
                }
    return best


def _footprint_set(footprint):
    if not footprint:
        return set()
    data = json.load(open(footprint))
    return {(r["chain"], r["resseq"]) for r in data.get("epitope_footprint", [])}


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("input_filepath", help="structure (.pdb or .cif)")
    ap.add_argument("output_folder")
    ap.add_argument("--group-a", default=None, help="CHAIN[:LO-HI][,...] (default: whole structure)")
    ap.add_argument("--group-b", default=None, help="CHAIN[:LO-HI][,...]; if set, only A×B bridges")
    ap.add_argument("--cb-min", type=float, default=3.0)
    ap.add_argument("--cb-max", type=float, default=4.5)
    ap.add_argument("--ca-max", type=float, default=7.5)
    ap.add_argument("--min-seq-sep", type=int, default=3)
    ap.add_argument("--exclude-footprint", default=None, help="contacts.json; drop pairs touching its epitope_footprint")
    ap.add_argument("--min-plddt", type=float, default=None, help="flag pairs with an anchor below this pLDDT")
    args = ap.parse_args(argv)

    result = suggest_disulphides(
        args.input_filepath, args.output_folder,
        group_a=args.group_a, group_b=args.group_b,
        cb_cb=(args.cb_min, args.cb_max), ca_ca_max=args.ca_max,
        min_seq_sep=args.min_seq_sep, exclude_footprint=args.exclude_footprint,
        min_plddt=args.min_plddt,
    )
    print(
        f"{Path(args.input_filepath).name}: {result['n_candidates']} disulphide candidate(s)"
        + (f", {result['n_excluded_by_footprint']} dropped by footprint" if result["exclude_footprint"] else "")
    )
    for c in result["candidates"][:10]:
        a, b = c["res_a"], c["res_b"]
        flag = "  [low pLDDT]" if c["low_plddt"] else ""
        print(
            f"  {a['resname']}{a['resseq']}/{a['chain']} (pLDDT {a['plddt']}) ↔ "
            f"{b['resname']}{b['resseq']}/{b['chain']} (pLDDT {b['plddt']})  "
            f"Sγ-Sγ={c['ss']} Å  χ3={c['chi3']}°  Cβ-Cβ={c['cb_cb']} Å{flag}"
        )


if __name__ == "__main__":
    main()
