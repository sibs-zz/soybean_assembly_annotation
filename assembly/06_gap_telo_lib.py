#!/usr/bin/env python3
"""Auto telomere extension + local gap fill helpers for 06_gap_telo.sh."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path


TELO_PLUS = "AAACCCT"
TELO_MINUS = "CCCTAAA"
TELO_TRI_PLUS = TELO_PLUS * 3
TELO_TRI_MINUS = TELO_MINUS * 3


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name = None
    chunks: list[str] = []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                name = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks).upper()
    return seqs


def write_fasta(path: Path, seqs: dict[str, str], width: int = 80) -> None:
    with path.open("w") as out:
        for name, seq in seqs.items():
            out.write(f">{name}\n")
            for i in range(0, len(seq), width):
                out.write(seq[i : i + width] + "\n")


def revcomp(seq: str) -> str:
    tbl = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(tbl)[::-1]


def parse_telo_info(path: Path) -> dict[str, dict]:
    """Parse quarTeT *.telo.info -> {chr: {status,leftnum,rightnum,...}}."""
    info: dict[str, dict] = {}
    with path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            if line.startswith("Chrid"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            chrid, chrlen, status, leftnum, leftdir, rightnum, rightdir = parts[:7]
            info[chrid] = {
                "chrlen": int(chrlen),
                "status": status,
                "leftnum": int(leftnum),
                "leftdir": leftdir,
                "rightnum": int(rightnum),
                "rightdir": rightdir,
                "has_left": int(leftnum) > 0,
                "has_right": int(rightnum) > 0,
            }
    return info


def top_scaffolds_by_length(seqs: dict[str, str], n: int) -> list[str]:
    return [k for k, _ in sorted(seqs.items(), key=lambda kv: -len(kv[1]))[:n]]


def count_telo(seq: str) -> tuple[int, int]:
    s = seq.upper()
    return s.count(TELO_PLUS), s.count(TELO_MINUS)


def parse_gaps_tsv(path: Path) -> list[dict]:
    gaps = []
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gaps.append(
                {
                    "scaffold": row["scaffold"],
                    "gap_index": int(row["gap_index"]),
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                    "length": int(row["length"]),
                }
            )
    return gaps


def write_audit_summary(
    out_path: Path,
    top: list[str],
    telo: dict[str, dict],
    gaps: list[dict],
    seqs: dict[str, str],
) -> None:
    gap_by = defaultdict(list)
    for g in gaps:
        gap_by[g["scaffold"]].append(g)

    both = 0
    missing_left = []
    missing_right = []
    with_gap = []

    lines = []
    lines.append(f"top_n\t{len(top)}")
    lines.append("scaffold\tlength\ttelo_status\thas_left\thas_right\tgap_n\tgap_bp")
    for name in top:
        t = telo.get(name)
        if t is None:
            status, hl, hr = "NA", False, False
        else:
            status, hl, hr = t["status"], t["has_left"], t["has_right"]
        gs = gap_by.get(name, [])
        gap_bp = sum(x["length"] for x in gs)
        if hl and hr:
            both += 1
        if not hl:
            missing_left.append(name)
        if not hr:
            missing_right.append(name)
        if gs:
            with_gap.append(name)
        lines.append(
            f"{name}\t{len(seqs[name])}\t{status}\t{int(hl)}\t{int(hr)}\t{len(gs)}\t{gap_bp}"
        )

    lines.append("")
    lines.append(f"top20_both_telomere\t{both}/{len(top)}\t{'YES' if both == len(top) else 'NO'}")
    lines.append(f"missing_left\t{','.join(missing_left) if missing_left else 'none'}")
    lines.append(f"missing_right\t{','.join(missing_right) if missing_right else 'none'}")
    lines.append(f"top20_with_gap\t{','.join(with_gap) if with_gap else 'none'}")
    lines.append(f"all_gap_count\t{len(gaps)}")
    lines.append(f"all_gap_bp\t{sum(g['length'] for g in gaps)}")

    out_path.write_text("\n".join(lines) + "\n")


def parse_paf(path: Path):
    rows = []
    if not path.exists() or path.stat().st_size == 0:
        return rows
    with path.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            rows.append(
                {
                    "qid": f[0],
                    "qlen": int(f[1]),
                    "qstart": int(f[2]),
                    "qend": int(f[3]),
                    "strand": f[4],
                    "tid": f[5],
                    "tlen": int(f[6]),
                    "tstart": int(f[7]),
                    "tend": int(f[8]),
                    "nmatch": int(f[9]),
                    "alen": int(f[10]),
                    "mapq": int(f[11]),
                }
            )
    return rows


def load_reads(path: Path, ids: set[str]) -> dict[str, str]:
    want = set(ids)
    got: dict[str, str] = {}
    if not want:
        return got
    name = None
    chunks: list[str] = []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None and name in want:
                    got[name] = "".join(chunks).upper()
                    if len(got) == len(want):
                        break
                name = line[1:].strip().split()[0]
                chunks = []
            else:
                if name in want:
                    chunks.append(line.strip())
        if name is not None and name in want and name not in got:
            got[name] = "".join(chunks).upper()
    return got


def best_telo_extension(
    paf_rows: list[dict],
    side: str,
    reads: dict[str, str],
    min_alen: int,
    min_ident: float,
    min_telo_copies: int,
    edge_bp: int,
    min_overhang: int = 50,
) -> dict | None:
    """Pick best overhang containing telomere motif for left/right end."""
    cands = []
    for r in paf_rows:
        if r["alen"] < min_alen:
            continue
        ident = r["nmatch"] / r["alen"] if r["alen"] else 0
        if ident < min_ident:
            continue
        seq = reads.get(r["qid"])
        if not seq:
            continue

        if side == "left":
            # alignment should touch near bait start
            if r["tstart"] > edge_bp:
                continue
            if r["strand"] == "+":
                overhang = seq[: r["qstart"]]
            else:
                # - strand: genome left corresponds to query after qend
                overhang = revcomp(seq[r["qend"] :])
        else:
            # right: alignment near bait end
            if r["tlen"] - r["tend"] > edge_bp:
                continue
            if r["strand"] == "+":
                overhang = seq[r["qend"] :]
            else:
                overhang = revcomp(seq[: r["qstart"]])

        if len(overhang) < min_overhang:
            continue
        plus_n, minus_n = count_telo(overhang)
        copies = max(plus_n, minus_n)
        if copies < min_telo_copies:
            continue
        # prefer more telomere copies, then longer overhang, then better identity
        cands.append(
            {
                "read": r["qid"],
                "strand": r["strand"],
                "ident": round(ident, 4),
                "alen": r["alen"],
                "overhang_len": len(overhang),
                "telo_plus": plus_n,
                "telo_minus": minus_n,
                "telo_copies": copies,
                "seq": overhang,
                "score": (copies, len(overhang), ident),
            }
        )
    if not cands:
        return None
    cands.sort(key=lambda x: x["score"], reverse=True)
    return cands[0]


def fill_telomeres(
    seqs: dict[str, str],
    telo_info: dict[str, dict],
    top: list[str],
    paf_path: Path,
    reads_fa: Path,
    min_alen: int,
    min_ident: float,
    min_telo_copies: int,
    missing_file: Path | None = None,
    edge_bp: int = 2000,
    min_overhang: int = 50,
) -> tuple[dict[str, str], list[dict]]:
    paf = parse_paf(paf_path)
    by_target: dict[str, list] = defaultdict(list)
    for r in paf:
        by_target[r["tid"]].append(r)

    need_ids = set()
    jobs: list[tuple[str, str, str]] = []

    if missing_file and missing_file.is_file() and missing_file.stat().st_size > 0:
        for line in missing_file.read_text().splitlines():
            if not line.strip():
                continue
            name, side = line.split("\t")[:2]
            jobs.append((name, side, f"{name}__{side}"))
    else:
        for name in top:
            t = telo_info.get(name, {})
            if not t.get("has_left", False):
                jobs.append((name, "left", f"{name}__left"))
            if not t.get("has_right", False):
                jobs.append((name, "right", f"{name}__right"))

    for _, _, tid in jobs:
        for r in by_target.get(tid, []):
            need_ids.add(r["qid"])
    reads = load_reads(reads_fa, need_ids)

    records = []
    out = dict(seqs)
    for name, side, tid in jobs:
        best = best_telo_extension(
            by_target.get(tid, []),
            side,
            reads,
            min_alen=min_alen,
            min_ident=min_ident,
            min_telo_copies=min_telo_copies,
            edge_bp=edge_bp,
            min_overhang=min_overhang,
        )
        if best is None and min_telo_copies > 3:
            # 放宽一轮：降低端粒重复要求
            best = best_telo_extension(
                by_target.get(tid, []),
                side,
                reads,
                min_alen=max(300, min_alen // 2),
                min_ident=max(0.55, min_ident - 0.1),
                min_telo_copies=3,
                edge_bp=edge_bp * 2,
                min_overhang=max(30, min_overhang // 2),
            )
        if best is None:
            records.append(
                {
                    "scaffold": name,
                    "side": side,
                    "status": "FAIL",
                    "read": "",
                    "added_bp": 0,
                    "telo_copies": 0,
                    "ident": "",
                    "note": "no_confident_overhang",
                }
            )
            continue
        if side == "left":
            out[name] = best["seq"] + out[name]
        else:
            out[name] = out[name] + best["seq"]
        records.append(
            {
                "scaffold": name,
                "side": side,
                "status": "OK",
                "read": best["read"],
                "added_bp": best["overhang_len"],
                "telo_copies": best["telo_copies"],
                "ident": best["ident"],
                "note": f"plus={best['telo_plus']};minus={best['telo_minus']}",
            }
        )
    return out, records


def best_gap_bridge(
    left_hits: list[dict],
    right_hits: list[dict],
    reads: dict[str, str],
    min_alen: int,
    min_ident: float,
    max_fill: int,
    gap_len: int = 0,
    allow_overlap: int = 0,
) -> dict | None:
    left_by = defaultdict(list)
    right_by = defaultdict(list)
    for r in left_hits:
        if r["alen"] >= min_alen and (r["nmatch"] / r["alen"]) >= min_ident:
            left_by[r["qid"]].append(r)
    for r in right_hits:
        if r["alen"] >= min_alen and (r["nmatch"] / r["alen"]) >= min_ident:
            right_by[r["qid"]].append(r)

    common = set(left_by) & set(right_by)
    cands = []
    for qid in common:
        seq = reads.get(qid)
        if not seq:
            continue
        for L in left_by[qid]:
            for R in right_by[qid]:
                # same strand relative orientation preferred
                if L["strand"] != R["strand"]:
                    continue
                if L["strand"] == "+":
                    # on read: left flank then fill then right flank
                    overlap = L["qend"] - R["qstart"]
                    if overlap > allow_overlap:
                        continue
                    fill = seq[L["qend"] : R["qstart"]]
                else:
                    # - strand: right comes first on read coordinates in genomic order
                    overlap = R["qend"] - L["qstart"]
                    if overlap > allow_overlap:
                        continue
                    fill = revcomp(seq[R["qend"] : L["qstart"]])
                if not fill or len(fill) > max_fill:
                    continue
                if gap_len > 200 and not (
                    max(1, gap_len // 3) <= len(fill) <= max(max_fill, gap_len * 3)
                ):
                    continue
                if "N" in fill:
                    continue
                ident = min(L["nmatch"] / L["alen"], R["nmatch"] / R["alen"])
                cands.append(
                    {
                        "read": qid,
                        "fill_len": len(fill),
                        "ident": round(ident, 4),
                        "seq": fill,
                        "score": (ident, -abs(len(fill)), len(fill)),
                    }
                )
    if not cands:
        return None
    cands.sort(key=lambda x: x["score"], reverse=True)
    return cands[0]


def best_gap_unified(
    paf_rows: list[dict],
    reads: dict[str, str],
    gap_start_in_bait: int,
    gap_end_in_bait: int,
    min_alen: int,
    min_ident: float,
    max_fill: int,
    slack: int = 200,
) -> dict | None:
    """Find read segment spanning gap region on unified bait (1-based bait coords)."""
    gs = gap_start_in_bait - 1
    ge = gap_end_in_bait
    cands = []
    for r in paf_rows:
        if r["alen"] < min_alen:
            continue
        ident = r["nmatch"] / r["alen"] if r["alen"] else 0
        if ident < min_ident:
            continue
        if r["tstart"] > gs + slack or r["tend"] < ge - slack:
            continue
        seq = reads.get(r["qid"])
        if not seq:
            continue
        if r["strand"] == "+":
            q1 = r["qstart"] + max(0, gs - r["tstart"])
            q2 = r["qstart"] + (ge - r["tstart"])
            fill = seq[q1:q2]
        else:
            q1 = r["qend"] - (ge - r["tstart"])
            q2 = r["qend"] - max(0, gs - r["tstart"])
            fill = revcomp(seq[q1:q2])
        if not fill or len(fill) > max_fill or "N" in fill:
            continue
        cands.append(
            {
                "read": r["qid"],
                "fill_len": len(fill),
                "ident": round(ident, 4),
                "seq": fill,
                "score": (ident, -abs(len(fill) - (ge - gs)), len(fill)),
            }
        )
    if not cands:
        return None
    cands.sort(key=lambda x: x["score"], reverse=True)
    return cands[0]


def fill_gaps(
    seqs: dict[str, str],
    gaps: list[dict],
    paf_path: Path,
    reads_fa: Path,
    min_alen: int,
    min_ident: float,
    max_fill: int,
) -> tuple[dict[str, str], list[dict]]:
    paf = parse_paf(paf_path)
    by_tid: dict[str, list] = defaultdict(list)
    for r in paf:
        by_tid[r["tid"]].append(r)

    need = set()
    for g in gaps:
        tid_l = f"{g['scaffold']}__gap{g['gap_index']}_L"
        tid_r = f"{g['scaffold']}__gap{g['gap_index']}_R"
        for r in by_tid.get(tid_l, []) + by_tid.get(tid_r, []):
            need.add(r["qid"])
    reads = load_reads(reads_fa, need)

    # apply from end of each scaffold so coordinates remain valid within scaffold edits
    by_scaff: dict[str, list] = defaultdict(list)
    for g in gaps:
        by_scaff[g["scaffold"]].append(g)

    out = dict(seqs)
    records = []
    for scaff, glist in by_scaff.items():
        seq = out[scaff]
        # sort descending by start
        for g in sorted(glist, key=lambda x: -x["start"]):
            tid_l = f"{scaff}__gap{g['gap_index']}_L"
            tid_r = f"{scaff}__gap{g['gap_index']}_R"
            best = best_gap_bridge(
                by_tid.get(tid_l, []),
                by_tid.get(tid_r, []),
                reads,
                min_alen=min_alen,
                min_ident=min_ident,
                max_fill=max_fill,
                gap_len=g["length"],
                allow_overlap=max(500, g["length"] * 2),
            )
            if best is None:
                best = best_gap_bridge(
                    by_tid.get(tid_l, []),
                    by_tid.get(tid_r, []),
                    reads,
                    min_alen=max(200, min_alen // 2),
                    min_ident=max(0.50, min_ident - 0.15),
                    max_fill=max(max_fill, g["length"] * 5),
                    gap_len=g["length"],
                    allow_overlap=max(2000, g["length"] * 5),
                )
            if best is None:
                records.append(
                    {
                        "scaffold": scaff,
                        "gap_index": g["gap_index"],
                        "old_start": g["start"],
                        "old_end": g["end"],
                        "old_length": g["length"],
                        "status": "FAIL",
                        "read": "",
                        "fill_bp": 0,
                        "ident": "",
                        "note": "no_spanning_read",
                    }
                )
                continue
            # replace N-run [start,end] 1-based inclusive
            left = seq[: g["start"] - 1]
            right = seq[g["end"] :]
            seq = left + best["seq"] + right
            records.append(
                {
                    "scaffold": scaff,
                    "gap_index": g["gap_index"],
                    "old_start": g["start"],
                    "old_end": g["end"],
                    "old_length": g["length"],
                    "status": "OK",
                    "read": best["read"],
                    "fill_bp": best["fill_len"],
                    "ident": best["ident"],
                    "note": "spanning_ont",
                }
            )
        out[scaff] = seq
    return out, records


def fill_gaps_unified(
    seqs: dict[str, str],
    gaps: list[dict],
    paf_path: Path,
    reads_fa: Path,
    flank: int,
    min_alen: int,
    min_ident: float,
    max_fill: int,
) -> tuple[dict[str, str], list[dict]]:
    paf = parse_paf(paf_path)
    by_tid: dict[str, list] = defaultdict(list)
    for r in paf:
        by_tid[r["tid"]].append(r)

    bait_meta: dict[str, dict] = {}
    for g in gaps:
        seq = seqs[g["scaffold"]]
        l1 = max(1, g["start"] - flank)
        r2 = min(len(seq), g["end"] + flank)
        tid = f"{g['scaffold']}__gap{g['gap_index']}"
        bait_meta[tid] = {
            "gap_start_in_bait": g["start"] - l1 + 1,
            "gap_end_in_bait": g["end"] - l1 + 1,
        }

    need = set()
    for tid in bait_meta:
        for r in by_tid.get(tid, []):
            need.add(r["qid"])
    reads = load_reads(reads_fa, need)

    by_scaff: dict[str, list] = defaultdict(list)
    for g in gaps:
        by_scaff[g["scaffold"]].append(g)

    out = dict(seqs)
    records = []
    for scaff, glist in by_scaff.items():
        seq = out[scaff]
        for g in sorted(glist, key=lambda x: -x["start"]):
            tid = f"{scaff}__gap{g['gap_index']}"
            meta = bait_meta[tid]
            best = best_gap_unified(
                by_tid.get(tid, []),
                reads,
                meta["gap_start_in_bait"],
                meta["gap_end_in_bait"],
                min_alen=min_alen,
                min_ident=min_ident,
                max_fill=max(max_fill, g["length"] * 5),
            )
            if best is None:
                best = best_gap_unified(
                    by_tid.get(tid, []),
                    reads,
                    meta["gap_start_in_bait"],
                    meta["gap_end_in_bait"],
                    min_alen=max(200, min_alen // 2),
                    min_ident=max(0.50, min_ident - 0.15),
                    max_fill=max(max_fill, g["length"] * 10),
                    slack=500,
                )
            if best is None:
                records.append(
                    {
                        "scaffold": scaff,
                        "gap_index": g["gap_index"],
                        "old_start": g["start"],
                        "old_end": g["end"],
                        "old_length": g["length"],
                        "status": "FAIL",
                        "read": "",
                        "fill_bp": 0,
                        "ident": "",
                        "note": "no_unified_span",
                    }
                )
                continue
            left = seq[: g["start"] - 1]
            right = seq[g["end"] :]
            seq = left + best["seq"] + right
            records.append(
                {
                    "scaffold": scaff,
                    "gap_index": g["gap_index"],
                    "old_start": g["start"],
                    "old_end": g["end"],
                    "old_length": g["length"],
                    "status": "OK",
                    "read": best["read"],
                    "fill_bp": best["fill_len"],
                    "ident": best["ident"],
                    "note": "unified_ont",
                }
            )
        out[scaff] = seq
    return out, records


def cmd_audit(args: argparse.Namespace) -> None:
    seqs = read_fasta(Path(args.fasta))
    top = top_scaffolds_by_length(seqs, args.top_n)
    Path(args.top_list).write_text("\n".join(top) + "\n")
    telo = parse_telo_info(Path(args.telo_info)) if args.telo_info else {}
    gaps = parse_gaps_tsv(Path(args.gaps_tsv)) if args.gaps_tsv else []
    # restrict gap list note: keep all gaps in genome for later fill, but summarize top20
    write_audit_summary(Path(args.summary), top, telo, gaps, seqs)
    # missing ends list for bait building
    miss = []
    for name in top:
        t = telo.get(name, {})
        if not t.get("has_left", False):
            miss.append(f"{name}\tleft")
        if not t.get("has_right", False):
            miss.append(f"{name}\tright")
    Path(args.missing_telo).write_text("\n".join(miss) + ("\n" if miss else ""))
    print(f"[audit] top{args.top_n} both_telo="
          f"{sum(1 for n in top if telo.get(n, {}).get('has_left') and telo.get(n, {}).get('has_right'))}"
          f"/{len(top)}; gaps={len(gaps)}; missing_ends={len(miss)}")


def cmd_make_telo_baits(args: argparse.Namespace) -> None:
    seqs = read_fasta(Path(args.fasta))
    bait_len = args.bait_len
    out_fa = Path(args.out_fa)
    with out_fa.open("w") as out:
        for line in Path(args.missing_telo).read_text().splitlines():
            if not line.strip():
                continue
            name, side = line.split("\t")
            seq = seqs[name]
            if side == "left":
                bait = seq[:bait_len]
                out.write(f">{name}__left\n")
            else:
                bait = seq[-bait_len:]
                out.write(f">{name}__right\n")
            for i in range(0, len(bait), 80):
                out.write(bait[i : i + 80] + "\n")


def cmd_make_gap_baits(args: argparse.Namespace) -> None:
    seqs = read_fasta(Path(args.fasta))
    gaps = parse_gaps_tsv(Path(args.gaps_tsv))
    flank = args.flank
    out_fa = Path(args.out_fa)
    unified = getattr(args, "unified", False)
    with out_fa.open("w") as out:
        for g in gaps:
            seq = seqs[g["scaffold"]]
            if unified:
                l1 = max(1, g["start"] - flank)
                r2 = min(len(seq), g["end"] + flank)
                bait = seq[l1 - 1 : r2]
                tid = f"{g['scaffold']}__gap{g['gap_index']}"
                if not bait:
                    continue
                out.write(f">{tid}\n")
                for i in range(0, len(bait), 80):
                    out.write(bait[i : i + 80] + "\n")
                continue
            l1 = max(1, g["start"] - flank)
            l2 = g["start"] - 1
            r1 = g["end"] + 1
            r2 = min(len(seq), g["end"] + flank)
            left = seq[l1 - 1 : l2]
            right = seq[r1 - 1 : r2]
            tid_l = f"{g['scaffold']}__gap{g['gap_index']}_L"
            tid_r = f"{g['scaffold']}__gap{g['gap_index']}_R"
            for tid, bait in ((tid_l, left), (tid_r, right)):
                if not bait:
                    continue
                out.write(f">{tid}\n")
                for i in range(0, len(bait), 80):
                    out.write(bait[i : i + 80] + "\n")


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})


def cmd_apply_telo(args: argparse.Namespace) -> None:
    seqs = read_fasta(Path(args.fasta))
    top = Path(args.top_list).read_text().splitlines()
    telo = parse_telo_info(Path(args.telo_info))
    out_seqs, recs = fill_telomeres(
        seqs,
        telo,
        top,
        Path(args.paf),
        Path(args.reads_fa),
        min_alen=args.min_alen,
        min_ident=args.min_ident,
        min_telo_copies=args.min_telo_copies,
        missing_file=Path(args.missing_telo) if args.missing_telo else None,
        edge_bp=args.edge_bp,
        min_overhang=args.min_overhang,
    )
    write_fasta(Path(args.out_fa), out_seqs)
    write_tsv(
        Path(args.record),
        recs,
        ["scaffold", "side", "status", "read", "added_bp", "telo_copies", "ident", "note"],
    )
    ok = sum(1 for r in recs if r["status"] == "OK")
    print(f"[telo_fill] ok={ok}/{len(recs)} -> {args.out_fa}")


def cmd_apply_gap(args: argparse.Namespace) -> None:
    seqs = read_fasta(Path(args.fasta))
    gaps = parse_gaps_tsv(Path(args.gaps_tsv))
    if getattr(args, "unified", False):
        out_seqs, recs = fill_gaps_unified(
            seqs,
            gaps,
            Path(args.paf),
            Path(args.reads_fa),
            flank=args.flank,
            min_alen=args.min_alen,
            min_ident=args.min_ident,
            max_fill=args.max_fill,
        )
    else:
        out_seqs, recs = fill_gaps(
            seqs,
            gaps,
            Path(args.paf),
            Path(args.reads_fa),
            min_alen=args.min_alen,
            min_ident=args.min_ident,
            max_fill=args.max_fill,
        )
    write_fasta(Path(args.out_fa), out_seqs)
    write_tsv(
        Path(args.record),
        recs,
        [
            "scaffold",
            "gap_index",
            "old_start",
            "old_end",
            "old_length",
            "status",
            "read",
            "fill_bp",
            "ident",
            "note",
        ],
    )
    ok = sum(1 for r in recs if r["status"] == "OK")
    print(f"[gap_fill] ok={ok}/{len(recs)} -> {args.out_fa}")


def main() -> None:
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)

    p = sp.add_parser("audit")
    p.add_argument("--fasta", required=True)
    p.add_argument("--telo-info", required=True)
    p.add_argument("--gaps-tsv", required=True)
    p.add_argument("--top-n", type=int, default=20)
    p.add_argument("--top-list", required=True)
    p.add_argument("--missing-telo", required=True)
    p.add_argument("--summary", required=True)
    p.set_defaults(func=cmd_audit)

    p = sp.add_parser("make-telo-baits")
    p.add_argument("--fasta", required=True)
    p.add_argument("--missing-telo", required=True)
    p.add_argument("--bait-len", type=int, default=50000)
    p.add_argument("--out-fa", required=True)
    p.set_defaults(func=cmd_make_telo_baits)

    p = sp.add_parser("make-gap-baits")
    p.add_argument("--fasta", required=True)
    p.add_argument("--gaps-tsv", required=True)
    p.add_argument("--flank", type=int, default=20000)
    p.add_argument("--out-fa", required=True)
    p.add_argument("--unified", action="store_true", help="one bait per gap including N-run")
    p.set_defaults(func=cmd_make_gap_baits)

    p = sp.add_parser("apply-telo")
    p.add_argument("--fasta", required=True)
    p.add_argument("--telo-info", required=True)
    p.add_argument("--top-list", required=True)
    p.add_argument("--missing-telo", default="")
    p.add_argument("--paf", required=True)
    p.add_argument("--reads-fa", required=True)
    p.add_argument("--out-fa", required=True)
    p.add_argument("--record", required=True)
    p.add_argument("--min-alen", type=int, default=500)
    p.add_argument("--min-ident", type=float, default=0.65)
    p.add_argument("--min-telo-copies", type=int, default=8)
    p.add_argument("--edge-bp", type=int, default=5000)
    p.add_argument("--min-overhang", type=int, default=50)
    p.set_defaults(func=cmd_apply_telo)

    p = sp.add_parser("apply-gap")
    p.add_argument("--fasta", required=True)
    p.add_argument("--gaps-tsv", required=True)
    p.add_argument("--paf", required=True)
    p.add_argument("--reads-fa", required=True)
    p.add_argument("--out-fa", required=True)
    p.add_argument("--record", required=True)
    p.add_argument("--flank", type=int, default=50000)
    p.add_argument("--unified", action="store_true")
    p.add_argument("--min-alen", type=int, default=500)
    p.add_argument("--min-ident", type=float, default=0.60)
    p.add_argument("--max-fill", type=int, default=5000)
    p.set_defaults(func=cmd_apply_gap)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
