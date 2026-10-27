#!/usr/bin/env python3
"""Helpers for 06_gap_telo_final_fix.sh (does NOT modify 06_gap_telo_lib.py).

Motif-based telomere calling + diagnose/trim/fill routing with ONT evidence
for buried tips (trim / keep / extend / unresolved).

Commands:
  call-motif, diagnose-orient, make-buried-baits, decide-buried-ont,
  trim-buried, accept-trial-trim, apply-revcomp, check-perfect
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from importlib import import_module

_lib = import_module("06_gap_telo_lib")  # noqa: E402

read_fasta = _lib.read_fasta
write_fasta = _lib.write_fasta
revcomp = _lib.revcomp
parse_telo_info = _lib.parse_telo_info
parse_gaps_tsv = _lib.parse_gaps_tsv
top_scaffolds_by_length = _lib.top_scaffolds_by_length
TELO_PLUS = _lib.TELO_PLUS
TELO_MINUS = _lib.TELO_MINUS

TELO_PLANT_MINUS = "AGGGTTT"


def motif_count(seq: str) -> int:
    s = seq.upper()
    return s.count(TELO_PLUS) + s.count(TELO_MINUS) + s.count(TELO_PLANT_MINUS)


def tip_copies(seq: str, tip_bp: int = 20000) -> tuple[int, int]:
    """Return (left_total, right_total) monomer counts."""
    s = seq.upper()
    left = s[:tip_bp]
    right = s[-tip_bp:] if len(s) >= tip_bp else s
    return motif_count(left), motif_count(right)


def resolve_top(seqs: dict[str, str], top_list: str, top_n: int) -> list[str]:
    if top_list and Path(top_list).is_file():
        names = [ln.strip() for ln in Path(top_list).read_text().splitlines() if ln.strip()]
        if names:
            return names
    return top_scaffolds_by_length(seqs, top_n)


def write_synthetic_telo_info(
    path: Path,
    seqs: dict[str, str],
    top: list[str],
    left_c: dict[str, int],
    right_c: dict[str, int],
    min_copies: int,
) -> None:
    """quarTeT-compatible .telo.info so apply-telo can still parse it."""
    lines = [
        "# synthetic motif-based telo.info (not from quarTeT)",
        "Chrid\tChrlen\tStatus\tLeftNum\tLeftDir\tRightNum\tRightDir",
    ]
    for name in top:
        if name not in seqs:
            continue
        hl = left_c[name] >= min_copies
        hr = right_c[name] >= min_copies
        if hl and hr:
            status = "both"
        elif hl:
            status = "left"
        elif hr:
            status = "right"
        else:
            status = "no"
        lines.append(
            f"{name}\t{len(seqs[name])}\t{status}\t"
            f"{left_c[name] if hl else 0}\t+\t{right_c[name] if hr else 0}\t-"
        )
    path.write_text("\n".join(lines) + "\n")


def cmd_call_motif(args: argparse.Namespace) -> None:
    seqs = read_fasta(Path(args.fasta))
    top = resolve_top(seqs, args.top_list, args.top_n)
    Path(args.out_top_list).write_text("\n".join(top) + "\n")

    left_c: dict[str, int] = {}
    right_c: dict[str, int] = {}
    rows = []
    miss = []
    both = 0
    for name in top:
        seq = seqs[name]
        lc, rc = tip_copies(seq, args.tip_bp)
        left_c[name], right_c[name] = lc, rc
        hl, hr = lc >= args.min_copies, rc >= args.min_copies
        if hl and hr:
            status, both = "both", both + 1
        elif hl:
            status = "left"
        elif hr:
            status = "right"
        else:
            status = "no"
        if not hl:
            miss.append(f"{name}\tleft")
        if not hr:
            miss.append(f"{name}\tright")
        rows.append(
            {
                "scaffold": name,
                "length": len(seq),
                "left_tip_copies": lc,
                "right_tip_copies": rc,
                "has_left": int(hl),
                "has_right": int(hr),
                "status": status,
            }
        )

    fields = [
        "scaffold",
        "length",
        "left_tip_copies",
        "right_tip_copies",
        "has_left",
        "has_right",
        "status",
    ]
    out = Path(args.out_tsv)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    Path(args.out_missing).write_text("\n".join(miss) + ("\n" if miss else ""))
    write_synthetic_telo_info(
        Path(args.out_telo_info), seqs, top, left_c, right_c, args.min_copies
    )

    gaps = parse_gaps_tsv(Path(args.gaps_tsv)) if args.gaps_tsv and Path(args.gaps_tsv).is_file() else []
    summary = [
        f"mode\tmotif",
        f"tip_bp\t{args.tip_bp}",
        f"min_copies\t{args.min_copies}",
        f"top_n\t{len(top)}",
        f"top_both_telomere\t{both}/{len(top)}\t{'YES' if both == len(top) else 'NO'}",
        f"missing_ends\t{len(miss)}",
        f"all_gap_count\t{len(gaps)}",
        f"missing_list\t{','.join(f'{a}:{b}' for a, b in (x.split(chr(9)) for x in miss)) if miss else 'none'}",
    ]
    Path(args.out_summary).write_text("\n".join(summary) + "\n")
    print(
        f"[call-motif] both={both}/{len(top)} missing_ends={len(miss)} "
        f"gaps={len(gaps)} -> {out}"
    )


def find_trim_for_buried(
    seq: str,
    side: str,
    tip_bp: int,
    tip_min_copies: int,
    internal_min_copies: int,
    max_trim_bp: int,
    step: int | None = None,
) -> dict:
    """
    Locate near-tip telomere island and propose trim of distal junk.
    Returns dict with keys: ok, action, trim_bp, keep_start, keep_end, ...
    """
    n = len(seq)
    s = seq.upper()
    step = step or max(tip_bp // 4, 5000)
    tip_c = tip_copies(seq, tip_bp)
    tip_now = tip_c[0] if side == "left" else tip_c[1]
    if tip_now >= tip_min_copies:
        return {
            "ok": False,
            "action": "skip",
            "trim_bp": 0,
            "note": f"tip_already_ok copies={tip_now}",
        }

    cands: list[tuple[int, int, int]] = []  # copies, start, end
    if side == "right":
        zone0 = max(0, n - max_trim_bp - tip_bp)
        for i in range(zone0, max(zone0, n - tip_bp) + 1, step):
            if i + tip_bp > n:
                break
            c = motif_count(s[i : i + tip_bp])
            if c >= internal_min_copies:
                cands.append((c, i, i + tip_bp))
        if not cands:
            return {
                "ok": False,
                "action": "report",
                "trim_bp": 0,
                "note": f"no_island_within_max_trim={max_trim_bp}",
            }
        # distal-most island edge (largest end), break ties by copies
        best = max(cands, key=lambda x: (x[2], x[0]))
        keep_start, keep_end = 0, best[2]
        trim_bp = n - keep_end
    else:
        zone1 = min(n, max_trim_bp + tip_bp)
        for i in range(0, max(0, zone1 - tip_bp) + 1, step):
            c = motif_count(s[i : i + tip_bp])
            if c >= internal_min_copies:
                cands.append((c, i, i + tip_bp))
        if not cands:
            return {
                "ok": False,
                "action": "report",
                "trim_bp": 0,
                "note": f"no_island_within_max_trim={max_trim_bp}",
            }
        best = min(cands, key=lambda x: (x[1], -x[0]))  # proximal-most start
        keep_start, keep_end = best[1], n
        trim_bp = keep_start

    if trim_bp <= 0:
        return {
            "ok": False,
            "action": "skip",
            "trim_bp": 0,
            "note": "island_already_at_tip",
        }
    if trim_bp > max_trim_bp:
        return {
            "ok": False,
            "action": "report",
            "trim_bp": trim_bp,
            "cluster_start": best[1],
            "cluster_end": best[2],
            "cluster_copies": best[0],
            "note": f"trim_bp={trim_bp}>{max_trim_bp} exceed_max_no_trim",
        }

    # preview tip after trim
    kept = s[keep_start:keep_end]
    new_tip = motif_count(kept[:tip_bp] if side == "left" else kept[-tip_bp:])
    return {
        "ok": True,
        "action": "trim",
        "trim_bp": trim_bp,
        "keep_start": keep_start,
        "keep_end": keep_end,
        "cluster_start": best[1],
        "cluster_end": best[2],
        "cluster_copies": best[0],
        "new_tip_copies": new_tip,
        "note": f"trim_{side}_{trim_bp}bp_island@{best[1]}-{best[2]}_c={best[0]}",
    }


def classify_missing_side(
    seq: str,
    side: str,
    tip_bp: int,
    tip_min_copies: int,
    internal_min_copies: int,
    near_tip_bp: int,
    max_trim_bp: int,
    trust_missing: bool = False,
) -> tuple[str, str, dict]:
    """
    Returns (class, note, trim_plan_or_empty).
    class: tip_has_telo | true_missing | buried_near_tip | deep_internal

    trust_missing=True：缺失列表来自 quarTeT 等权威源时，不因 tip motif
    达标而 skip（仍可继续判 buried / true_missing）。
    """
    lp, rp = tip_copies(seq, tip_bp)
    tip_c = lp if side == "left" else rp
    if tip_c >= tip_min_copies and not trust_missing:
        return "tip_has_telo", f"tip_copies={tip_c}", {"action": "skip"}

    n = len(seq)
    win = tip_bp
    step = max(tip_bp // 2, 10000)
    near_best = deep_best = 0
    near_pos = deep_pos = -1
    for i in range(0, max(0, n - win) + 1, step):
        c = motif_count(seq[i : i + win].upper())
        if c <= 0:
            continue
        mid = i + win // 2
        dist = mid if side == "left" else n - mid
        if dist <= near_tip_bp:
            if c > near_best:
                near_best, near_pos = c, i
        else:
            if c > deep_best:
                deep_best, deep_pos = c, i

    if tip_c < tip_min_copies and near_best >= internal_min_copies:
        plan = find_trim_for_buried(
            seq, side, tip_bp, tip_min_copies, internal_min_copies, max_trim_bp
        )
        note = f"tip_copies={tip_c}; near_tip_best={near_best}@{near_pos}"
        if plan.get("action") == "trim":
            return "buried_near_tip", note + f"; will_trim={plan['trim_bp']}", plan
        if plan.get("action") == "report":
            return "buried_near_tip", note + f"; {plan.get('note', '')}", plan
        return "buried_near_tip", note, plan

    if deep_best >= internal_min_copies:
        return (
            "deep_internal",
            f"tip_copies={tip_c}; deep_best={deep_best}@{deep_pos}",
            {"action": "report"},
        )
    return (
        "true_missing",
        f"tip_copies={tip_c}; near={near_best}; deep={deep_best}",
        {"action": "fill"},
    )


def cmd_diagnose_orient(args: argparse.Namespace) -> None:
    seqs = read_fasta(Path(args.fasta))
    top = resolve_top(seqs, args.top_list, args.top_n)
    quartet = (
        parse_telo_info(Path(args.telo_info))
        if args.telo_info and Path(args.telo_info).is_file()
        else {}
    )

    # Prefer motif missing list if provided; else motif-call on the fly
    missing_jobs: list[tuple[str, str]] = []
    if args.missing_telo and Path(args.missing_telo).is_file():
        for line in Path(args.missing_telo).read_text().splitlines():
            if not line.strip():
                continue
            name, side = line.split("\t")[:2]
            missing_jobs.append((name, side))
    else:
        for name in top:
            lc, rc = tip_copies(seqs[name], args.tip_bp)
            if lc < args.tip_min_copies:
                missing_jobs.append((name, "left"))
            if rc < args.tip_min_copies:
                missing_jobs.append((name, "right"))

    fields = [
        "scaffold",
        "length",
        "missing_side",
        "quartet_status",
        "left_tip_copies",
        "right_tip_copies",
        "class",
        "action",
        "trim_bp",
        "suggest_revcomp",
        "note",
    ]
    rows = []
    revcomp_targets: list[str] = []
    trim_rows = []
    fill_miss: list[str] = []

    for name, side in missing_jobs:
        if name not in seqs:
            continue
        seq = seqs[name]
        lc, rc = tip_copies(seq, args.tip_bp)
        tip_c = lc if side == "left" else rc
        q = quartet.get(name, {})
        cls, note, plan = classify_missing_side(
            seq,
            side,
            tip_bp=args.tip_bp,
            tip_min_copies=args.tip_min_copies,
            internal_min_copies=args.internal_min_copies,
            near_tip_bp=args.near_tip_bp,
            max_trim_bp=args.max_trim_bp,
            trust_missing=bool(args.trust_missing),
        )
        # quarTeT 缺端但 tip motif 已高：仍尝试 fill
        if (
            args.trust_missing
            and tip_c >= args.tip_min_copies
            and cls != "buried_near_tip"
            and cls != "deep_internal"
        ):
            cls = "true_missing"
            note = f"quartet_missing_motif_tip={tip_c}; {note}"
            plan = {"action": "fill"}
        action = plan.get("action", "report")
        if cls == "tip_has_telo":
            action = "skip"
        elif cls == "true_missing":
            action = "fill"
            fill_miss.append(f"{name}\t{side}")
        elif cls == "buried_near_tip":
            if plan.get("action") == "trim" and plan.get("ok"):
                action = "buried_check"
                trim_rows.append(
                    {
                        "scaffold": name,
                        "side": side,
                        "trim_bp": plan["trim_bp"],
                        "keep_start": plan["keep_start"],
                        "keep_end": plan["keep_end"],
                        "cluster_start": plan.get("cluster_start", ""),
                        "cluster_end": plan.get("cluster_end", ""),
                        "cluster_copies": plan.get("cluster_copies", ""),
                        "new_tip_copies": plan.get("new_tip_copies", ""),
                        "note": plan.get("note", ""),
                    }
                )
            elif plan.get("action") == "report":
                action = "unresolved"
            else:
                action = plan.get("action", "report")
            # do NOT auto-fill buried tails
        elif cls == "deep_internal":
            action = "unresolved"

        suggest = 0
        if (
            args.flip_onesided
            and side == "left"
            and cls == "true_missing"
            and lc < args.tip_min_copies
            and rc >= args.tip_min_copies
        ):
            suggest = 1
            revcomp_targets.append(name)
            note = note + "; flip_onesided_right_to_left"

        rows.append(
            {
                "scaffold": name,
                "length": len(seq),
                "missing_side": side,
                "quartet_status": q.get("status", "NA"),
                "left_tip_copies": lc,
                "right_tip_copies": rc,
                "class": cls,
                "action": action,
                "trim_bp": plan.get("trim_bp", 0)
                if action in ("buried_check", "trim", "unresolved", "report")
                else 0,
                "suggest_revcomp": suggest,
                "note": note,
            }
        )

    out = Path(args.out_tsv)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    trim_path = Path(args.trim_plan)
    trim_fields = [
        "scaffold",
        "side",
        "trim_bp",
        "keep_start",
        "keep_end",
        "cluster_start",
        "cluster_end",
        "cluster_copies",
        "new_tip_copies",
        "note",
    ]
    with trim_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=trim_fields, delimiter="\t")
        w.writeheader()
        for r in trim_rows:
            w.writerow(r)

    Path(args.fill_missing).write_text("\n".join(fill_miss) + ("\n" if fill_miss else ""))
    targets = Path(args.revcomp_list)
    uniq = sorted(set(revcomp_targets))
    targets.write_text("\n".join(uniq) + ("\n" if uniq else ""))

    print(f"[diagnose-orient] rows={len(rows)} -> {out}")
    print(f"[diagnose-orient] bury_candidates={len(trim_rows)} -> {trim_path}")
    print(f"[diagnose-orient] fill_missing={len(fill_miss)} -> {args.fill_missing}")
    for r in rows:
        print(
            f"  {r['scaffold']}\t{r['missing_side']}\t{r['class']}\t"
            f"action={r['action']}\t{r['note']}"
        )


def cmd_make_buried_baits(args: argparse.Namespace) -> None:
    """Extract island-edge + distal-tail baits for ONT evidence mapping."""
    seqs = read_fasta(Path(args.fasta))
    cands = list(csv.DictReader(Path(args.candidates).open(), delimiter="\t"))
    anchor = args.anchor_bp
    bait_seqs: dict[str, str] = {}
    coord_rows = []
    for r in cands:
        name, side = r["scaffold"], r["side"]
        if name not in seqs:
            continue
        seq = seqs[name]
        n = len(seq)
        keep_start, keep_end = int(r["keep_start"]), int(r["keep_end"])
        trim_bp = int(r["trim_bp"])
        tid = f"{name}__{side}_buried"
        if side == "right":
            bait_start = max(0, keep_end - anchor)
            bait_end = n
            junction = keep_end - bait_start
        else:
            bait_start = 0
            bait_end = min(n, keep_start + anchor)
            junction = keep_start
        if bait_end - bait_start < 1000:
            continue
        bait_seqs[tid] = seq[bait_start:bait_end]
        coord_rows.append(
            {
                "scaffold": name,
                "side": side,
                "tid": tid,
                "bait_start": bait_start,
                "bait_end": bait_end,
                "junction": junction,
                "trim_bp": trim_bp,
                "keep_start": keep_start,
                "keep_end": keep_end,
                "cluster_start": r.get("cluster_start", ""),
                "cluster_end": r.get("cluster_end", ""),
                "cluster_copies": r.get("cluster_copies", ""),
                "new_tip_copies": r.get("new_tip_copies", ""),
                "note": r.get("note", ""),
            }
        )
    write_fasta(Path(args.out_fa), bait_seqs)
    fields = [
        "scaffold",
        "side",
        "tid",
        "bait_start",
        "bait_end",
        "junction",
        "trim_bp",
        "keep_start",
        "keep_end",
        "cluster_start",
        "cluster_end",
        "cluster_copies",
        "new_tip_copies",
        "note",
    ]
    with Path(args.out_coords).open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for row in coord_rows:
            w.writerow(row)
    print(f"[make-buried-baits] n={len(bait_seqs)} -> {args.out_fa}")


def _paf_ident(r: dict) -> float:
    return r["nmatch"] / r["alen"] if r["alen"] > 0 else 0.0


def _query_tip_overhang(r: dict, side: str) -> int:
    """Overhang bp beyond current assembly tip (distal end of buried bait)."""
    if side == "right":
        if r["tlen"] - r["tend"] > 200:
            return 0
        return (r["qlen"] - r["qend"]) if r["strand"] == "+" else r["qstart"]
    if r["tstart"] > 200:
        return 0
    return r["qstart"] if r["strand"] == "+" else (r["qlen"] - r["qend"])


def cmd_decide_buried_ont(args: argparse.Namespace) -> None:
    """Classify buried candidates by ONT PAF: trim / keep / extend / unresolved."""
    coords = {
        r["tid"]: r
        for r in csv.DictReader(Path(args.coords).open(), delimiter="\t")
    }
    cands = list(csv.DictReader(Path(args.candidates).open(), delimiter="\t"))
    paf = _lib.parse_paf(Path(args.paf))
    by_tid: dict[str, list] = {}
    for r in paf:
        by_tid.setdefault(r["tid"], []).append(r)

    need_ids = set()
    for tid, rows in by_tid.items():
        meta = coords.get(tid)
        if not meta:
            continue
        side = meta["side"]
        for r in rows:
            if r["alen"] < args.min_alen or _paf_ident(r) < args.min_ident:
                continue
            if _query_tip_overhang(r, side) >= args.min_overhang:
                need_ids.add(r["qid"])
    reads = _lib.load_reads(Path(args.reads_fa), need_ids) if args.reads_fa else {}

    decide_rows = []
    trim_rows = []
    extend_lines = []
    unresolved_rows = []

    for c in cands:
        name, side = c["scaffold"], c["side"]
        tid = f"{name}__{side}_buried"
        meta = coords.get(tid, {})
        junction = int(meta["junction"]) if meta else 0
        trim_bp = int(c["trim_bp"])
        rows = by_tid.get(tid, [])

        span_reads: set[str] = set()
        distal_reads: set[str] = set()
        extend_qids: set[str] = set()

        for r in rows:
            if r["alen"] < args.min_alen or _paf_ident(r) < args.min_ident:
                continue
            if r["tstart"] < junction - 500 and r["tend"] > junction + 500:
                span_reads.add(r["qid"])
            if side == "right":
                if r["tend"] - junction >= max(2000, int(trim_bp * 0.2)):
                    distal_reads.add(r["qid"])
            else:
                if junction - r["tstart"] >= max(2000, int(trim_bp * 0.2)):
                    distal_reads.add(r["qid"])

            oh = _query_tip_overhang(r, side)
            if oh >= args.min_overhang and r["qid"] in reads:
                q = reads[r["qid"]]
                if side == "right":
                    oh_seq = q[-oh:] if r["strand"] == "+" else revcomp(q[:oh])
                else:
                    oh_seq = q[:oh] if r["strand"] == "+" else revcomp(q[-oh:])
                if motif_count(oh_seq) >= args.min_telo_copies:
                    extend_qids.add(r["qid"])

        n_span, n_distal, n_ext = len(span_reads), len(distal_reads), len(extend_qids)
        distal_supported = (
            n_span >= args.min_support_reads or n_distal >= args.min_support_reads
        )

        # 有远端 ONT 支持：默认试 trim 一次（需 new_tip motif 达标且 trim≤max）
        #   脚本侧再用 quarTeT 验收；不达标则回滚 → unresolved
        # 无远端支持但 tip 有端粒 overhang：允许 extend（脚本侧只试 1 轮）
        # 完全无支持：trim（AUTO_TRIM=1）
        if distal_supported:
            new_tip = int(float(c.get("new_tip_copies") or 0))
            ks = c.get("keep_start", "")
            ke = c.get("keep_end", "")
            can_trial = (
                args.allow_distal_trial_trim
                and trim_bp > 0
                and trim_bp <= args.max_trim_bp
                and new_tip >= args.min_new_tip_copies
                and str(ks) != ""
                and str(ke) != ""
            )
            if can_trial:
                decision = "trial_trim"
                note = (
                    f"span={n_span};distal={n_distal};extend_reads={n_ext};"
                    f"distal_supported_trial_trim;new_tip={new_tip}"
                )
            else:
                decision = "unresolved"
                note = (
                    f"span={n_span};distal={n_distal};extend_reads={n_ext};"
                    f"distal_supported_no_viable_trim;"
                    f"new_tip={new_tip};trim_bp={trim_bp}"
                )
        elif n_ext >= args.min_extend_reads:
            decision = "extend"
            note = f"span={n_span};distal={n_distal};extend_reads={n_ext};tip_overhang_only"
        elif n_span == 0 and n_distal == 0:
            decision = "trim" if args.auto_trim else "unresolved"
            note = "span=0;distal=0;no_ont_support"
        else:
            decision = "unresolved"
            note = f"span={n_span};distal={n_distal};extend_reads={n_ext};ambiguous"

        decide_rows.append(
            {
                "scaffold": name,
                "side": side,
                "decision": decision,
                "span_reads": n_span,
                "distal_reads": n_distal,
                "extend_reads": n_ext,
                "trim_bp": trim_bp,
                "note": note,
            }
        )
        if decision in ("trim", "trial_trim"):
            trim_rows.append(c)
        elif decision == "extend":
            extend_lines.append(f"{name}\t{side}")
        elif decision == "unresolved":
            unresolved_rows.append(
                {
                    "scaffold": name,
                    "side": side,
                    "status": "unresolved",
                    "span_reads": n_span,
                    "distal_reads": n_distal,
                    "extend_reads": n_ext,
                    "trim_bp": trim_bp,
                    "note": note,
                }
            )

    trim_fields = [
        "scaffold",
        "side",
        "trim_bp",
        "keep_start",
        "keep_end",
        "cluster_start",
        "cluster_end",
        "cluster_copies",
        "new_tip_copies",
        "note",
    ]
    with Path(args.out_trim_plan).open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=trim_fields, delimiter="\t")
        w.writeheader()
        for r in trim_rows:
            w.writerow({k: r.get(k, "") for k in trim_fields})

    Path(args.out_extend).write_text(
        "\n".join(extend_lines) + ("\n" if extend_lines else "")
    )

    ufields = [
        "scaffold",
        "side",
        "status",
        "span_reads",
        "distal_reads",
        "extend_reads",
        "trim_bp",
        "note",
    ]
    with Path(args.out_unresolved).open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=ufields, delimiter="\t")
        w.writeheader()
        for r in unresolved_rows:
            w.writerow(r)

    dfields = [
        "scaffold",
        "side",
        "decision",
        "span_reads",
        "distal_reads",
        "extend_reads",
        "trim_bp",
        "note",
    ]
    with Path(args.out_decide).open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=dfields, delimiter="\t")
        w.writeheader()
        for r in decide_rows:
            w.writerow(r)

    print(
        f"[decide-buried-ont] trim={sum(1 for r in decide_rows if r['decision'] == 'trim')} "
        f"trial_trim={sum(1 for r in decide_rows if r['decision'] == 'trial_trim')} "
        f"extend={sum(1 for r in decide_rows if r['decision'] == 'extend')} "
        f"unresolved={len(unresolved_rows)} -> {args.out_decide}"
    )


def cmd_trim_buried(args: argparse.Namespace) -> None:
    seqs = read_fasta(Path(args.fasta))
    plan_path = Path(args.trim_plan)
    out_clipped_fa = getattr(args, "out_clipped_fa", "") or ""
    out_clipped_tsv = getattr(args, "out_clipped_tsv", "") or ""
    if not plan_path.is_file() or plan_path.stat().st_size == 0:
        write_fasta(Path(args.out_fa), seqs)
        Path(args.record).write_text(
            "scaffold\tside\tstatus\ttrim_bp\told_len\tnew_len\tnote\n"
        )
        if out_clipped_fa:
            Path(out_clipped_fa).write_text("")
        if out_clipped_tsv:
            Path(out_clipped_tsv).write_text(
                "scaffold\tside\ttrim_bp\tclipped_start\tclipped_end\tclipped_len\tnote\n"
            )
        print(f"[trim-buried] no plan, copy -> {args.out_fa}")
        return

    # one scaffold may appear twice (both sides) — apply sequentially carefully
    # Group by scaffold: apply left trim first (coords from original), then right on result
    rows = list(csv.DictReader(plan_path.open(), delimiter="\t"))
    by_scaf: dict[str, list[dict]] = {}
    for r in rows:
        by_scaf.setdefault(r["scaffold"], []).append(r)

    records = []
    clipped_recs = []
    clipped_fa_lines: list[str] = []
    n_ok = 0
    for name, plans in by_scaf.items():
        if name not in seqs:
            records.append(
                {
                    "scaffold": name,
                    "side": ",".join(p["side"] for p in plans),
                    "status": "FAIL",
                    "trim_bp": 0,
                    "old_len": "",
                    "new_len": "",
                    "note": "scaffold_missing",
                }
            )
            continue
        seq = seqs[name]
        old_len = len(seq)
        # Apply left then right using absolute coords from ORIGINAL sequence,
        # by converting to keep slice once.
        keep_l = 0
        keep_r = old_len
        total_trim = 0
        sides = []
        notes = []
        for p in plans:
            side = p["side"]
            ks, ke = int(p["keep_start"]), int(p["keep_end"])
            tb = int(p["trim_bp"])
            if tb > args.max_trim_bp:
                notes.append(f"{side}:skip_exceed")
                continue
            if side == "left":
                if ks > 0:
                    clip = seq[:ks]
                    clipped_fa_lines.append(
                        f">{name}__left_clip len={len(clip)}\n{clip}\n"
                    )
                    clipped_recs.append(
                        {
                            "scaffold": name,
                            "side": "left",
                            "trim_bp": ks,
                            "clipped_start": 0,
                            "clipped_end": ks,
                            "clipped_len": len(clip),
                            "note": p.get("note", ""),
                        }
                    )
                keep_l = max(keep_l, ks)
            else:
                if ke < old_len:
                    clip = seq[ke:]
                    clipped_fa_lines.append(
                        f">{name}__right_clip len={len(clip)}\n{clip}\n"
                    )
                    clipped_recs.append(
                        {
                            "scaffold": name,
                            "side": "right",
                            "trim_bp": old_len - ke,
                            "clipped_start": ke,
                            "clipped_end": old_len,
                            "clipped_len": len(clip),
                            "note": p.get("note", ""),
                        }
                    )
                keep_r = min(keep_r, ke)
            total_trim += tb
            sides.append(side)
            notes.append(p.get("note", ""))
        if keep_l >= keep_r:
            records.append(
                {
                    "scaffold": name,
                    "side": ",".join(sides) or "NA",
                    "status": "FAIL",
                    "trim_bp": total_trim,
                    "old_len": old_len,
                    "new_len": old_len,
                    "note": "invalid_keep_range",
                }
            )
            continue
        new_seq = seq[keep_l:keep_r]
        seqs[name] = new_seq
        n_ok += 1
        records.append(
            {
                "scaffold": name,
                "side": ",".join(sides),
                "status": "OK",
                "trim_bp": old_len - len(new_seq),
                "old_len": old_len,
                "new_len": len(new_seq),
                "note": ";".join(notes),
            }
        )

    write_fasta(Path(args.out_fa), seqs)
    with Path(args.record).open("w", newline="") as fh:
        fields = ["scaffold", "side", "status", "trim_bp", "old_len", "new_len", "note"]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in records:
            w.writerow(r)
    if out_clipped_fa:
        Path(out_clipped_fa).write_text("".join(clipped_fa_lines))
    if out_clipped_tsv:
        cfields = [
            "scaffold",
            "side",
            "trim_bp",
            "clipped_start",
            "clipped_end",
            "clipped_len",
            "note",
        ]
        with Path(out_clipped_tsv).open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cfields, delimiter="\t")
            w.writeheader()
            for r in clipped_recs:
                w.writerow(r)
    print(f"[trim-buried] ok={n_ok}/{len(by_scaf)} -> {args.out_fa}")
    if clipped_recs:
        print(f"[trim-buried] clipped_ends={len(clipped_recs)}")


def cmd_accept_trial_trim(args: argparse.Namespace) -> None:
    """Keep trial_trim only if motif tip AND quarTeT side pass; else restore pre-trim."""
    pre = read_fasta(Path(args.pre_fa))
    post = read_fasta(Path(args.post_fa))
    telo = (
        parse_telo_info(Path(args.telo_info))
        if args.telo_info and Path(args.telo_info).is_file()
        else {}
    )
    decide = list(csv.DictReader(Path(args.decide).open(), delimiter="\t"))
    trials = [r for r in decide if r.get("decision") == "trial_trim"]

    accept_rows = []
    reject_rows = []
    restore_names: set[str] = set()

    for r in trials:
        name, side = r["scaffold"], r["side"]
        if name not in post:
            reject_rows.append(
                {
                    "scaffold": name,
                    "side": side,
                    "status": "unresolved",
                    "span_reads": r.get("span_reads", ""),
                    "distal_reads": r.get("distal_reads", ""),
                    "extend_reads": r.get("extend_reads", ""),
                    "trim_bp": r.get("trim_bp", ""),
                    "note": "trial_trim_scaffold_missing_post",
                }
            )
            continue
        seq = post[name]
        lc, rc = tip_copies(seq, args.tip_bp)
        tip_c = lc if side == "left" else rc
        t = telo.get(name, {})
        q_ok = bool(t.get("has_left") if side == "left" else t.get("has_right"))
        motif_ok = tip_c >= args.min_copies
        if motif_ok and q_ok:
            accept_rows.append(
                {
                    "scaffold": name,
                    "side": side,
                    "status": "accepted",
                    "tip_copies": tip_c,
                    "quartet_ok": 1,
                    "trim_bp": r.get("trim_bp", ""),
                    "note": r.get("note", ""),
                }
            )
        else:
            restore_names.add(name)
            reject_rows.append(
                {
                    "scaffold": name,
                    "side": side,
                    "status": "unresolved",
                    "span_reads": r.get("span_reads", ""),
                    "distal_reads": r.get("distal_reads", ""),
                    "extend_reads": r.get("extend_reads", ""),
                    "trim_bp": r.get("trim_bp", ""),
                    "note": (
                        f"trial_trim_reject;motif={tip_c}"
                        f"{'PASS' if motif_ok else 'FAIL'};"
                        f"quartet={'PASS' if q_ok else 'FAIL'};"
                        f"{r.get('note', '')}"
                    ),
                }
            )

    for name in restore_names:
        if name in pre:
            post[name] = pre[name]

    # Keep only accepted clipped fragments
    if args.clipped_tsv and Path(args.clipped_tsv).is_file():
        rej_sides = {(x["scaffold"], x["side"]) for x in reject_rows}
        kept_clip = []
        with Path(args.clipped_tsv).open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            fields = list(reader.fieldnames or [])
            for row in reader:
                if row.get("scaffold") in restore_names:
                    continue
                if (row.get("scaffold"), row.get("side")) in rej_sides:
                    continue
                kept_clip.append(row)
        out_ct = Path(args.out_clipped_tsv or args.clipped_tsv)
        if fields:
            with out_ct.open("w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
                w.writeheader()
                for row in kept_clip:
                    w.writerow(row)

    if args.clipped_fa and Path(args.clipped_fa).is_file() and restore_names:
        # rewrite clipped fa without restored scaffolds
        kept_fa: list[str] = []
        cur_hdr = None
        cur_seq: list[str] = []
        drop = False

        def flush() -> None:
            nonlocal cur_hdr, cur_seq, drop
            if cur_hdr and not drop:
                kept_fa.append(cur_hdr + "\n" + "".join(cur_seq) + "\n")
            cur_hdr, cur_seq, drop = None, [], False

        with Path(args.clipped_fa).open() as fh:
            for line in fh:
                if line.startswith(">"):
                    flush()
                    cur_hdr = line.rstrip("\n")
                    scaf = cur_hdr[1:].split("__", 1)[0]
                    drop = scaf in restore_names
                else:
                    cur_seq.append(line.rstrip("\n"))
            flush()
        Path(args.out_clipped_fa or args.clipped_fa).write_text("".join(kept_fa))

    write_fasta(Path(args.out_fa), post)

    af = [
        "scaffold",
        "side",
        "status",
        "tip_copies",
        "quartet_ok",
        "trim_bp",
        "note",
    ]
    with Path(args.out_accept).open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=af, delimiter="\t")
        w.writeheader()
        for row in accept_rows:
            w.writerow(row)

    ufields = [
        "scaffold",
        "side",
        "status",
        "span_reads",
        "distal_reads",
        "extend_reads",
        "trim_bp",
        "note",
    ]
    out_u = Path(args.out_unresolved)
    existing = []
    if out_u.is_file() and out_u.stat().st_size > 0:
        with out_u.open() as fh:
            existing = list(csv.DictReader(fh, delimiter="\t"))
    with out_u.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=ufields, delimiter="\t")
        w.writeheader()
        for row in existing:
            w.writerow({k: row.get(k, "") for k in ufields})
        for row in reject_rows:
            w.writerow(row)

    print(
        f"[accept-trial-trim] accepted={len(accept_rows)} "
        f"rejected={len(reject_rows)} restored_scaffolds={len(restore_names)} "
        f"-> {args.out_fa}"
    )
    for row in accept_rows:
        print(
            f"  ACCEPT {row['scaffold']}\t{row['side']}\ttip={row['tip_copies']}\t"
            f"trim_bp={row['trim_bp']}"
        )
    for row in reject_rows:
        print(f"  REJECT {row['scaffold']}\t{row['side']}\t{row['note']}")


def cmd_apply_revcomp(args: argparse.Namespace) -> None:
    seqs = read_fasta(Path(args.fasta))
    names = [
        ln.strip()
        for ln in Path(args.revcomp_list).read_text().splitlines()
        if ln.strip() and not ln.startswith("#")
    ]
    n = 0
    for name in names:
        if name not in seqs:
            print(f"[revcomp] WARN missing scaffold: {name}", file=sys.stderr)
            continue
        seqs[name] = revcomp(seqs[name])
        n += 1
    write_fasta(Path(args.out_fa), seqs)
    print(f"[revcomp] flipped={n} -> {args.out_fa}")


def _motif_both_miss(
    seqs: dict[str, str], top: list[str], tip_bp: int, min_copies: int
) -> tuple[int, list[str]]:
    both = 0
    miss = []
    for name in top:
        lc, rc = tip_copies(seqs[name], tip_bp)
        if lc >= min_copies and rc >= min_copies:
            both += 1
        else:
            miss.append(name)
    return both, miss


def _quartet_both_miss(telo: dict, top: list[str]) -> tuple[int, list[str]]:
    both = 0
    miss = []
    for name in top:
        t = telo.get(name, {})
        if t.get("has_left") and t.get("has_right"):
            both += 1
        else:
            miss.append(name)
    return both, miss


def cmd_check_perfect(args: argparse.Namespace) -> None:
    gaps = parse_gaps_tsv(Path(args.gaps_tsv)) if Path(args.gaps_tsv).is_file() else []
    gap_n = len(gaps)
    seqs = read_fasta(Path(args.fasta))
    top = resolve_top(seqs, args.top_list, args.top_n)
    telo_for_unres = (
        parse_telo_info(Path(args.telo_info))
        if args.telo_info and Path(args.telo_info).is_file()
        else {}
    )

    unresolved_n = 0
    unresolved_names: list[str] = []
    if args.unresolved and Path(args.unresolved).is_file():
        with Path(args.unresolved).open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                name = row.get("scaffold", "")
                side = row.get("side", "")
                if not name or name not in seqs:
                    unresolved_n += 1
                    unresolved_names.append(f"{name}:{side}")
                    continue
                # quartet 模式：该侧已有 quarTeT 端粒则视为已解决
                if args.mode == "quartet" and telo_for_unres:
                    t = telo_for_unres.get(name, {})
                    ok_side = (
                        t.get("has_left") if side == "left" else t.get("has_right")
                    )
                    if ok_side:
                        continue
                elif args.mode != "quartet":
                    lc, rc = tip_copies(seqs[name], args.tip_bp)
                    tip_c = lc if side == "left" else rc
                    if tip_c >= args.min_copies:
                        continue
                unresolved_n += 1
                unresolved_names.append(f"{name}:{side}")

    if args.mode == "motif":
        both, miss = _motif_both_miss(seqs, top, args.tip_bp, args.min_copies)
        ok = both == len(top) and gap_n == 0 and unresolved_n == 0
        print(
            f"[check-perfect:motif] both_telo={both}/{len(top)} gaps={gap_n} "
            f"unresolved={unresolved_n} min_copies={args.min_copies} tip_bp={args.tip_bp} "
            f"{'PASS' if ok else 'FAIL'} missing={','.join(miss) if miss else 'none'}"
        )
        sys.exit(0 if ok else 1)

    if args.mode == "quartet":
        if not args.telo_info or not Path(args.telo_info).is_file():
            print("[check-perfect:quartet] FAIL missing --telo-info", file=sys.stderr)
            sys.exit(1)
        telo = parse_telo_info(Path(args.telo_info))
        both, miss = _quartet_both_miss(telo, top)
        ok = both == len(top) and gap_n == 0 and unresolved_n == 0
        print(
            f"[check-perfect:quartet] both_telo={both}/{len(top)} gaps={gap_n} "
            f"unresolved={unresolved_n} "
            f"{'PASS' if ok else 'FAIL'} missing={','.join(miss) if miss else 'none'}"
        )
        sys.exit(0 if ok else 1)

    # both / either: combine motif and quarTeT
    m_both, m_miss = _motif_both_miss(seqs, top, args.tip_bp, args.min_copies)
    motif_ok = m_both == len(top)
    quartet_ok = False
    q_both, q_miss = 0, list(top)
    if args.telo_info and Path(args.telo_info).is_file():
        telo = parse_telo_info(Path(args.telo_info))
        q_both, q_miss = _quartet_both_miss(telo, top)
        quartet_ok = q_both == len(top)
    elif args.mode == "both":
        print("[check-perfect:both] FAIL missing --telo-info", file=sys.stderr)
        sys.exit(1)

    if args.mode == "both":
        ok = gap_n == 0 and unresolved_n == 0 and motif_ok and quartet_ok
        tag = "both"
        combo = "AND"
    else:
        ok = gap_n == 0 and unresolved_n == 0 and (motif_ok or quartet_ok)
        tag = "either"
        combo = "OR"

    print(
        f"[check-perfect:{tag}] gaps={gap_n} unresolved={unresolved_n} "
        f"[{','.join(unresolved_names) if unresolved_names else 'none'}] "
        f"motif={m_both}/{len(top)}{'PASS' if motif_ok else 'FAIL'} "
        f"missing_motif={','.join(m_miss) if m_miss else 'none'} {combo} "
        f"quartet={q_both}/{len(top)}{'PASS' if quartet_ok else 'FAIL'} "
        f"missing_quartet={','.join(q_miss) if q_miss else 'none'} | "
        f"{'PASS' if ok else 'FAIL'}"
    )
    sys.exit(0 if ok else 1)


def cmd_merge_fill_records(args: argparse.Namespace) -> None:
    """Merge base (06) + new (final_fix) fill records for 07 polish.

    Dedup:
      telo key = scaffold + side  → prefer new over base
      gap  key = scaffold + gap_index + old_start + old_end → prefer new
    Rows from either side are kept; conflicting keys keep the newer row.
    """
    if args.kind == "telo":
        fields = [
            "scaffold",
            "side",
            "status",
            "read",
            "added_bp",
            "telo_copies",
            "ident",
            "note",
        ]

        def key_fn(r: dict) -> tuple:
            return (r.get("scaffold", ""), r.get("side", ""))

    else:
        fields = [
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
        ]

        def key_fn(r: dict) -> tuple:
            return (
                r.get("scaffold", ""),
                r.get("gap_index", ""),
                r.get("old_start", ""),
                r.get("old_end", ""),
            )

    def load_rows(path: Path, source: str) -> list[dict]:
        if not path.is_file() or path.stat().st_size == 0:
            return []
        rows = []
        with path.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for r in reader:
                if not r.get("scaffold"):
                    continue
                note = r.get("note") or ""
                if "source=" not in note:
                    note = f"source={source};{note}" if note else f"source={source}"
                out = {k: r.get(k, "") for k in fields}
                out["note"] = note
                rows.append(out)
        return rows

    base_rows = load_rows(Path(args.base), "06")
    new_rows = load_rows(Path(args.new), "final_fix")
    merged: dict[tuple, dict] = {}
    order: list[tuple] = []
    for r in base_rows:
        k = key_fn(r)
        if k not in merged:
            order.append(k)
        merged[k] = r
    for r in new_rows:
        k = key_fn(r)
        if k not in merged:
            order.append(k)
        merged[k] = r  # prefer final_fix

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for k in order:
            w.writerow(merged[k])
    n_ok = sum(1 for k in order if merged[k].get("status") == "OK")
    print(
        f"[merge-fill-records:{args.kind}] base={len(base_rows)} "
        f"new={len(new_rows)} merged={len(order)} ok={n_ok} -> {out}"
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    sp = ap.add_subparsers(dest="cmd", required=True)

    p = sp.add_parser("call-motif")
    p.add_argument("--fasta", required=True)
    p.add_argument("--top-list", default="")
    p.add_argument("--top-n", type=int, default=20)
    p.add_argument("--tip-bp", type=int, default=20000)
    p.add_argument("--min-copies", type=int, default=50)
    p.add_argument("--gaps-tsv", default="")
    p.add_argument("--out-tsv", required=True)
    p.add_argument("--out-missing", required=True)
    p.add_argument("--out-telo-info", required=True)
    p.add_argument("--out-top-list", required=True)
    p.add_argument("--out-summary", required=True)
    p.set_defaults(func=cmd_call_motif)

    p = sp.add_parser("diagnose-orient")
    p.add_argument("--fasta", required=True)
    p.add_argument("--telo-info", default="")  # optional quarTeT for comparison
    p.add_argument("--missing-telo", default="")  # preferred: quarTeT missing
    p.add_argument("--top-list", default="")
    p.add_argument("--top-n", type=int, default=20)
    p.add_argument("--out-tsv", required=True)
    p.add_argument("--revcomp-list", required=True)
    p.add_argument("--trim-plan", required=True)
    p.add_argument("--fill-missing", required=True)
    p.add_argument(
        "--trust-missing",
        action="store_true",
        help="缺失列表来自 quarTeT：不因 tip motif 达标而 skip",
    )
    p.add_argument("--tip-bp", type=int, default=20000)
    p.add_argument("--tip-min-copies", type=int, default=50)
    p.add_argument("--internal-min-copies", type=int, default=50)
    p.add_argument("--near-tip-bp", type=int, default=500000)
    p.add_argument("--max-trim-bp", type=int, default=500000)
    p.add_argument("--flip-onesided", action="store_true")
    p.set_defaults(func=cmd_diagnose_orient)

    p = sp.add_parser("make-buried-baits")
    p.add_argument("--fasta", required=True)
    p.add_argument("--candidates", required=True)
    p.add_argument("--out-fa", required=True)
    p.add_argument("--out-coords", required=True)
    p.add_argument("--anchor-bp", type=int, default=50000)
    p.set_defaults(func=cmd_make_buried_baits)

    p = sp.add_parser("decide-buried-ont")
    p.add_argument("--candidates", required=True)
    p.add_argument("--coords", required=True)
    p.add_argument("--paf", required=True)
    p.add_argument("--reads-fa", required=True)
    p.add_argument("--out-trim-plan", required=True)
    p.add_argument("--out-extend", required=True)
    p.add_argument("--out-unresolved", required=True)
    p.add_argument("--out-decide", required=True)
    p.add_argument("--min-alen", type=int, default=500)
    p.add_argument("--min-ident", type=float, default=0.55)
    p.add_argument("--min-support-reads", type=int, default=3)
    p.add_argument("--min-extend-reads", type=int, default=1)
    p.add_argument("--min-overhang", type=int, default=50)
    p.add_argument("--min-telo-copies", type=int, default=5)
    p.add_argument(
        "--auto-trim",
        type=int,
        default=1,
        choices=(0, 1),
        help="1=trim when no ONT support (default); 0=mark unresolved instead",
    )
    p.add_argument(
        "--allow-distal-trial-trim",
        type=int,
        default=1,
        choices=(0, 1),
        help="1=distal_supported 时试 trim（需 new_tip 达标）",
    )
    p.add_argument(
        "--min-new-tip-copies",
        type=int,
        default=50,
        help="trial_trim 要求预测新 tip motif copies",
    )
    p.add_argument("--max-trim-bp", type=int, default=500000)
    p.set_defaults(func=cmd_decide_buried_ont)

    p = sp.add_parser("trim-buried")
    p.add_argument("--fasta", required=True)
    p.add_argument("--trim-plan", required=True)
    p.add_argument("--out-fa", required=True)
    p.add_argument("--record", required=True)
    p.add_argument("--max-trim-bp", type=int, default=500000)
    p.add_argument("--out-clipped-fa", default="")
    p.add_argument("--out-clipped-tsv", default="")
    p.set_defaults(func=cmd_trim_buried)

    p = sp.add_parser("accept-trial-trim")
    p.add_argument("--pre-fa", required=True)
    p.add_argument("--post-fa", required=True)
    p.add_argument("--decide", required=True)
    p.add_argument("--telo-info", required=True)
    p.add_argument("--out-fa", required=True)
    p.add_argument("--out-accept", required=True)
    p.add_argument("--out-unresolved", required=True)
    p.add_argument("--tip-bp", type=int, default=20000)
    p.add_argument("--min-copies", type=int, default=50)
    p.add_argument("--clipped-tsv", default="")
    p.add_argument("--clipped-fa", default="")
    p.add_argument("--out-clipped-tsv", default="")
    p.add_argument("--out-clipped-fa", default="")
    p.set_defaults(func=cmd_accept_trial_trim)

    p = sp.add_parser("apply-revcomp")
    p.add_argument("--fasta", required=True)
    p.add_argument("--revcomp-list", required=True)
    p.add_argument("--out-fa", required=True)
    p.set_defaults(func=cmd_apply_revcomp)

    p = sp.add_parser("check-perfect")
    p.add_argument("--fasta", required=True)
    p.add_argument("--gaps-tsv", required=True)
    p.add_argument("--top-list", default="")
    p.add_argument("--top-n", type=int, default=20)
    p.add_argument(
        "--mode",
        choices=("both", "either", "motif", "quartet"),
        default="quartet",
        help="quartet=quarTeT only (default); both/either/motif kept for debug",
    )
    p.add_argument("--telo-info", default="")
    p.add_argument("--unresolved", default="")
    p.add_argument("--tip-bp", type=int, default=20000)
    p.add_argument("--min-copies", type=int, default=50)
    p.set_defaults(func=cmd_check_perfect)

    p = sp.add_parser("merge-fill-records")
    p.add_argument("--kind", required=True, choices=("telo", "gap"))
    p.add_argument("--base", required=True, help="原 06/04_final 记录")
    p.add_argument("--new", required=True, help="final_fix 本次记录")
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_merge_fill_records)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
