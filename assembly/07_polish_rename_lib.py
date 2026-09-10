#!/usr/bin/env python3
"""Helpers for 07_polish_rename: fill-window polish patch + Chr rename/orient."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path


def revcomp(seq: str) -> str:
    trans = str.maketrans(
        "ACGTRYMKBDHVacgtrymkbdhvNn",
        "TGCAYRKMVHDBtgcayrkmvhdbNn",
    )
    return seq.translate(trans)[::-1]


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


def write_fasta_record(handle, name: str, seq: str, width: int = 80) -> None:
    handle.write(f">{name}\n")
    for i in range(0, len(seq), width):
        handle.write(seq[i : i + width] + "\n")


def parse_scaffold_number(name: str) -> int | None:
    m = re.fullmatch(r"scaffold_(\d+)", name)
    return int(m.group(1)) if m else None


def parse_chr_number(name: str) -> int | None:
    m = re.search(r"Chr0*(\d+)", name, re.IGNORECASE)
    if not m:
        return None
    n = int(m.group(1))
    return n if 1 <= n <= 20 else None


def load_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def collect_windows(
    seqs: dict[str, str],
    telo_tsv: Path,
    gap_tsv: Path,
    flank: int,
) -> list[dict]:
    """Build 1-based inclusive windows on the FINAL assembly for OK fills."""
    windows: list[dict] = []

    for row in load_tsv(telo_tsv):
        if row.get("status") != "OK":
            continue
        name = row["scaffold"]
        side = row["side"]
        added = int(float(row.get("added_bp") or 0))
        if added <= 0 or name not in seqs:
            continue
        L = len(seqs[name])
        if side == "left":
            start = 1
            end = min(L, added + flank)
        else:
            start = max(1, L - added - flank + 1)
            end = L
        wid = f"{name}__telo_{side}_{start}_{end}"
        windows.append(
            {
                "id": wid,
                "scaffold": name,
                "kind": f"telo_{side}",
                "start": start,
                "end": end,
                "note": f"added_bp={added};read={row.get('read','')}",
            }
        )

    # Gap fills: telo only extends ends, so internal old_start is still valid on final
    # if a single gap per scaffold (true for KD35). Multiple gaps: apply descending
    # was done at fill time; old_start in record is pre-fill on after_telo coords —
    # with one gap/scaffold this matches final.
    for row in load_tsv(gap_tsv):
        if row.get("status") != "OK":
            continue
        name = row["scaffold"]
        fill_bp = int(float(row.get("fill_bp") or 0))
        old_start = int(float(row.get("old_start") or 0))
        if fill_bp <= 0 or old_start <= 0 or name not in seqs:
            continue
        L = len(seqs[name])
        start = max(1, old_start - flank)
        end = min(L, old_start + fill_bp - 1 + flank)
        if end < start:
            continue
        gidx = row.get("gap_index", "0")
        wid = f"{name}__gap{gidx}_{start}_{end}"
        windows.append(
            {
                "id": wid,
                "scaffold": name,
                "kind": f"gap{gidx}",
                "start": start,
                "end": end,
                "note": f"fill_bp={fill_bp};read={row.get('read','')}",
            }
        )

    return windows


def cmd_make_windows(args: argparse.Namespace) -> None:
    seqs = read_fasta(Path(args.fasta))
    windows = collect_windows(
        seqs, Path(args.telo_record), Path(args.gap_record), args.flank
    )
    out_fa = Path(args.out_fa)
    out_bed = Path(args.out_bed)
    with out_fa.open("w") as fa, out_bed.open("w") as bed:
        bed.write("window_id\tscaffold\tstart\tend\tkind\tnote\n")
        for w in windows:
            seq = seqs[w["scaffold"]][w["start"] - 1 : w["end"]]
            fa.write(f">{w['id']}\n")
            for i in range(0, len(seq), 80):
                fa.write(seq[i : i + 80] + "\n")
            bed.write(
                f"{w['id']}\t{w['scaffold']}\t{w['start']}\t{w['end']}\t{w['kind']}\t{w['note']}\n"
            )
    print(f"[make-windows] n={len(windows)} -> {out_fa} ; {out_bed}")


def cmd_apply_windows(args: argparse.Namespace) -> None:
    seqs = read_fasta(Path(args.fasta))
    polished = read_fasta(Path(args.polished_windows))
    bed_rows = load_tsv(Path(args.windows_bed))
    record_path = Path(args.record)
    out = dict(seqs)
    rows = []

    # apply from scaffold end so coordinates stay valid within a scaffold
    by_scaff: dict[str, list] = defaultdict(list)
    for r in bed_rows:
        by_scaff[r["scaffold"]].append(r)

    for scaff, glist in by_scaff.items():
        seq = out[scaff]
        for r in sorted(glist, key=lambda x: -int(x["start"])):
            wid = r["window_id"]
            start = int(r["start"])
            end = int(r["end"])
            old = seq[start - 1 : end]
            if wid not in polished:
                rows.append(
                    {
                        "window_id": wid,
                        "scaffold": scaff,
                        "start": start,
                        "end": end,
                        "old_len": len(old),
                        "new_len": len(old),
                        "status": "KEEP",
                        "note": "missing_in_racon_output",
                    }
                )
                continue
            new = polished[wid].upper()
            # racon may change length slightly; allow modest delta
            if abs(len(new) - len(old)) > max(500, int(0.05 * len(old)) + 1):
                rows.append(
                    {
                        "window_id": wid,
                        "scaffold": scaff,
                        "start": start,
                        "end": end,
                        "old_len": len(old),
                        "new_len": len(new),
                        "status": "SKIP_LEN",
                        "note": "length_delta_too_large",
                    }
                )
                continue
            seq = seq[: start - 1] + new + seq[end:]
            rows.append(
                {
                    "window_id": wid,
                    "scaffold": scaff,
                    "start": start,
                    "end": end,
                    "old_len": len(old),
                    "new_len": len(new),
                    "status": "OK",
                    "note": "racon",
                }
            )
        out[scaff] = seq

    # preserve original key order
    ordered = {k: out[k] for k in seqs}
    write_fasta(Path(args.out_fa), ordered)
    with record_path.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=[
                "window_id",
                "scaffold",
                "start",
                "end",
                "old_len",
                "new_len",
                "status",
                "note",
            ],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(rows)
    ok = sum(1 for r in rows if r["status"] == "OK")
    print(f"[apply-windows] ok={ok}/{len(rows)} -> {args.out_fa}")


def accumulate_support(
    paf_path: Path,
    seqs: dict[str, str],
    min_mapq: int,
    min_block: int,
    restrict_to_first_20: bool,
) -> dict:
    support: dict = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    with paf_path.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            qname = parts[0]
            strand = parts[4]
            tname = parts[5]
            try:
                alen = int(parts[10])
                mapq = int(parts[11])
            except ValueError:
                continue
            if qname not in seqs:
                continue
            sn = parse_scaffold_number(qname)
            if sn is None:
                continue
            if restrict_to_first_20 and sn > 20:
                continue
            cn = parse_chr_number(tname)
            if cn is None or strand not in "+-":
                continue
            if mapq < min_mapq or alen < min_block:
                continue
            tags = parts[12:]
            tp = [t for t in tags if t.startswith("tp:A:")]
            if tp and "tp:A:P" not in tp:
                continue
            support[qname][cn][strand] += alen
    return support


def judge_mapping(
    seqs: dict[str, str],
    support: dict,
    restrict_to_first_20: bool,
    min_best_chr_ratio: float,
    min_strand_ratio: float,
) -> list[dict]:
    rows = []
    names = []
    for qname in seqs:
        sn = parse_scaffold_number(qname)
        if sn is None:
            continue
        if restrict_to_first_20 and sn > 20:
            continue
        names.append(qname)
    names.sort(key=lambda x: parse_scaffold_number(x) or 0)

    for qname in names:
        if qname not in support:
            rows.append(
                {
                    "old_name": qname,
                    "new_name": "",
                    "chr_num": "",
                    "strand": "",
                    "best_chr_support": 0,
                    "second_chr_support": 0,
                    "total_support": 0,
                    "best_chr_ratio": 0.0,
                    "plus_support": 0,
                    "minus_support": 0,
                    "strand_ratio": 0.0,
                    "status": "NO_HIT",
                }
            )
            continue
        chr_support = {
            c: support[qname][c].get("+", 0) + support[qname][c].get("-", 0)
            for c in support[qname]
        }
        total = sum(chr_support.values())
        ranked = sorted(chr_support.items(), key=lambda x: -x[1])
        best_c, best_s = ranked[0]
        second_s = ranked[1][1] if len(ranked) > 1 else 0
        ratio = best_s / total if total else 0.0
        plus = support[qname][best_c].get("+", 0)
        minus = support[qname][best_c].get("-", 0)
        strand = "+" if plus >= minus else "-"
        strand_s = plus if strand == "+" else minus
        strand_ratio = strand_s / best_s if best_s else 0.0
        status = "PASS"
        if ratio < min_best_chr_ratio:
            status = "AMBIGUOUS_CHROMOSOME"
        if strand_ratio < min_strand_ratio:
            status = "AMBIGUOUS_ORIENTATION"
        rows.append(
            {
                "old_name": qname,
                "new_name": f"Chr{best_c:02d}",
                "chr_num": best_c,
                "strand": strand,
                "best_chr_support": best_s,
                "second_chr_support": second_s,
                "total_support": total,
                "best_chr_ratio": ratio,
                "plus_support": plus,
                "minus_support": minus,
                "strand_ratio": strand_ratio,
                "status": status,
            }
        )
    return rows


def write_mapping_tsv(path: Path, rows: list[dict]) -> None:
    fields = [
        "old_name",
        "new_name",
        "chr_num",
        "strand",
        "best_chr_support",
        "second_chr_support",
        "total_support",
        "best_chr_ratio",
        "plus_support",
        "minus_support",
        "strand_ratio",
        "status",
    ]
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in rows:
            out = dict(r)
            out["best_chr_ratio"] = f'{float(r["best_chr_ratio"]):.6f}'
            out["strand_ratio"] = f'{float(r["strand_ratio"]):.6f}'
            w.writerow(out)


def cmd_rename(args: argparse.Namespace) -> None:
    seqs = read_fasta(Path(args.fasta))
    primary = accumulate_support(
        Path(args.primary_paf),
        seqs,
        args.min_mapq,
        args.min_block,
        args.restrict_to_first_20,
    )
    rows = judge_mapping(
        seqs,
        primary,
        args.restrict_to_first_20,
        args.min_best_chr_ratio,
        args.min_strand_ratio,
    )
    write_mapping_tsv(Path(args.mapping_out), rows)

    # ZH13 check
    if args.check_paf:
        check = accumulate_support(
            Path(args.check_paf),
            seqs,
            args.min_mapq,
            args.min_block,
            args.restrict_to_first_20,
        )
        check_rows = judge_mapping(
            seqs,
            check,
            args.restrict_to_first_20,
            args.min_best_chr_ratio,
            args.min_strand_ratio,
        )
        write_mapping_tsv(Path(args.check_out), check_rows)
        primary_pass = {r["old_name"]: r for r in rows if r["status"] == "PASS"}
        check_pass = {r["old_name"]: r for r in check_rows if r["status"] == "PASS"}
        with Path(args.check_report).open("w") as fh:
            fh.write("old_name\twm82_chr\twm82_strand\tzh13_chr\tzh13_strand\tconcordant\n")
            for name in sorted(primary_pass, key=lambda x: parse_scaffold_number(x) or 0):
                a = primary_pass[name]
                b = check_pass.get(name)
                if b is None:
                    fh.write(
                        f"{name}\t{a['new_name']}\t{a['strand']}\tNA\tNA\tNO_ZH13_PASS\n"
                    )
                    continue
                ok = a["chr_num"] == b["chr_num"] and a["strand"] == b["strand"]
                fh.write(
                    f"{name}\t{a['new_name']}\t{a['strand']}\t{b['new_name']}\t{b['strand']}\t"
                    f"{'YES' if ok else 'CONFLICT'}\n"
                )
                if not ok:
                    print(
                        f"[zh13-check] CONFLICT {name}: "
                        f"Wm82 {a['new_name']}{a['strand']} vs ZH13 {b['new_name']}{b['strand']}",
                        file=sys.stderr,
                    )

    best = {r["old_name"]: r for r in rows if r["status"] == "PASS"}
    chr_to = defaultdict(list)
    for name, r in best.items():
        chr_to[int(r["chr_num"])].append(name)

    dups = {c: v for c, v in chr_to.items() if len(v) > 1}
    missing = [c for c in range(1, 21) if c not in chr_to]
    if dups:
        for c, names in sorted(dups.items()):
            print(f"ERROR duplicate Chr{c:02d}: {', '.join(names)}", file=sys.stderr)
        sys.exit(1)
    if args.require_complete_20 and (len(best) != 20 or missing):
        print(
            f"ERROR incomplete mapping: pass={len(best)} missing={missing}",
            file=sys.stderr,
        )
        sys.exit(1)

    used = set()
    main = {}
    for name, r in best.items():
        seq = seqs[name]
        if r["strand"] == "-":
            seq = revcomp(seq)
        main[int(r["chr_num"])] = (r["new_name"], seq, name)
        used.add(name)

    unanchored = [n for n in seqs if n not in used]

    def ukey(n: str):
        sn = parse_scaffold_number(n)
        return (0, sn, n) if sn is not None else (1, 0, n)

    unanchored.sort(key=ukey)

    in_len = sum(len(s) for s in seqs.values())
    out_len = 0
    with Path(args.out_fa).open("w") as fout:
        for c in sorted(main):
            new_name, seq, _ = main[c]
            write_fasta_record(fout, new_name, seq)
            out_len += len(seq)
        for i, old in enumerate(unanchored, start=1):
            new_name = f"Unanchored_{i:04d}_original_{old}"
            write_fasta_record(fout, new_name, seqs[old])
            out_len += len(seqs[old])

    if out_len != in_len:
        sys.exit(f"ERROR length mismatch in={in_len} out={out_len}")
    print(
        f"[rename] Chr={len(main)} unanchored={len(unanchored)} -> {args.out_fa}"
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)

    p = sp.add_parser("make-windows")
    p.add_argument("--fasta", required=True)
    p.add_argument("--telo-record", required=True)
    p.add_argument("--gap-record", required=True)
    p.add_argument("--flank", type=int, default=20000)
    p.add_argument("--out-fa", required=True)
    p.add_argument("--out-bed", required=True)
    p.set_defaults(func=cmd_make_windows)

    p = sp.add_parser("apply-windows")
    p.add_argument("--fasta", required=True)
    p.add_argument("--polished-windows", required=True)
    p.add_argument("--windows-bed", required=True)
    p.add_argument("--out-fa", required=True)
    p.add_argument("--record", required=True)
    p.set_defaults(func=cmd_apply_windows)

    p = sp.add_parser("rename")
    p.add_argument("--fasta", required=True)
    p.add_argument("--primary-paf", required=True)
    p.add_argument("--check-paf", default="")
    p.add_argument("--out-fa", required=True)
    p.add_argument("--mapping-out", required=True)
    p.add_argument("--check-out", default="")
    p.add_argument("--check-report", default="")
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--min-block", type=int, default=100000)
    p.add_argument("--min-best-chr-ratio", type=float, default=0.60)
    p.add_argument("--min-strand-ratio", type=float, default=0.60)
    p.add_argument("--restrict-to-first-20", action="store_true", default=True)
    p.add_argument("--no-restrict-to-first-20", action="store_false", dest="restrict_to_first_20")
    p.add_argument("--require-complete-20", action="store_true", default=True)
    p.add_argument("--no-require-complete-20", action="store_false", dest="require_complete_20")
    p.set_defaults(func=cmd_rename)

    args = ap.parse_args()
    if args.cmd == "rename":
        if args.check_paf and not args.check_out:
            args.check_out = str(Path(args.mapping_out).with_name("zh13_mapping.tsv"))
        if args.check_paf and not args.check_report:
            args.check_report = str(Path(args.mapping_out).with_name("zh13_check.tsv"))
    args.func(args)


if __name__ == "__main__":
    main()
