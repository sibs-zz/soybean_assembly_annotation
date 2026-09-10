#!/usr/bin/env python3
"""Helpers for 08_fa_qc: PAF/coords dotplots + summary aggregation."""
from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path


def parse_chr_key(name: str) -> tuple:
    m = re.search(r"(?:Chr|chr|CHR)[_\s.-]*0*(\d+)", name)
    if m:
        return (0, int(m.group(1)), name)
    return (1, 0, name)


def load_fai_lengths(fai: Path) -> dict[str, int]:
    lens: dict[str, int] = {}
    if not fai.is_file():
        return lens
    with fai.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                lens[parts[0]] = int(parts[1])
    return lens


def read_paf(path: Path, min_alen: int = 10000, min_mapq: int = 20):
    rows = []
    with path.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 12:
                continue
            try:
                qname, qlen, qs, qe = p[0], int(p[1]), int(p[2]), int(p[3])
                strand = p[4]
                tname, tlen, ts, te = p[5], int(p[6]), int(p[7]), int(p[8])
                alen = int(p[10])
                mapq = int(p[11])
            except ValueError:
                continue
            if alen < min_alen or mapq < min_mapq:
                continue
            tags = p[12:]
            tp = [t for t in tags if t.startswith("tp:A:")]
            if tp and "tp:A:P" not in tp:
                continue
            rows.append(
                {
                    "q": qname,
                    "qlen": qlen,
                    "qs": qs,
                    "qe": qe,
                    "strand": strand,
                    "t": tname,
                    "tlen": tlen,
                    "ts": ts,
                    "te": te,
                    "alen": alen,
                }
            )
    return rows


def read_coords(path: Path, min_alen: int = 10000):
    """Parse show-coords -T -H output (tab, no header)."""
    rows = []
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("[") or line.startswith("="):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            try:
                # S1 E1 S2 E2 LEN1 LEN2 %IDY TAG1?  or with -r: ref query
                # -T H without -r: S1 E1 S2 E2 LEN1 LEN2 %IDY ref query
                s1, e1, s2, e2 = int(p[0]), int(p[1]), int(p[2]), int(p[3])
                len1, len2 = int(p[4]), int(p[5])
                tname, qname = p[-2], p[-1]
            except ValueError:
                continue
            alen = min(len1, len2)
            if alen < min_alen:
                continue
            strand = "+" if (s2 <= e2) == (s1 <= e1) or (s2 < e2) else "-"
            # normalize
            ts, te = (s1, e1) if s1 <= e1 else (e1, s1)
            qs, qe = (s2, e2) if s2 <= e2 else (e2, s2)
            if s2 > e2:
                strand = "-"
            else:
                strand = "+"
            rows.append(
                {
                    "q": qname,
                    "qlen": 0,
                    "qs": qs,
                    "qe": qe,
                    "strand": strand,
                    "t": tname,
                    "tlen": 0,
                    "ts": ts,
                    "te": te,
                    "alen": alen,
                }
            )
    return rows


def cumulative_offsets(names: list[str], lengths: dict[str, int]) -> dict[str, int]:
    off = {}
    cur = 0
    for n in names:
        off[n] = cur
        cur += lengths.get(n, 0)
    return off


def cmd_dotplot(args: argparse.Namespace) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    paf = Path(args.paf) if args.paf else None
    coords = Path(args.coords) if args.coords else None
    if paf and paf.is_file():
        hits = read_paf(paf, min_alen=args.min_alen, min_mapq=args.min_mapq)
        q_lens = {h["q"]: h["qlen"] for h in hits}
        t_lens = {h["t"]: h["tlen"] for h in hits}
    elif coords and coords.is_file():
        hits = read_coords(coords, min_alen=args.min_alen)
        q_lens = load_fai_lengths(Path(args.query_fai)) if args.query_fai else {}
        t_lens = load_fai_lengths(Path(args.ref_fai)) if args.ref_fai else {}
        for h in hits:
            if h["q"] not in q_lens:
                q_lens[h["q"]] = max(h["qe"], q_lens.get(h["q"], 0))
            if h["t"] not in t_lens:
                t_lens[h["t"]] = max(h["te"], t_lens.get(h["t"], 0))
    else:
        raise SystemExit("Need --paf or --coords")

    if args.query_fai:
        q_lens.update(load_fai_lengths(Path(args.query_fai)))
    if args.ref_fai:
        t_lens.update(load_fai_lengths(Path(args.ref_fai)))

    # keep chromosome-like only
    def is_chr(n: str) -> bool:
        return bool(re.search(r"(?:Chr|chr)[_\s.-]*\d+", n))

    hits = [h for h in hits if is_chr(h["q"]) and is_chr(h["t"])]
    if not hits:
        raise SystemExit("No chromosome hits to plot")

    q_names = sorted({h["q"] for h in hits}, key=parse_chr_key)
    t_names = sorted({h["t"] for h in hits}, key=parse_chr_key)
    q_off = cumulative_offsets(q_names, q_lens)
    t_off = cumulative_offsets(t_names, t_lens)
    q_total = sum(q_lens.get(n, 0) for n in q_names) / 1e6
    t_total = sum(t_lens.get(n, 0) for n in t_names) / 1e6

    segs_fwd, segs_rev = [], []
    for h in hits:
        x0 = (t_off[h["t"]] + h["ts"]) / 1e6
        x1 = (t_off[h["t"]] + h["te"]) / 1e6
        y0 = (q_off[h["q"]] + h["qs"]) / 1e6
        y1 = (q_off[h["q"]] + h["qe"]) / 1e6
        if h["strand"] == "-":
            segs_rev.append([(x0, y1), (x1, y0)])
        else:
            segs_fwd.append([(x0, y0), (x1, y1)])

    fig_w = max(8.0, min(14.0, t_total / 80.0 + 6.0))
    fig_h = max(8.0, min(14.0, q_total / 80.0 + 6.0))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    if segs_fwd:
        ax.add_collection(
            LineCollection(segs_fwd, colors="#1f4e79", linewidths=0.35, alpha=0.85)
        )
    if segs_rev:
        ax.add_collection(
            LineCollection(segs_rev, colors="#b22222", linewidths=0.35, alpha=0.85)
        )

    # chromosome grid
    for n in t_names[1:]:
        ax.axvline(t_off[n] / 1e6, color="#bbbbbb", lw=0.4, zorder=0)
    for n in q_names[1:]:
        ax.axhline(q_off[n] / 1e6, color="#bbbbbb", lw=0.4, zorder=0)

    ax.set_xlim(0, t_total)
    ax.set_ylim(0, q_total)
    ax.set_xlabel(f"{args.ref_label} (Mb)", fontsize=11)
    ax.set_ylabel(f"{args.query_label} (Mb)", fontsize=11)
    ax.set_title(args.title, fontsize=13, pad=10)
    ax.set_aspect("equal", adjustable="box")

    # tick labels at chromosome midpoints
    ax.set_xticks([(t_off[n] + t_lens.get(n, 0) / 2) / 1e6 for n in t_names])
    ax.set_xticklabels(
        [re.sub(r".*?(Chr0*)", "Chr", n) for n in t_names],
        rotation=90,
        fontsize=7,
    )
    ax.set_yticks([(q_off[n] + q_lens.get(n, 0) / 2) / 1e6 for n in q_names])
    ax.set_yticklabels(
        [re.sub(r".*?(Chr0*)", "Chr", n) for n in q_names],
        fontsize=7,
    )
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_linewidth(0.6)

    from matplotlib.lines import Line2D

    ax.legend(
        handles=[
            Line2D([0], [0], color="#1f4e79", lw=2, label="Forward"),
            Line2D([0], [0], color="#b22222", lw=2, label="Reverse"),
        ],
        loc="upper left",
        frameon=False,
        fontsize=9,
    )
    fig.tight_layout()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"[dotplot] hits={len(hits)} fwd={len(segs_fwd)} rev={len(segs_rev)} -> {out}")


def _read_text(path: Path) -> str:
    return path.read_text() if path.is_file() else ""


def cmd_summarize(args: argparse.Namespace) -> None:
    sample = args.sample
    out_tsv = Path(args.out_tsv)
    out_txt = Path(args.out_txt)
    rows = []

    # seqkit stats (optional precomputed)
    stats_file = Path(args.stats)
    if stats_file.is_file():
        with stats_file.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for r in reader:
                rows.append(("assembly", "num_seqs", r.get("num_seqs", "")))
                rows.append(("assembly", "sum_len", r.get("sum_len", "").replace(",", "")))
                rows.append(("assembly", "N50", r.get("N50", "").replace(",", "")))
                rows.append(("assembly", "GC(%)", r.get("GC(%)", "")))
                rows.append(("assembly", "sum_gap", r.get("sum_gap", "")))

    # Merqury QV + k-mer completeness
    # Genome-level *_merqury.qv (NOT per-scaffold *_merqury.*.qv):
    #   Assembly  k-mers_unique  k-mers_total  QV  Error
    # completeness.stats:
    #   Assembly  Region  Found  Total  Completeness(%)
    merq_dir = Path(args.merqury_dir) if args.merqury_dir else None
    if merq_dir and merq_dir.is_dir():
        qv_cands = sorted(merq_dir.glob("*_merqury.qv"))
        qv_path = qv_cands[0] if qv_cands else None
        if qv_path is None:
            for cand in sorted(merq_dir.glob("*.qv")):
                lines = [
                    ln
                    for ln in cand.read_text().splitlines()
                    if ln.strip() and not ln.startswith("#")
                ]
                if len(lines) == 1 and len(lines[0].split()) >= 5:
                    qv_path = cand
                    break
        if qv_path is not None:
            with qv_path.open() as fh:
                for line in fh:
                    if not line.strip() or line.startswith("#"):
                        continue
                    p = line.split()
                    if len(p) >= 5 and p[0] != "Assembly":
                        rows.append(("merqury", "assembly", p[0]))
                        rows.append(("merqury", "kmers_unique", p[1]))
                        rows.append(("merqury", "total_kmers", p[2]))
                        rows.append(("merqury", "QV", p[3]))
                        rows.append(("merqury", "error_rate", p[4]))
                        break

        comp = merq_dir / "completeness.stats"
        if comp.is_file():
            with comp.open() as fh:
                for line in fh:
                    if not line.strip() or line.startswith("#"):
                        continue
                    p = line.split()
                    if len(p) >= 5 and p[1] == "all":
                        rows.append(("merqury", "completeness(%)", p[4]))
                        rows.append(("merqury", "completeness_found", p[2]))
                        rows.append(("merqury", "completeness_total", p[3]))
                        rows.append(("merqury", "completeness_region", p[1]))
                        break

        spectra = sorted(merq_dir.glob("*spectra-cn*.png"))
        if spectra:
            rows.append(
                (
                    "merqury",
                    "spectra_cn_plots",
                    ",".join(p.name for p in spectra[:6]),
                )
            )

    # BUSCO short summary
    if args.busco_dir:
        busco_dir = Path(args.busco_dir)
        shorts = list(busco_dir.rglob("short_summary*.txt"))
        if shorts:
            txt = _read_text(shorts[0])
            m = re.search(
                r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)",
                txt,
            )
            if m:
                rows.append(("busco", "C(%)", m.group(1)))
                rows.append(("busco", "S(%)", m.group(2)))
                rows.append(("busco", "D(%)", m.group(3)))
                rows.append(("busco", "F(%)", m.group(4)))
                rows.append(("busco", "M(%)", m.group(5)))
                rows.append(("busco", "n", m.group(6)))
            m2 = re.search(r"The lineage dataset is:\s+(\S+)", txt)
            if m2:
                rows.append(("busco", "lineage", m2.group(1)))

    # QUAST report.tsv
    if args.quast_dir:
        rep = Path(args.quast_dir) / "report.tsv"
        if rep.is_file():
            with rep.open() as fh:
                reader = csv.reader(fh, delimiter="\t")
                header = next(reader, None)
                want = {
                    "Total length",
                    "N50",
                    "NG50",
                    "L50",
                    "GC (%)",
                    "# misassemblies",
                    "Genome fraction (%)",
                    "Duplication ratio",
                    "# N's per 100 kbp",
                    "Largest alignment",
                    "Total aligned length",
                }
                for line in reader:
                    if len(line) < 2:
                        continue
                    if line[0] in want:
                        rows.append(("quast", line[0], line[1]))

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "section", "metric", "value"])
        for sec, met, val in rows:
            w.writerow([sample, sec, met, val])

    with out_txt.open("w") as fh:
        fh.write(f"========== {sample} genome QC ==========\n")
        cur = None
        for sec, met, val in rows:
            if sec != cur:
                fh.write(f"\n----- {sec} -----\n")
                cur = sec
            fh.write(f"{met}\t{val}\n")
    print(f"[summarize] {out_tsv} ; {out_txt}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    sp = ap.add_subparsers(dest="cmd", required=True)

    p = sp.add_parser("dotplot")
    p.add_argument("--paf", default="")
    p.add_argument("--coords", default="")
    p.add_argument("--query-fai", default="")
    p.add_argument("--ref-fai", default="")
    p.add_argument("--out", required=True)
    p.add_argument("--title", default="Synteny dotplot")
    p.add_argument("--query-label", default="Query")
    p.add_argument("--ref-label", default="Reference")
    p.add_argument("--min-alen", type=int, default=20000)
    p.add_argument("--min-mapq", type=int, default=20)
    p.set_defaults(func=cmd_dotplot)

    p = sp.add_parser("summarize")
    p.add_argument("--sample", required=True)
    p.add_argument("--stats", default="")
    p.add_argument("--merqury-dir", default="")
    p.add_argument("--busco-dir", default="")
    p.add_argument("--quast-dir", default="")
    p.add_argument("--out-tsv", required=True)
    p.add_argument("--out-txt", required=True)
    p.set_defaults(func=cmd_summarize)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
