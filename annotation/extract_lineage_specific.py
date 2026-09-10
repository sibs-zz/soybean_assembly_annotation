#!/usr/bin/env python3
"""Extract cultivar-specific (ZH13-divergent) annotation files from 12_renamed.

Filter: Chr01–20 + confidence >= medium + TE overlap <= 0.5
        + expression_supported=yes + ZH13 v3 protein identity < 50
        (no ZH13 hit counts as identity 0).
Disjoint from the HC set (HC requires ZH13 pident >= 50).
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from build_publication_annotation import (  # noqa: E402
    annotate_gff,
    confidence_pass,
    extract_proteins,
    is_nuclear_chrom,
    load_summary,
    parse_attr,
    zh13_pident,
)


def lineage_pass(row: dict[str, str], min_conf: str, max_te: float, max_zh13: float) -> bool:
    if not is_nuclear_chrom(row.get("chrom", "")):
        return False
    if not confidence_pass(row, min_conf, max_te, require_expr=True):
        return False
    return zh13_pident(row) < max_zh13


def tx2gene_from_gff(gff: str) -> dict[str, str]:
    out: dict[str, str] = {}
    with open(gff) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] not in ("mRNA", "transcript"):
                continue
            a = parse_attr(p[8])
            tid, parent = a.get("ID"), a.get("Parent", "")
            if tid and parent:
                out[tid] = parent.split(",")[0]
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dir", required=True, help="12_renamed directory")
    ap.add_argument("--prefix", required=True, help="sample prefix, e.g. KD35")
    ap.add_argument("--min-confidence", default="medium_confidence")
    ap.add_argument("--max-te", type=float, default=0.5)
    ap.add_argument("--max-zh13", type=float, default=50.0,
                    help="keep genes with ZH13 v3 pident strictly below this")
    args = ap.parse_args()

    d = args.dir
    p = args.prefix
    summary = os.path.join(d, f"{p}.renamed.annotation.summary.tsv")
    gff = os.path.join(d, f"{p}.renamed.gff3")
    tx_fa = os.path.join(d, f"{p}.renamed.transcripts.fa")
    cds_fa = os.path.join(d, f"{p}.renamed.cds.fa")
    prot_fa = os.path.join(d, f"{p}.renamed.proteins.fa")
    for path in (summary, gff, tx_fa, cds_fa, prot_fa):
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            print(f"[ERROR] missing or empty: {path}", file=sys.stderr)
            return 1

    stem = os.path.join(d, f"{p}.renamed.lineage_specific")
    out_list = f"{stem}.genes.list"
    out_gff = f"{stem}.annotated.gff3"
    out_tx = f"{stem}.transcripts.fa"
    out_cds = f"{stem}.cds.fa"
    out_prot = f"{stem}.proteins.fa"
    out_readme = f"{stem}.README.txt"

    rows = load_summary(summary)
    genes = {
        gid
        for gid, r in rows.items()
        if lineage_pass(r, args.min_confidence, args.max_te, args.max_zh13)
    }
    ordered = sorted(genes, key=lambda x: (rows[x].get("chrom", ""), int(rows[x].get("start", 0) or 0)))
    with open(out_list, "w") as f:
        for gid in ordered:
            f.write(f"{gid}\n")

    n_gff = annotate_gff(gff, out_gff, rows, gene_filter=genes)
    tx2gene = tx2gene_from_gff(gff)
    n_tx = extract_proteins(tx_fa, genes, tx2gene, out_tx)
    n_cds = extract_proteins(cds_fa, genes, tx2gene, out_cds)
    n_prot = extract_proteins(prot_fa, genes, tx2gene, out_prot)

    n_high = sum(1 for g in genes if rows[g].get("confidence_level") == "high_confidence")
    n_nohit = sum(1 for g in genes if zh13_pident(rows[g]) == 0.0)
    n_mid = len(genes) - n_nohit

    with open(out_readme, "w") as f:
        f.write(
            f"""Lineage-specific (ZH13-divergent) annotation for {p}
================================================================

Not a replacement for the HC set. These loci are kept in the full nuclear GFF;
this subset is for cultivar-specific / distant-homolog analyses and supplements.

Filter (all required):
  Chr01–20
  confidence >= {args.min_confidence}
  TE overlap ratio <= {args.max_te}
  expression_supported = yes
  ZH13 v3 protein identity < {args.max_zh13}  (no hit = 0)

Counts:
  genes              {len(genes)}
  GFF gene features  {n_gff}
  transcripts/mRNA   {n_tx}
  CDS                {n_cds}
  proteins           {n_prot}
  high_confidence    {n_high}
  no ZH13 v3 hit     {n_nohit}
  ZH13 hit but <{args.max_zh13:.0f}  {n_mid}

Files:
  {os.path.basename(out_list)}
  {os.path.basename(out_gff)}
  {os.path.basename(out_tx)}
  {os.path.basename(out_cds)}
  {os.path.basename(out_prot)}

Source: {summary}
"""
        )

    print(f"[INFO] lineage-specific genes: {out_list} ({len(genes)})")
    print(f"[INFO] lineage-specific GFF3: {out_gff} ({n_gff} genes)")
    print(f"[INFO] lineage-specific transcripts: {out_tx} ({n_tx})")
    print(f"[INFO] lineage-specific CDS: {out_cds} ({n_cds})")
    print(f"[INFO] lineage-specific proteins: {out_prot} ({n_prot})")
    print(f"[INFO] README: {out_readme}")
    if n_gff != len(genes):
        print(f"[WARN] GFF gene count {n_gff} != list {len(genes)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
