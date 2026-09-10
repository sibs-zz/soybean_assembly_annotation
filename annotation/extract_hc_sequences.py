#!/usr/bin/env python3
"""Subset HC mRNA (transcripts) and CDS FASTA from 12_renamed using the HC gene list."""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from build_publication_annotation import extract_proteins, parse_attr  # noqa: E402


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


def load_list(path: str) -> set[str]:
    genes: set[str] = set()
    with open(path) as f:
        for line in f:
            gid = line.strip()
            if gid:
                genes.add(gid)
    return genes


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dir", required=True, help="12_renamed directory")
    ap.add_argument("--prefix", required=True, help="sample prefix, e.g. KD35")
    args = ap.parse_args()

    d = args.dir
    p = args.prefix
    gene_list = os.path.join(d, f"{p}.renamed.high_confidence.genes.list")
    gff = os.path.join(d, f"{p}.renamed.gff3")
    tx_fa = os.path.join(d, f"{p}.renamed.transcripts.fa")
    cds_fa = os.path.join(d, f"{p}.renamed.cds.fa")
    out_tx = os.path.join(d, f"{p}.renamed.high_confidence.transcripts.fa")
    out_cds = os.path.join(d, f"{p}.renamed.high_confidence.cds.fa")

    for path in (gene_list, gff, tx_fa, cds_fa):
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            print(f"[ERROR] missing or empty: {path}", file=sys.stderr)
            return 1

    genes = load_list(gene_list)
    tx2gene = tx2gene_from_gff(gff)
    n_tx = extract_proteins(tx_fa, genes, tx2gene, out_tx)
    n_cds = extract_proteins(cds_fa, genes, tx2gene, out_cds)
    print(f"[INFO] HC transcripts/mRNA: {out_tx} ({n_tx})")
    print(f"[INFO] HC CDS: {out_cds} ({n_cds})")
    print(f"[INFO] HC genes in list: {len(genes)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
