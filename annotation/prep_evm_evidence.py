#!/usr/bin/env python3
"""Normalize gene / transcript / protein GFF into EVM 2.x input tracks."""

from __future__ import annotations

import argparse
import sys


def parse_attr(s: str) -> dict[str, str]:
    d: dict[str, str] = {}
    for item in s.strip().split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
        elif " " in item:
            k, v = item.split(" ", 1)
            d[k] = v.strip().strip('"')
    return d


def fmt_attr(d: dict[str, str]) -> str:
    return ";".join(f"{k}={v}" for k, v in d.items() if v)


def rewrite_source(path: str, source: str, keep_types: set[str], out) -> int:
    n = 0
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            if p[2] not in keep_types:
                continue
            p[1] = source
            if p[2] == "transcript":
                p[2] = "mRNA"
            out.write("\t".join(p) + "\n")
            n += 1
    return n


def miniprot_to_protein(path: str, source: str, out) -> int:
    n = 0
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "CDS":
                continue
            a = parse_attr(p[8])
            parent = a.get("Parent", a.get("ID", f"mp{n}"))
            target = a.get("Target", parent)
            p[1] = source
            p[2] = "nucleotide_to_protein_match"
            p[8] = fmt_attr({"ID": a.get("ID", f"{parent}.cds{n}"), "Target": target})
            out.write("\t".join(p) + "\n")
            n += 1
    return n


def pasa_to_transcript(path: str, source: str, out) -> int:
    n = 0
    keep = {
        "cDNA_match",
        "EST_match",
        "expressed_sequence_match",
        "nucleotide_to_EST_match",
        "match",
        "match_part",
        "exon",
    }
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            if p[2] not in keep:
                continue
            a = parse_attr(p[8])
            tid = a.get("ID") or a.get("Parent") or a.get("Target") or f"pasa{n}"
            target = a.get("Target", tid)
            p[1] = source
            p[2] = "nucleotide_to_EST_match"
            p[8] = fmt_attr({"ID": a.get("ID", tid), "Target": target})
            out.write("\t".join(p) + "\n")
            n += 1
    return n


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", action="append", nargs=2, metavar=("SOURCE", "FILE"), default=[])
    ap.add_argument("--protein", nargs=2, metavar=("SOURCE", "FILE"))
    ap.add_argument("--transcript", nargs=2, metavar=("SOURCE", "FILE"))
    ap.add_argument("--gene-out", default="")
    ap.add_argument("--protein-out", default="")
    ap.add_argument("--transcript-out", default="")
    args = ap.parse_args()

    gene_types = {"gene", "mRNA", "transcript", "exon", "CDS"}
    if args.gene_out:
        n = 0
        with open(args.gene_out, "w") as out:
            out.write("##gff-version 3\n")
            for source, path in args.gene:
                n += rewrite_source(path, source, gene_types, out)
        print(f"[INFO] gene_predictions lines={n} -> {args.gene_out}", file=sys.stderr)
    if args.protein and args.protein_out:
        source, path = args.protein
        with open(args.protein_out, "w") as out:
            out.write("##gff-version 3\n")
            n = miniprot_to_protein(path, source, out)
        print(f"[INFO] protein_alignments lines={n} -> {args.protein_out}", file=sys.stderr)
    if args.transcript and args.transcript_out:
        source, path = args.transcript
        with open(args.transcript_out, "w") as out:
            out.write("##gff-version 3\n")
            n = pasa_to_transcript(path, source, out)
        print(f"[INFO] transcript_alignments lines={n} -> {args.transcript_out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
