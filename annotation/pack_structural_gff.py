#!/usr/bin/env python3
"""Normalize annotation to a ZH13-like structural GFF3.

Keeps: gene / mRNA / exon / CDS / five_prime_UTR / three_prime_UTR
Drops: intron, start_codon, stop_codon, and extra attributes.
Optional: nuclear-only (Chr01–Chr20) and minimum protein length.
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field

KEEP = {
    "gene",
    "mRNA",
    "transcript",
    "exon",
    "CDS",
    "five_prime_UTR",
    "five_prime_utr",
    "three_prime_UTR",
    "three_prime_utr",
}
UTR5 = {"five_prime_UTR", "five_prime_utr", "5'-UTR", "five_prime_utr"}
UTR3 = {"three_prime_UTR", "three_prime_utr", "3'-UTR", "three_prime_utr"}
DROP = {"intron", "start_codon", "stop_codon"}
GFF_ATTR = re.compile(r"(\w+)=([^;]*)")
GTF_ATTR = re.compile(r'(\w+)\s+"([^"]*)"')
TX_SUFFIX = re.compile(r"\.[tm]\d+")
TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")
AA = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


@dataclass
class Feature:
    chrom: str
    source: str
    ftype: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attrs: dict = field(default_factory=dict)


def parse_attr_blob(s: str) -> dict[str, str]:
    if "=" in s and "transcript_id \"" not in s[:20]:
        d = {k: v for k, v in GFF_ATTR.findall(s)}
        if d:
            return d
    d = {k: v for k, v in GTF_ATTR.findall(s)}
    if d:
        return d
    s = s.strip().rstrip(";")
    if s and " " not in s and "=" not in s:
        return {"ID": s}
    return {}


def gene_id_of(tid: str) -> str:
    return TX_SUFFIX.split(tid)[0]


def load_fasta(path: str) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    name = None
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def revcomp(seq: str) -> str:
    return seq.translate(TABLE)[::-1]


def translate_cds(
    genome: dict[str, str], chrom: str, strand: str, cds: list[tuple[int, int]]
) -> str | None:
    if chrom not in genome or not cds:
        return None
    seq = genome[chrom]
    chunks = []
    ordered = sorted(cds, reverse=(strand == "-"))
    for s, e in ordered:
        chunk = seq[s - 1 : e]
        if strand == "-":
            chunk = revcomp(chunk)
        chunks.append(chunk)
    nt = "".join(chunks).upper()
    nt = re.sub(r"[^ACGT]", "N", nt)
    cut = len(nt) - (len(nt) % 3)
    if cut < 3:
        return None
    aa = []
    for i in range(0, cut, 3):
        aa.append(AA.get(nt[i : i + 3], "X"))
    prot = "".join(aa)
    if "*" in prot[:-1]:
        return None
    if prot.endswith("*"):
        prot = prot[:-1]
    return prot or None


def parse_gff_or_gtf(path: str) -> list[Feature]:
    feats: list[Feature] = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            ft = p[2]
            if ft in DROP:
                continue
            if ft not in KEEP:
                continue
            attrs = parse_attr_blob(p[8])
            if ft == "gene":
                gid = attrs.get("ID") or attrs.get("gene_id") or p[8].strip().rstrip(";")
                attrs = {"ID": gid}
            elif ft in ("mRNA", "transcript"):
                tid = attrs.get("ID") or attrs.get("transcript_id") or p[8].strip().rstrip(";")
                gid = attrs.get("Parent") or attrs.get("gene_id") or gene_id_of(tid)
                attrs = {"ID": tid, "Parent": gid}
                ft = "mRNA"
            else:
                tid = attrs.get("Parent") or attrs.get("transcript_id") or attrs.get("ID", "")
                if ft in UTR5:
                    ft = "five_prime_UTR"
                elif ft in UTR3:
                    ft = "three_prime_UTR"
                attrs = {"Parent": tid}
            feats.append(
                Feature(p[0], p[1], ft, int(p[3]), int(p[4]), p[5], p[6], p[7], attrs)
            )
    return feats


def is_nuclear(chrom: str) -> bool:
    return bool(re.match(r"^Chr\d+$", chrom, re.I))


def rebuild(feats: list[Feature]) -> list[Feature]:
    genes: dict[str, Feature] = {}
    mrnas: dict[str, Feature] = {}
    kids: dict[str, list[Feature]] = defaultdict(list)

    for ft in feats:
        if ft.ftype == "gene":
            genes[ft.attrs["ID"]] = ft
        elif ft.ftype == "mRNA":
            mrnas[ft.attrs["ID"]] = ft
        else:
            kids[ft.attrs.get("Parent", "")].append(ft)

    out: list[Feature] = []
    for tid, mf in mrnas.items():
        gid = mf.attrs.get("Parent", gene_id_of(tid))
        child = kids.get(tid, [])
        if not child:
            continue
        for j, c in enumerate(child, start=1):
            if "ID" not in c.attrs or not c.attrs["ID"]:
                suffix = {"exon": "exon", "CDS": "cds", "five_prime_UTR": "utr5", "three_prime_UTR": "utr3"}.get(
                    c.ftype, c.ftype
                )
                c.attrs["ID"] = f"{tid}.{suffix}{j}"
            c.attrs["Parent"] = tid
        mf.start = min(c.start for c in child)
        mf.end = max(c.end for c in child)
        if gid not in genes:
            genes[gid] = Feature(mf.chrom, mf.source, "gene", mf.start, mf.end, ".", mf.strand, ".", {"ID": gid})
        else:
            g = genes[gid]
            g.start = min(g.start, mf.start)
            g.end = max(g.end, mf.end)

    for gid, g in genes.items():
        g.ftype = "gene"
        g.attrs = {"ID": gid}
        out.append(g)
        for tid, mf in mrnas.items():
            if mf.attrs.get("Parent") != gid:
                continue
            mf.ftype = "mRNA"
            mf.attrs = {"ID": tid, "Parent": gid}
            out.append(mf)
            for c in sorted(kids.get(tid, []), key=lambda x: (x.start, x.end, x.ftype)):
                out.append(c)
    return out


def write_gff3(path: str, feats: list[Feature]) -> None:
    with open(path, "w") as f:
        f.write("##gff-version 3\n")
        for ft in feats:
            attr = ";".join(f"{k}={v}" for k, v in ft.attrs.items() if v)
            f.write(
                "\t".join(
                    [
                        ft.chrom,
                        ft.source,
                        ft.ftype,
                        str(ft.start),
                        str(ft.end),
                        ft.score or ".",
                        ft.strand or ".",
                        ft.phase or ".",
                        attr,
                    ]
                )
                + "\n"
            )


def filter_nuclear(feats: list[Feature]) -> list[Feature]:
    return [ft for ft in feats if is_nuclear(ft.chrom)]


def filter_min_aa(feats: list[Feature], genome: dict[str, str], min_aa: int) -> list[Feature]:
    cds_by_tx: dict[str, list[tuple[int, int]]] = defaultdict(list)
    tx_meta: dict[str, Feature] = {}
    for ft in feats:
        if ft.ftype == "mRNA":
            tx_meta[ft.attrs["ID"]] = ft
        elif ft.ftype == "CDS":
            cds_by_tx[ft.attrs.get("Parent", "")].append((ft.start, ft.end))
    keep_tx = set()
    for tid, cds in cds_by_tx.items():
        mf = tx_meta.get(tid)
        if not mf:
            continue
        prot = translate_cds(genome, mf.chrom, mf.strand, cds)
        if prot and len(prot) >= min_aa:
            keep_tx.add(tid)
    keep_gene = {tx_meta[t].attrs.get("Parent") for t in keep_tx if t in tx_meta}
    out = []
    for ft in feats:
        if ft.ftype == "gene" and ft.attrs.get("ID") in keep_gene:
            out.append(ft)
        elif ft.ftype == "mRNA" and ft.attrs.get("ID") in keep_tx:
            out.append(ft)
        elif ft.ftype != "gene" and ft.ftype != "mRNA" and ft.attrs.get("Parent") in keep_tx:
            out.append(ft)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--nuclear-out", default="")
    ap.add_argument("--min-aa", type=int, default=0)
    ap.add_argument("--genome", default="")
    args = ap.parse_args()

    feats = rebuild(parse_gff_or_gtf(args.input))
    write_gff3(args.output, feats)
    print(f"[INFO] Structural GFF3: {args.output}", file=sys.stderr)
    if args.nuclear_out:
        nuc = filter_nuclear(feats)
        if args.min_aa > 0:
            if not args.genome:
                print("[ERROR] --min-aa requires --genome", file=sys.stderr)
                return 1
            nuc = filter_min_aa(nuc, load_fasta(args.genome), args.min_aa)
        nuc = rebuild(nuc)
        write_gff3(args.nuclear_out, nuc)
        print(f"[INFO] Nuclear GFF3: {args.nuclear_out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
