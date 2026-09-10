#!/usr/bin/env python3
"""Rename gene loci to publication IDs (e.g. KD35_01G000100 / .m1) by chrom + coordinate.

Gene IDs follow ZH13 v3 (PREFIX_chrG + 6 digits, step 100). Transcripts from EVM
(evm.model.*), PASA (novel_model_*/temp_model_*), or BRAKER (g1.t1) become
PREFIX_01G000100.m1; exons/CDS/UTRs become .m1.exon1 / .cds1 / .utr5p1 / .utr3p1.
"""

from __future__ import annotations

import csv
import os
import re
import sys
from collections import defaultdict


def parse_attr(s: str) -> dict[str, str]:
    d: dict[str, str] = {}
    for item in s.strip().split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d


def fmt_attr(d: dict[str, str]) -> str:
    return ";".join(f"{k}={v}" for k, v in d.items() if v)


def chrom_sort_key(chrom: str) -> tuple:
    m = re.match(r"^Chr(\d+)$", chrom, re.I)
    if m:
        return (0, int(m.group(1)), chrom)
    m = re.match(r"^Unanchored_(\d+)$", chrom, re.I)
    if m:
        return (1, int(m.group(1)), chrom)
    return (2, chrom, chrom)


def chrom_code(chrom: str) -> str:
    m = re.match(r"^Chr(\d+)$", chrom, re.I)
    if m:
        return f"{int(m.group(1)):02d}"
    m = re.match(r"^Unanchored_(\d+)$", chrom, re.I)
    if m:
        return f"U{int(m.group(1)):02d}"
    safe = re.sub(r"[^A-Za-z0-9]", "", chrom)[:6]
    return safe or "XX"


def merge_alias(attrs: dict[str, str], old_gene_id: str) -> None:
    alias_val = old_gene_id
    existing = attrs.get("Alias", "")
    if existing and old_gene_id not in existing.split(","):
        attrs["Alias"] = f"{existing},{alias_val}"
    else:
        attrs["Alias"] = alias_val
    note = attrs.get("Note", "")
    braker_note = f"BRAKER_ID={old_gene_id}"
    if note and braker_note not in note:
        attrs["Note"] = f"{note}; {braker_note}"
    elif not note:
        attrs["Note"] = braker_note


def build_gene_rename_map(
    gff_path: str,
    id_prefix: str,
    gene_digits: int,
    gene_step: int = 100,
) -> tuple[dict[str, str], dict[str, dict[str, str]]]:
    """Return old_gene_id -> new_gene_id and metadata per old gene."""
    genes: list[tuple[str, int, int, str, str]] = []
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "gene":
                continue
            a = parse_attr(p[8])
            gid = a.get("ID")
            if not gid:
                continue
            genes.append((p[0], int(p[3]), int(p[4]), p[6], gid))

    genes.sort(key=lambda x: (chrom_sort_key(x[0]), x[1], x[4]))

    gene_map: dict[str, str] = {}
    meta: dict[str, dict[str, str]] = {}
    counters: dict[str, int] = defaultdict(int)

    step = max(1, int(gene_step))
    for chrom, start, end, strand, old_gid in genes:
        cc = chrom_code(chrom)
        counters[cc] += 1
        rank = counters[cc]
        locus_num = rank * step
        new_gid = f"{id_prefix}_{cc}G{locus_num:0{gene_digits}d}"
        gene_map[old_gid] = new_gid
        meta[old_gid] = {
            "old_gene_id": old_gid,
            "new_gene_id": new_gid,
            "chrom": chrom,
            "chrom_code": cc,
            "start": str(start),
            "end": str(end),
            "strand": strand,
            "rank_on_chrom": str(rank),
            "locus_number": str(locus_num),
            "gene_step": str(step),
        }
    return gene_map, meta


_BRAKER_ID_RE = re.compile(r"^(g\d+)(.*)$")
_TX_TYPES = {"mRNA", "transcript"}
_CHILD_SUFFIX = {
    "exon": "exon",
    "CDS": "cds",
    "five_prime_UTR": "utr5p",
    "five_prime_utr": "utr5p",
    "three_prime_UTR": "utr3p",
    "three_prime_utr": "utr3p",
}


def build_feature_id_map(gff_path: str, gene_map: dict[str, str]) -> dict[str, str]:
    """Map gene/mRNA/exon/CDS/UTR IDs onto ZH13-style PREFIX_chrG######.mN.* names."""
    id_map = dict(gene_map)
    mrnas_by_gene: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    strand_by_mrna: dict[str, str] = {}
    kids_by_mrna: dict[str, list[tuple[str, int, int, str]]] = defaultdict(list)

    with open(gff_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            a = parse_attr(p[8])
            fid = a.get("ID")
            parent = (a.get("Parent") or "").split(",")[0].strip()
            ftype = p[2]
            start, end = int(p[3]), int(p[4])
            if ftype in _TX_TYPES and fid and parent in gene_map:
                mrnas_by_gene[parent].append((start, end, fid))
                strand_by_mrna[fid] = p[6]
            elif ftype in _CHILD_SUFFIX and fid and parent:
                kids_by_mrna[parent].append((ftype, start, end, fid))

    for gid, mrnas in mrnas_by_gene.items():
        new_gene = gene_map[gid]
        mrnas.sort(key=lambda x: (x[0], x[1], x[2]))
        for i, (_s, _e, old_tid) in enumerate(mrnas, start=1):
            new_tid = f"{new_gene}.m{i}"
            id_map[old_tid] = new_tid
            reverse = strand_by_mrna.get(old_tid, "+") == "-"
            counters: dict[str, int] = defaultdict(int)
            kids = kids_by_mrna.get(old_tid, [])
            kids.sort(key=lambda x: (-x[1], -x[2], x[0]) if reverse else (x[1], x[2], x[0]))
            for ftype, _cs, _ce, old_cid in kids:
                tag = _CHILD_SUFFIX[ftype]
                counters[tag] += 1
                id_map[old_cid] = f"{new_tid}.{tag}{counters[tag]}"
    return id_map


def remap_feature_id(
    old_id: str,
    id_map: dict[str, str],
    gene_map: dict[str, str] | None = None,
) -> str:
    """Map EVM/PASA/BRAKER feature IDs onto PREFIX_chrG / .m1 names."""
    if old_id in id_map:
        return id_map[old_id]
    gene_map = gene_map or {}
    m = _BRAKER_ID_RE.match(old_id)
    if m:
        base, suffix = m.group(1), m.group(2)
        if base in gene_map:
            suffix = re.sub(r"\.t(\d+)", r".m\1", suffix)
            return gene_map[base] + suffix
    if old_id in gene_map:
        return gene_map[old_id]
    return old_id


def remap_parent_field(parent: str, id_map: dict[str, str], gene_map: dict[str, str]) -> str:
    parts = []
    for pid in parent.split(","):
        pid = pid.strip()
        if pid:
            parts.append(remap_feature_id(pid, id_map, gene_map))
    return ",".join(parts)


def rename_gff(
    in_gff: str,
    out_gff: str,
    gene_map: dict[str, str],
    id_map: dict[str, str],
) -> int:
    n_genes = 0
    rev_gene = {v: k for k, v in gene_map.items()}
    with open(in_gff) as fin, open(out_gff, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            if not line.strip():
                fout.write(line)
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                fout.write(line + "\n")
                continue

            a = parse_attr(p[8])
            if "ID" in a:
                a["ID"] = remap_feature_id(a["ID"], id_map, gene_map)
            if "Parent" in a:
                a["Parent"] = remap_parent_field(a["Parent"], id_map, gene_map)
            for key in ("Name",):
                if key in a:
                    a[key] = remap_feature_id(a[key], id_map, gene_map)

            if p[2] == "gene":
                old_gene = rev_gene.get(a.get("ID", ""))
                if old_gene:
                    merge_alias(a, old_gene)
                    n_genes += 1

            p[8] = fmt_attr(a)
            fout.write("\t".join(p) + "\n")
    return n_genes


def rename_fasta(
    in_fa: str,
    out_fa: str,
    id_map: dict[str, str],
    gene_map: dict[str, str],
) -> int:
    n = 0
    with open(in_fa) as fin, open(out_fa, "w") as fout:
        for line in fin:
            if not line.startswith(">"):
                fout.write(line)
                continue
            head = line[1:].rstrip("\n")
            sid = head.split()[0]
            new_sid = remap_feature_id(sid, id_map, gene_map)
            rest = head[len(sid) :]
            fout.write(f">{new_sid}{rest}\n")
            n += 1
    return n


def rename_tsv_gene_column(
    in_tsv: str,
    out_tsv: str,
    gene_map: dict[str, str],
    gene_col: str = "gene_id",
) -> int:
    n = 0
    with open(in_tsv) as fin, open(out_tsv, "w", newline="") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        if not reader.fieldnames:
            return 0
        fields = list(reader.fieldnames)
        w = csv.DictWriter(fout, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        rows = list(reader)
        if gene_col in fields and rows and "chrom" in fields and "start" in fields:
            rows.sort(
                key=lambda r: (
                    chrom_sort_key(r.get("chrom", "")),
                    int(r.get("start", 0) or 0),
                )
            )
        for row in rows:
            old = row.get(gene_col, "")
            if old in gene_map:
                row[gene_col] = gene_map[old]
                n += 1
            w.writerow(row)
    return n


def write_mapping_table(
    path: str,
    gene_map: dict[str, str],
    meta: dict[str, dict[str, str]],
) -> None:
    header = [
        "old_gene_id",
        "new_gene_id",
        "chrom",
        "chrom_code",
        "start",
        "end",
        "strand",
        "rank_on_chrom",
        "locus_number",
        "gene_step",
    ]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header, delimiter="\t")
        w.writeheader()
        for old_gid in sorted(
            gene_map,
            key=lambda g: (
                chrom_sort_key(meta[g]["chrom"]),
                int(meta[g]["start"]),
            ),
        ):
            w.writerow(meta[old_gid])


def rename_gene_list(in_list: str, out_list: str, gene_map: dict[str, str], meta: dict[str, dict[str, str]]) -> int:
    ids = []
    with open(in_list) as f:
        for line in f:
            gid = line.strip()
            if gid and gid in gene_map:
                ids.append(gene_map[gid])
    ids.sort(
        key=lambda g: next(
            (
                (
                    chrom_sort_key(meta[og]["chrom"]),
                    int(meta[og]["start"]),
                )
                for og, ng in gene_map.items()
                if ng == g
            ),
            ((), 0),
        )
    )
    with open(out_list, "w") as f:
        for gid in ids:
            f.write(f"{gid}\n")
    return len(ids)


def main() -> int:
    in_gff = os.environ["RENAME_IN_GFF"]
    out_gff = os.environ["RENAME_OUT_GFF"]
    id_prefix = os.environ.get("RENAME_ID_PREFIX", os.environ.get("RENAME_PREFIX", "KD35"))
    gene_digits = int(os.environ.get("RENAME_GENE_DIGITS", "6"))
    gene_step = int(os.environ.get("RENAME_GENE_ID_STEP", "100"))
    mapping_out = os.environ["RENAME_MAPPING_TSV"]

    gene_map, meta = build_gene_rename_map(in_gff, id_prefix, gene_digits, gene_step)
    if not gene_map:
        print("[ERROR] No gene features found in GFF3", file=sys.stderr)
        return 1

    id_map = build_feature_id_map(in_gff, gene_map)
    n_tx = sum(1 for v in id_map.values() if re.search(r"\.m\d+$", v))
    print(
        f"[INFO] Assigned {len(gene_map)} gene loci "
        f"({id_prefix}_{{chrom}}G, step={gene_step}, e.g. ...G{gene_step:0{gene_digits}d})"
    )
    print(f"[INFO] Mapped {n_tx} transcripts to {{gene}}.mN (EVM/PASA/BRAKER)")
    write_mapping_table(mapping_out, gene_map, meta)

    n = rename_gff(in_gff, out_gff, gene_map, id_map)
    print(f"[INFO] Renamed GFF3: {out_gff} ({n} genes)")

    pairs = [
        ("RENAME_IN_GENE_TABLE", "RENAME_OUT_GENE_TABLE"),
        ("RENAME_IN_SUMMARY", "RENAME_OUT_SUMMARY"),
    ]
    for env_in, env_out in pairs:
        if env_in in os.environ and env_out in os.environ and os.path.isfile(os.environ[env_in]):
            n_rows = rename_tsv_gene_column(os.environ[env_in], os.environ[env_out], gene_map)
            print(f"[INFO] Renamed table: {os.environ[env_out]} ({n_rows} rows)")

    for env_in, env_out in (
        ("RENAME_IN_TRANSCRIPT_FA", "RENAME_OUT_TRANSCRIPT_FA"),
        ("RENAME_IN_CDS_FA", "RENAME_OUT_CDS_FA"),
        ("RENAME_IN_PROTEIN_FA", "RENAME_OUT_PROTEIN_FA"),
    ):
        if env_in in os.environ and env_out in os.environ and os.path.isfile(os.environ[env_in]):
            n_seq = rename_fasta(os.environ[env_in], os.environ[env_out], id_map, gene_map)
            print(f"[INFO] Renamed FASTA: {os.environ[env_out]} ({n_seq} sequences)")

    if all(k in os.environ for k in ("RENAME_IN_HC_GFF", "RENAME_OUT_HC_GFF")):
        if os.path.isfile(os.environ["RENAME_IN_HC_GFF"]):
            n_hc = rename_gff(os.environ["RENAME_IN_HC_GFF"], os.environ["RENAME_OUT_HC_GFF"], gene_map, id_map)
            print(f"[INFO] Renamed HC GFF3: {os.environ['RENAME_OUT_HC_GFF']} ({n_hc} genes)")

    if all(k in os.environ for k in ("RENAME_IN_HC_PROTEIN", "RENAME_OUT_HC_PROTEIN")):
        if os.path.isfile(os.environ["RENAME_IN_HC_PROTEIN"]):
            n_hp = rename_fasta(
                os.environ["RENAME_IN_HC_PROTEIN"],
                os.environ["RENAME_OUT_HC_PROTEIN"],
                id_map,
                gene_map,
            )
            print(f"[INFO] Renamed HC proteins: {os.environ['RENAME_OUT_HC_PROTEIN']} ({n_hp} seqs)")

    if all(k in os.environ for k in ("RENAME_IN_HC_LIST", "RENAME_OUT_HC_LIST")):
        if os.path.isfile(os.environ["RENAME_IN_HC_LIST"]):
            n_hl = rename_gene_list(
                os.environ["RENAME_IN_HC_LIST"],
                os.environ["RENAME_OUT_HC_LIST"],
                gene_map,
                meta,
            )
            print(f"[INFO] Renamed HC gene list: {os.environ['RENAME_OUT_HC_LIST']} ({n_hl} genes)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
