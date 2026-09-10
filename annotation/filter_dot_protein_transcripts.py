#!/usr/bin/env python3
"""Drop transcripts whose gffread protein sequences contain '.' (invalid for DIAMOND)."""

from __future__ import annotations

import os
import shutil
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


def bad_protein_ids(protein_fa: str) -> set[str]:
    bad: set[str] = set()
    sid: str | None = None
    chunks: list[str] = []
    with open(protein_fa) as f:
        for line in f:
            if line.startswith(">"):
                if sid is not None and "." in "".join(chunks):
                    bad.add(sid)
                sid = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if sid is not None and "." in "".join(chunks):
            bad.add(sid)
    return bad


def filter_fasta(in_fa: str, out_fa: str, drop: set[str]) -> tuple[int, int]:
    kept = dropped = 0
    with open(in_fa) as fin, open(out_fa, "w") as fout:
        write = False
        for line in fin:
            if line.startswith(">"):
                sid = line[1:].split()[0]
                write = sid not in drop
                if write:
                    kept += 1
                    fout.write(line)
                else:
                    dropped += 1
            elif write:
                fout.write(line)
    return kept, dropped


def filter_gff(in_gff: str, out_gff: str, drop_tx: set[str]) -> tuple[int, int]:
    gene2tx: dict[str, set[str]] = defaultdict(set)
    tx_lines: dict[str, list[str]] = {}

    with open(in_gff) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            feature = p[2]
            a = parse_attr(p[8])
            if feature == "gene":
                gid = a.get("ID")
                if gid:
                    gene2tx.setdefault(gid, set())
            elif feature in ("mRNA", "transcript"):
                tid = a.get("ID")
                parent = a.get("Parent", "").split(",")[0]
                if tid:
                    if parent:
                        gene2tx[parent].add(tid)
                    else:
                        gene2tx.setdefault(tid, set())

    drop_genes = {g for g, txs in gene2tx.items() if txs and txs.issubset(drop_tx)}

    n_drop_lines = 0
    n_keep = 0
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
                fout.write(line)
                continue
            a = parse_attr(p[8])
            fid = a.get("ID", "")
            parents = [x.strip() for x in a.get("Parent", "").split(",") if x.strip()]

            skip = False
            if p[2] == "gene" and fid in drop_genes:
                skip = True
            elif p[2] in ("mRNA", "transcript") and fid in drop_tx:
                skip = True
            else:
                for pid in parents:
                    if pid in drop_tx or pid in drop_genes:
                        skip = True
                        break
                    if pid in gene2tx and gene2tx[pid] <= drop_tx:
                        skip = True
                        break

            if skip:
                n_drop_lines += 1
            else:
                n_keep += 1
                fout.write(line)
    return n_keep, n_drop_lines


def main() -> int:
    protein_fa = os.environ["QC_PROTEIN_FA"]
    transcript_fa = os.environ["QC_TRANSCRIPT_FA"]
    cds_fa = os.environ["QC_CDS_FA"]
    in_gff = os.environ["QC_IN_GFF"]
    out_gff = os.environ["QC_OUT_GFF"]
    dropped_list = os.environ["QC_DROPPED_LIST"]
    marker = os.environ["QC_MARKER"]

    drop = bad_protein_ids(protein_fa)
    if not drop:
        print("[INFO] No proteins with '.' found; passing files through unchanged.")
        shutil.copyfile(in_gff, out_gff)
        with open(dropped_list, "w") as f:
            f.write("transcript_id\treason\n")
        open(marker, "w").write("ok\n")
        return 0

    print(f"[INFO] Dropping {len(drop)} transcript(s) with '.' in protein sequence:")
    for tid in sorted(drop):
        print(f"       {tid}")

    with open(dropped_list, "w") as f:
        f.write("transcript_id\treason\n")
        for tid in sorted(drop):
            f.write(f"{tid}\tgffread_protein_contains_dot\n")

    tmp_prot = protein_fa + ".tmp"
    tmp_tx = transcript_fa + ".tmp"
    tmp_cds = cds_fa + ".tmp"
    kp, dp = filter_fasta(protein_fa, tmp_prot, drop)
    kt, dt = filter_fasta(transcript_fa, tmp_tx, drop)
    kc, dc = filter_fasta(cds_fa, tmp_cds, drop)
    os.replace(tmp_prot, protein_fa)
    os.replace(tmp_tx, transcript_fa)
    os.replace(tmp_cds, cds_fa)

    kg, dg = filter_gff(in_gff, out_gff, drop)
    print(f"[INFO] FASTA kept/dropped headers: protein {kp}/{dp}, transcript {kt}/{dt}, cds {kc}/{dc}")
    print(f"[INFO] GFF3 features kept/dropped lines: {kg}/{dg}")
    print(f"[INFO] QC GFF3: {out_gff}")
    print(f"[INFO] Dropped list: {dropped_list}")

    open(marker, "w").write("ok\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
