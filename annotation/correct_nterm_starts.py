#!/usr/bin/env python3
"""
Correct late CDS start ATGs using ZH13v3 / Wm82 ortholog evidence.

Logic (approach A):
  1. Read diamond blastp vs ZH13 and Wm82 (precomputed).
  2. Flag transcripts where query aligns from ~N-term but subject starts mid-protein
     (sstart >> 1) with high identity — classic N-terminal truncation.
  3. On the genome, search upstream in-frame ATGs with no intervening stop.
  4. Prefer the ATG whose extension length matches missing subject N-term and
     whose translated peptide matches the subject N-terminus.
  5. Rewrite GFF (gene/mRNA/exon/CDS/start_codon) and write evidence TSV.

Environment / CLI via argparse.
"""
from __future__ import annotations

import argparse
import os
import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

TABLE = {
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


def parse_attr(s: str) -> Dict[str, str]:
    d = {}
    for item in s.strip().split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def translate(cds: str) -> str:
    cds = cds.upper().replace("U", "T")
    aa = []
    for i in range(0, len(cds) - 2, 3):
        codon = cds[i : i + 3]
        if len(codon) < 3:
            break
        a = TABLE.get(codon, "X")
        if a == "*":
            break
        aa.append(a)
    return "".join(aa)


def read_fasta(path: str) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    name = None
    chunks: List[str] = []
    with open(path) as f:
        for line in f:
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


def faidx_fetch(genome: Dict[str, str], chrom: str, start: int, end: int) -> str:
    """1-based inclusive."""
    if chrom not in genome:
        # try common aliases
        for k in genome:
            if k == chrom or k.split()[0] == chrom:
                chrom = k
                break
        else:
            raise KeyError(f"chrom {chrom} not in genome")
    if start < 1:
        start = 1
    seq = genome[chrom]
    if end > len(seq):
        end = len(seq)
    return seq[start - 1 : end]


@dataclass
class Hit:
    ref_db: str
    ref_id: str
    pident: float
    length: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    qlen: int
    slen: int


def load_blast(
    path: str,
    ref_db: str,
    proteins: Dict[str, str],
    refs: Dict[str, str],
    min_pident: float,
    min_missing_aa: int = 15,
) -> Dict[str, Hit]:
    """Keep best bitscore hit per query that looks like N-term truncation."""
    best: Dict[str, Hit] = {}
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return best
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 12:
                continue
            qid, sid = p[0], p[1]
            pident = float(p[2])
            length = int(p[3])
            qstart, qend = int(p[6]), int(p[7])
            sstart, send = int(p[8]), int(p[9])
            evalue = float(p[10])
            bitscore = float(p[11])
            if pident < min_pident:
                continue
            qlen = len(proteins.get(qid, ""))
            slen = len(refs.get(sid, ""))
            if qlen == 0 or slen == 0:
                continue
            # subject N-term missing relative to query alignment start
            # require query near its own N-term
            if qstart > 5:
                continue
            if sstart - 1 < min_missing_aa:
                continue
            # alignment should cover most of the (short) query
            if length < max(50, int(0.7 * qlen)):
                continue
            hit = Hit(
                ref_db=ref_db,
                ref_id=sid,
                pident=pident,
                length=length,
                qstart=qstart,
                qend=qend,
                sstart=sstart,
                send=send,
                evalue=evalue,
                bitscore=bitscore,
                qlen=qlen,
                slen=slen,
            )
            prev = best.get(qid)
            if prev is None or hit.bitscore > prev.bitscore:
                best[qid] = hit
    return best


@dataclass
class Feat:
    chrom: str
    source: str
    ftype: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attr: str
    d: Dict[str, str]
    line_no: int


def load_gff(path: str) -> List[Feat]:
    feats = []
    with open(path) as f:
        for i, line in enumerate(f):
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            d = parse_attr(p[8])
            feats.append(
                Feat(
                    chrom=p[0],
                    source=p[1],
                    ftype=p[2],
                    start=int(p[3]),
                    end=int(p[4]),
                    score=p[5],
                    strand=p[6],
                    phase=p[7],
                    attr=p[8],
                    d=d,
                    line_no=i,
                )
            )
    return feats


def identity(a: str, b: str) -> float:
    n = min(len(a), len(b))
    if n == 0:
        return 0.0
    return sum(1 for i in range(n) if a[i] == b[i]) / n


def find_extension(
    genome: Dict[str, str],
    chrom: str,
    strand: str,
    cds_start: int,
    cds_end: int,
    missing_aa: int,
    ref_nterm: str,
    max_extra_bp: int = 90,
    min_nterm_ident: float = 0.85,
) -> Optional[Tuple[int, int, str, float]]:
    """
    Return (new_cds_start, new_cds_end, extension_aa, nterm_identity) or None.
    Only shifts the 5' end; 3' end unchanged.
    """
    need_bp = missing_aa * 3
    search_bp = need_bp + max_extra_bp

    if strand == "+":
        win_start = max(1, cds_start - search_bp)
        win_end = cds_start - 1
        if win_end < win_start:
            return None
        upstream = faidx_fetch(genome, chrom, win_start, win_end)
        # positions in upstream corresponding to genomic coords
        # try ATGs in same frame as current start: (cds_start - atg) % 3 == 0
        candidates = []
        for i in range(0, len(upstream) - 2):
            if upstream[i : i + 3] != "ATG":
                continue
            atg = win_start + i
            if (cds_start - atg) % 3 != 0:
                continue
            if atg >= cds_start:
                continue
            ext = upstream[i:]  # ATG .. just before old start
            # translate extension only (should have no stop)
            aa_ext = []
            stop = False
            for j in range(0, len(ext) - 2, 3):
                codon = ext[j : j + 3]
                if len(codon) < 3:
                    break
                a = TABLE.get(codon, "X")
                if a == "*":
                    stop = True
                    break
                aa_ext.append(a)
            if stop:
                continue
            pep = "".join(aa_ext)
            # compare to ref N-term of length len(pep)
            ref_part = ref_nterm[: len(pep)]
            if len(ref_part) < max(5, missing_aa - 5):
                # too short relative to expected missing
                if abs(len(pep) - missing_aa) > 5:
                    continue
            ident = identity(pep, ref_part) if ref_part else 0.0
            # prefer length close to missing_aa and high identity
            length_pen = abs(len(pep) - missing_aa)
            candidates.append((ident, -length_pen, atg, pep))
        if not candidates:
            return None
        candidates.sort(reverse=True)
        ident, _, atg, pep = candidates[0]
        if ident < min_nterm_ident:
            return None
        # require roughly expected length
        if abs(len(pep) - missing_aa) > max(5, int(0.25 * missing_aa) + 1):
            return None
        return atg, cds_end, pep, ident

    # strand == "-"
    win_start = cds_end + 1
    win_end = cds_end + search_bp
    downstream = faidx_fetch(genome, chrom, win_start, win_end)  # genomic + direction
    # on - strand, CDS 5' is at cds_end; upstream in transcript sense is higher coords
    candidates = []
    for i in range(0, len(downstream) - 2):
        # genomic ATG on - strand is reverse-complement of codon at that position
        # codon occupying genomic [pos, pos+2] on - strand is revcomp(seq[pos:pos+3])
        # For start, ATG in transcript sense means genomic triplet revcomp == ATG
        # i.e. genomic == CAT (revcomp of ATG)
        # Actually: reading on - strand from cds_end toward higher coords:
        # transcript 5' extension is revcomp of genomic[cds_end+1 : new_end]
        pass

    # Simpler for - strand: for each candidate new_end > cds_end with (new_end-cds_end)%3==0
    # take genomic segment [cds_start, new_end], revcomp, translate from start
    for ext_bp in range(3, search_bp + 1, 3):
        new_end = cds_end + ext_bp
        try:
            seg = faidx_fetch(genome, chrom, cds_start, new_end)
        except Exception:
            break
        cds = revcomp(seg)
        if not cds.startswith("ATG"):
            continue
        # check no stop in the extension part only (first ext_bp bases)
        ext = cds[:ext_bp]
        aa_ext = []
        stop = False
        for j in range(0, len(ext) - 2, 3):
            codon = ext[j : j + 3]
            a = TABLE.get(codon, "X")
            if a == "*":
                stop = True
                break
            aa_ext.append(a)
        if stop:
            continue
        pep = "".join(aa_ext)
        ref_part = ref_nterm[: len(pep)]
        ident = identity(pep, ref_part) if ref_part else 0.0
        length_pen = abs(len(pep) - missing_aa)
        candidates.append((ident, -length_pen, new_end, pep))

    if not candidates:
        return None
    candidates.sort(reverse=True)
    ident, _, new_end, pep = candidates[0]
    if ident < min_nterm_ident:
        return None
    if abs(len(pep) - missing_aa) > max(5, int(0.25 * missing_aa) + 1):
        return None
    return cds_start, new_end, pep, ident


def apply_fix(
    feats: List[Feat],
    tid: str,
    gid: str,
    strand: str,
    old_start: int,
    old_end: int,
    new_start: int,
    new_end: int,
) -> None:
    """In-place update of features belonging to transcript/gene."""
    for ft in feats:
        fid = ft.d.get("ID", "")
        parents = [x for x in ft.d.get("Parent", "").split(",") if x]

        if ft.ftype == "gene" and fid == gid:
            if strand == "+":
                ft.start = min(ft.start, new_start)
            else:
                ft.end = max(ft.end, new_end)
            continue

        if ft.ftype in ("mRNA", "transcript") and fid == tid:
            if strand == "+":
                ft.start = min(ft.start, new_start)
            else:
                ft.end = max(ft.end, new_end)
            continue

        if tid not in parents:
            continue

        if ft.ftype in ("exon", "CDS"):
            if strand == "+":
                if ft.start == old_start:
                    ft.start = new_start
                    if ft.ftype == "CDS":
                        ft.phase = "0"
            else:
                if ft.end == old_end:
                    ft.end = new_end
                    if ft.ftype == "CDS":
                        ft.phase = "0"
        elif ft.ftype == "start_codon":
            if strand == "+":
                ft.start = new_start
                ft.end = new_start + 2
            else:
                ft.end = new_end
                ft.start = new_end - 2
        elif ft.ftype == "five_prime_UTR":
            if strand == "+":
                if ft.start >= new_start:
                    ft.start = 0  # delete
                elif ft.end >= new_start:
                    ft.end = new_start - 1
            else:
                if ft.end <= new_end:
                    ft.start = 0
                elif ft.start <= new_end:
                    ft.start = new_end + 1


def write_gff(path: str, header_lines: List[str], feats: List[Feat]) -> None:
    with open(path, "w") as w:
        for h in header_lines:
            w.write(h)
        for ft in feats:
            if ft.start <= 0 or ft.end < ft.start:
                continue
            w.write(
                f"{ft.chrom}\t{ft.source}\t{ft.ftype}\t{ft.start}\t{ft.end}\t"
                f"{ft.score}\t{ft.strand}\t{ft.phase}\t{ft.attr}\n"
            )


def main() -> int:
    ap = argparse.ArgumentParser(description="Ortholog-guided CDS start correction")
    ap.add_argument("--genome", required=True)
    ap.add_argument("--gff-in", required=True)
    ap.add_argument("--gff-out", required=True)
    ap.add_argument("--proteins", required=True, help="query protein fasta (gffread -y)")
    ap.add_argument("--zh13-fa", required=True)
    ap.add_argument("--wm82-fa", required=True)
    ap.add_argument("--zh13-blast", required=True)
    ap.add_argument("--wm82-blast", required=True)
    ap.add_argument("--evidence-tsv", required=True)
    ap.add_argument("--min-pident", type=float, default=95.0)
    ap.add_argument("--min-nterm-ident", type=float, default=0.95)
    ap.add_argument("--min-missing-aa", type=int, default=15)
    args = ap.parse_args()

    print(f"[INFO] Loading genome: {args.genome}", flush=True)
    genome = read_fasta(args.genome)
    print(f"[INFO] Loading query proteins: {args.proteins}", flush=True)
    proteins = read_fasta(args.proteins)
    print(f"[INFO] Loading ZH13 proteins: {args.zh13_fa}", flush=True)
    zh13 = read_fasta(args.zh13_fa)
    print(f"[INFO] Loading Wm82 proteins: {args.wm82_fa}", flush=True)
    wm82 = read_fasta(args.wm82_fa)

    hits_zh = load_blast(
        args.zh13_blast, "ZH13v3", proteins, zh13, args.min_pident, args.min_missing_aa
    )
    hits_wm = load_blast(
        args.wm82_blast, "Wm82_gnm6", proteins, wm82, args.min_pident, args.min_missing_aa
    )

    # merge: prefer higher bitscore across refs
    candidates: Dict[str, Hit] = {}
    for qid, hit in list(hits_zh.items()) + list(hits_wm.items()):
        prev = candidates.get(qid)
        if prev is None or hit.bitscore > prev.bitscore:
            candidates[qid] = hit

    print(f"[INFO] N-term truncation candidates: {len(candidates)}", flush=True)

    # parse GFF with headers
    header_lines: List[str] = []
    with open(args.gff_in) as f:
        for line in f:
            if line.startswith("#"):
                header_lines.append(line)
            else:
                break
    feats = load_gff(args.gff_in)

    # index CDS by transcript Parent
    cds_by_tid: Dict[str, List[Feat]] = defaultdict(list)
    mrna_by_id: Dict[str, Feat] = {}
    gene_of_mrna: Dict[str, str] = {}
    for ft in feats:
        if ft.ftype in ("mRNA", "transcript"):
            mrna_by_id[ft.d.get("ID", "")] = ft
            gene_of_mrna[ft.d.get("ID", "")] = ft.d.get("Parent", "")
        if ft.ftype == "CDS":
            for p in ft.d.get("Parent", "").split(","):
                if p:
                    cds_by_tid[p].append(ft)

    evidence_rows = []
    n_corrected = 0

    for tid, hit in sorted(candidates.items()):
        cds_list = sorted(cds_by_tid.get(tid, []), key=lambda x: (x.start, x.end))
        if not cds_list:
            evidence_rows.append(
                {
                    "transcript_id": tid,
                    "status": "skipped_no_cds",
                    "note": "no CDS features",
                    "hit": hit,
                }
            )
            continue
        strand = cds_list[0].strand
        chrom = cds_list[0].chrom
        if strand == "+":
            old_start = min(c.start for c in cds_list)
            old_end = max(c.end for c in cds_list)
        else:
            old_start = min(c.start for c in cds_list)
            old_end = max(c.end for c in cds_list)

        missing_aa = hit.sstart - 1
        ref_db = zh13 if hit.ref_db == "ZH13v3" else wm82
        ref_seq = ref_db.get(hit.ref_id, "")
        ref_nterm = ref_seq[:missing_aa]
        if len(ref_nterm) < 5:
            evidence_rows.append(
                {
                    "transcript_id": tid,
                    "status": "skipped_short_ref_nterm",
                    "note": f"missing_aa={missing_aa}",
                    "hit": hit,
                }
            )
            continue

        try:
            res = find_extension(
                genome,
                chrom,
                strand,
                old_start,
                old_end,
                missing_aa,
                ref_nterm,
                min_nterm_ident=args.min_nterm_ident,
            )
        except KeyError as e:
            evidence_rows.append(
                {
                    "transcript_id": tid,
                    "status": "skipped_genome_chrom",
                    "note": str(e),
                    "hit": hit,
                }
            )
            continue

        if res is None:
            evidence_rows.append(
                {
                    "transcript_id": tid,
                    "gene_id": gene_of_mrna.get(tid, ""),
                    "chrom": chrom,
                    "strand": strand,
                    "status": "skipped_no_upstream_atg",
                    "note": f"missing_aa={missing_aa}; no matching in-frame ATG",
                    "hit": hit,
                    "old_cds_start": old_start,
                    "old_cds_end": old_end,
                }
            )
            continue

        new_start, new_end, pep, ident = res
        gid = gene_of_mrna.get(tid, "")
        apply_fix(feats, tid, gid, strand, old_start, old_end, new_start, new_end)
        n_corrected += 1
        evidence_rows.append(
            {
                "transcript_id": tid,
                "gene_id": gid,
                "chrom": chrom,
                "strand": strand,
                "status": "corrected",
                "note": "upstream in-frame ATG matches ortholog N-term",
                "hit": hit,
                "old_cds_start": old_start,
                "old_cds_end": old_end,
                "new_cds_start": new_start,
                "new_cds_end": new_end,
                "extension_aa": pep,
                "extension_bp": abs(
                    (new_start - old_start) if strand == "+" else (new_end - old_end)
                ),
                "nterm_identity": f"{ident:.4f}",
                "missing_aa_expected": missing_aa,
                "upstream_atg": new_start if strand == "+" else new_end - 2,
            }
        )

    write_gff(args.gff_out, header_lines, feats)

    cols = [
        "transcript_id",
        "gene_id",
        "chrom",
        "strand",
        "status",
        "ref_db",
        "ref_id",
        "pident",
        "qstart",
        "qend",
        "sstart",
        "send",
        "qlen",
        "slen",
        "bitscore",
        "evalue",
        "missing_aa_expected",
        "old_cds_start",
        "old_cds_end",
        "new_cds_start",
        "new_cds_end",
        "upstream_atg",
        "extension_bp",
        "extension_aa",
        "nterm_identity",
        "note",
    ]

    def cell(row, k, default=""):
        v = row.get(k, default)
        return "" if v is None else str(v)

    with open(args.evidence_tsv, "w") as w:
        w.write("\t".join(cols) + "\n")
        for row in evidence_rows:
            hit = row.get("hit")
            if hit is None:
                line = [
                    cell(row, "transcript_id"),
                    cell(row, "gene_id"),
                    cell(row, "chrom"),
                    cell(row, "strand"),
                    cell(row, "status"),
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    cell(row, "missing_aa_expected"),
                    cell(row, "old_cds_start"),
                    cell(row, "old_cds_end"),
                    cell(row, "new_cds_start"),
                    cell(row, "new_cds_end"),
                    cell(row, "upstream_atg"),
                    cell(row, "extension_bp"),
                    cell(row, "extension_aa"),
                    cell(row, "nterm_identity"),
                    cell(row, "note"),
                ]
            else:
                line = [
                    cell(row, "transcript_id"),
                    cell(row, "gene_id"),
                    cell(row, "chrom"),
                    cell(row, "strand"),
                    cell(row, "status"),
                    hit.ref_db,
                    hit.ref_id,
                    f"{hit.pident:.2f}",
                    str(hit.qstart),
                    str(hit.qend),
                    str(hit.sstart),
                    str(hit.send),
                    str(hit.qlen),
                    str(hit.slen),
                    f"{hit.bitscore:.1f}",
                    str(hit.evalue),
                    cell(row, "missing_aa_expected", str(hit.sstart - 1)),
                    cell(row, "old_cds_start"),
                    cell(row, "old_cds_end"),
                    cell(row, "new_cds_start"),
                    cell(row, "new_cds_end"),
                    cell(row, "upstream_atg"),
                    cell(row, "extension_bp"),
                    cell(row, "extension_aa"),
                    cell(row, "nterm_identity"),
                    cell(row, "note"),
                ]
            w.write("\t".join(line) + "\n")

    print(f"[INFO] Corrected transcripts: {n_corrected}", flush=True)
    print(f"[INFO] Evidence: {args.evidence_tsv}", flush=True)
    print(f"[INFO] GFF out: {args.gff_out}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
