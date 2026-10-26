#!/usr/bin/env python3
"""Revise short high-confidence models in place inside 13_final.

13_final keeps one GFF and one set of FASTA files. A proposed CDS replaces
the current model only when it has an ATG, a terminal stop, no internal stop,
and a length divisible by 3. expressed_zh13_below_50 is not revised: those
genes are defined by ZH13 v3 identity below 50, so a full ZH13 protein is not
evidence for rewriting them. Locus IDs are not renumbered.
"""

import argparse
import bisect
import csv
import os
import re
import subprocess
import sys
from collections import defaultdict

STOP = {"TAA", "TAG", "TGA"}
ZH13_GENE = re.compile(r"(SoyZH13V3_\d{2}G\d+)")
TX_SUFFIX = re.compile(r"\.[mt]\d+$")
NUC = re.compile(r"^Chr(0[1-9]|1[0-9]|20)$")
STAMP = "结构修订已写入"
RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CODON = {
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


def fail(msg):
    print(f"[ERROR] {msg}", file=sys.stderr)
    sys.exit(1)


def attr(text, key):
    match = re.search(rf"(?:^|;){key}=([^;\s]+)", text)
    return match.group(1) if match else ""


def zh13_gene(name):
    match = ZH13_GENE.search(name or "")
    return match.group(1) if match else ""


def gene_of_tx(tid):
    return TX_SUFFIX.sub("", tid)


def tx_key(tid):
    match = re.search(r"\.m(\d+)$", tid)
    return int(match.group(1)) if match else 10**9


def locus_key(gid):
    match = re.search(r"_(\d+)G(\d+)$", gid)
    if not match:
        return (999, 10**12, gid)
    return (int(match.group(1)), int(match.group(2)), gid)


def fasta_lengths(path):
    lengths = {}
    name = None
    size = 0
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name:
                    lengths[name] = size
                name = line[1:].split()[0]
                size = 0
            else:
                size += len(line.strip())
        if name:
            lengths[name] = size
    return lengths


def fasta_seqs(path, names):
    wanted = set(names)
    seqs = {}
    name = None
    keep = False
    chunks = []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if keep and name:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                keep = name in wanted
                chunks = []
            elif keep:
                chunks.append(line.strip())
        if keep and name:
            seqs[name] = "".join(chunks)
    return seqs


def load_genome(path):
    seqs = {}
    name = None
    chunks = []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if name:
            seqs[name] = "".join(chunks).upper()
    return seqs


def rc(seq):
    return seq.translate(RC_TABLE)[::-1]


def classify_cds(seq):
    seq = (seq or "").upper().replace("U", "T")
    if not seq:
        return ["无CDS"]
    if any(base not in "ACGT" for base in seq):
        return ["含非ACGT碱基"]
    if len(seq) % 3 != 0:
        return ["CDS长度不能被3整除"]
    reasons = []
    if not seq.startswith("ATG"):
        reasons.append("缺少起始密码子ATG")
    codons = [seq[i:i + 3] for i in range(0, len(seq), 3)]
    if any(codon in STOP for codon in codons[:-1]):
        reasons.append("内部终止密码子")
    if not codons or codons[-1] not in STOP:
        reasons.append("缺少终止密码子")
    return reasons


def translate(seq):
    """Translate a CDS, dropping a terminal stop when one is present."""
    aa = []
    codons = [seq[i:i + 3] for i in range(0, len(seq), 3)]
    if codons and codons[-1] in STOP:
        codons = codons[:-1]
    for codon in codons:
        aa.append(CODON.get(codon, "X"))
    return "".join(aa)


def splice(cds_parts, chrom_seq):
    """Concatenate miniprot CDS intervals in transcription order.

    Miniprot writes a phase column, but the CDS intervals already contain the
    coding bases, including split codons. Skipping those phase bases drops
    real codon sequence. The revised GFF therefore records phase 0 and the
    FASTA is this concatenation.
    """
    strand = cds_parts[0][3]
    ordered = sorted(cds_parts, key=lambda item: item[0], reverse=(strand == "-"))
    parts = []
    for start, end, _phase, part_strand in ordered:
        piece = chrom_seq[start - 1:end]
        if part_strand == "-":
            piece = rc(piece)
        parts.append(piece)
    return "".join(parts).upper()


def merge_intervals(intervals):
    if not intervals:
        return []
    ordered = sorted(intervals)
    out = [list(ordered[0])]
    for start, end in ordered[1:]:
        if start <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], end)
        else:
            out.append([start, end])
    return [(start, end) for start, end in out]


def intersect_len(left, right):
    i = j = total = 0
    while i < len(left) and j < len(right):
        start = max(left[i][0], right[j][0])
        end = min(left[i][1], right[j][1])
        if start <= end:
            total += end - start + 1
        if left[i][1] < right[j][1]:
            i += 1
        else:
            j += 1
    return total


def wrap(seq):
    return "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60))


def load_banned(renamed, prefix):
    banned = set()
    drop_tsv = os.path.join(renamed, f"{prefix}.dropped.tsv")
    drop_list = os.path.join(renamed, f"{prefix}.dropped_cds_nonmod3.genes.list")
    if os.path.isfile(drop_tsv) and os.path.getsize(drop_tsv):
        with open(drop_tsv) as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row.get("删除范围") == "基因" and row.get("基因ID"):
                    banned.add(row["基因ID"].strip())
    if os.path.isfile(drop_list) and os.path.getsize(drop_list):
        with open(drop_list) as handle:
            for line in handle:
                gene = line.strip()
                if gene:
                    banned.add(gene)
    return banned


def load_mapping(path):
    mapping = {}
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            mapping[row["old_gene_id"]] = row["new_gene_id"]
    return mapping


def tx_to_gene(tx, mapping):
    gene = "evm.TU." + tx[len("evm.model."):]
    if gene in mapping:
        return mapping[gene]
    gene = re.sub(r"\.\d+$", "", gene)
    return mapping.get(gene)


def load_hc(path):
    genes = {}
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["class"] == "high_confidence":
                genes[row["gene_id"]] = row
    return genes


def load_zh13_hits(path):
    hits = {}
    key = "ZH13v3.0_mRNA_nuclear_20250729_pro_hit"
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            hits[row["gene_id"]] = zh13_gene(row.get(key, ""))
    return hits


def select_short(down, prefix, hc, zh13_len):
    renamed = os.path.join(down, "12_renamed")
    mapping = load_mapping(os.path.join(renamed, f"{prefix}.gene_id_mapping.tsv"))
    query_len = fasta_lengths(os.path.join(down, "02_sequences", f"{prefix}.protein.fa"))
    diamond = os.path.join(
        down, "05_diamond", f"{prefix}.vs_ZH13v3.0_mRNA_nuclear_20250729_pro.tsv"
    )
    if not os.path.isfile(diamond):
        fail(f"missing diamond table: {diamond}")
    best = {}
    with open(diamond) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            bits = float(parts[11])
            old = best.get(parts[0])
            if old is None or bits > old[0]:
                best[parts[0]] = (
                    bits, parts[1], float(parts[2]),
                    int(parts[6]), int(parts[7]), int(parts[8]), int(parts[9]),
                )
    per_gene = {}
    for tx, rec in best.items():
        gid = tx_to_gene(tx, mapping)
        if not gid or gid not in hc:
            continue
        old = per_gene.get(gid)
        if old is None or rec[0] > old[0][0]:
            per_gene[gid] = (rec, tx)
    short = {}
    for gid, (rec, tx) in per_gene.items():
        _bits, sid, pid, qstart, qend, sstart, send = rec
        qlen = query_len.get(tx)
        slen = zh13_len.get(sid)
        if not qlen or not slen:
            continue
        qcov = (abs(qend - qstart) + 1) / qlen
        scov = (abs(send - sstart) + 1) / slen
        if qcov >= 0.5 and scov < 0.5:
            chrom, start, end, strand = (
                hc[gid]["chrom"], int(hc[gid]["start"]), int(hc[gid]["end"]), hc[gid]["strand"]
            )
            short[gid] = {
                "chrom": chrom, "start": start, "end": end, "strand": strand,
                "zh13_id": sid, "zh13_gene": zh13_gene(sid),
                "query_aa": qlen, "zh13_aa": slen,
                "qcov": qcov, "scov": scov, "pident": pid,
            }
    return short


def ensure_miniprot(args, short, zh13_fasta, zh13_len):
    if args.miniprot_gff:
        if not os.path.isfile(args.miniprot_gff):
            fail(f"missing miniprot gff: {args.miniprot_gff}")
        return args.miniprot_gff
    work = os.path.join(args.downstream, "tmp_curate")
    os.makedirs(work, exist_ok=True)
    gff = os.path.join(work, f"{args.prefix}.miniprot.gff")
    queries = os.path.join(work, "zh13_queries.faa")
    stamp = os.path.join(work, "query_ids.txt")
    wanted = sorted({row["zh13_id"] for row in short.values()})
    previous = []
    if os.path.isfile(stamp):
        with open(stamp) as handle:
            previous = [line.strip() for line in handle if line.strip()]
    if wanted == previous and os.path.isfile(gff) and os.path.getsize(gff) > 0 and not args.force:
        print(f"[INFO] reuse miniprot gff: {gff}", flush=True)
        return gff
    if not args.miniprot or not os.path.isfile(args.miniprot):
        fail(f"miniprot not found: {args.miniprot}")
    have = fasta_seqs(zh13_fasta, wanted)
    missing = [name for name in wanted if name not in have]
    if missing:
        fail(f"{len(missing)} ZH13 proteins missing from the reference fasta, e.g. {missing[:3]}")
    with open(queries, "w") as handle:
        for name in wanted:
            handle.write(f">{name}\n{wrap(have[name])}\n")
    print(f"[INFO] miniprot {len(wanted)} proteins", flush=True)
    log_path = os.path.join(work, "miniprot.log")
    cmd = [
        args.miniprot, "--gff", "-t", str(args.threads),
        "--outs=0", "--outc=0.5", "-N", "30", args.genome, queries,
    ]
    with open(gff, "w") as stdout, open(log_path, "w") as stderr:
        subprocess.run(cmd, stdout=stdout, stderr=stderr, check=True)
    with open(stamp, "w") as handle:
        handle.write("\n".join(wanted) + "\n")
    return gff


def parse_miniprot(path):
    mrnas = {}
    order = []
    pending = None
    current = None
    with open(path) as handle:
        for line in handle:
            if line.startswith("##PAF"):
                tags = {}
                for field in line.rstrip("\n").split("\t"):
                    if field.startswith("fs:i:") or field.startswith("st:i:"):
                        key, value = field.split(":i:")
                        tags[key] = int(value)
                pending = tags
                current = None
                continue
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[2] == "mRNA":
                mid = attr(parts[8], "ID")
                target = re.search(r"Target=(\S+)\s+(\d+)\s+(\d+)", parts[8])
                if not mid or not target:
                    current = None
                    pending = None
                    continue
                rec = {
                    "id": mid, "chrom": parts[0], "start": int(parts[3]), "end": int(parts[4]),
                    "strand": parts[6], "ident": float(attr(parts[8], "Identity") or 0),
                    "target": target.group(1), "t1": int(target.group(2)), "t2": int(target.group(3)),
                    "fs": None if pending is None else pending.get("fs"),
                    "st": None if pending is None else pending.get("st"),
                    "cds": [], "stops": [],
                }
                mrnas[mid] = rec
                order.append(rec)
                current = rec
                pending = None
            elif current is not None and parts[2] == "CDS" and attr(parts[8], "Parent") == current["id"]:
                phase = int(parts[7]) if parts[7] in {"0", "1", "2"} else 0
                current["cds"].append((int(parts[3]), int(parts[4]), phase, parts[6]))
            elif current is not None and parts[2] == "stop_codon" and attr(parts[8], "Parent") == current["id"]:
                current["stops"].append((int(parts[3]), int(parts[4])))
    by_chrom = defaultdict(list)
    for rec in order:
        by_chrom[rec["chrom"]].append(rec)
    for chrom in by_chrom:
        by_chrom[chrom].sort(key=lambda rec: rec["start"])
    return mrnas, by_chrom


def load_gff(path):
    header = []
    blocks = []
    current = None
    tx_of = defaultdict(list)
    cds_of = defaultdict(list)
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                header.append(line if line.endswith("\n") else line + "\n")
                continue
            if not line.strip():
                continue
            if not line.endswith("\n"):
                line += "\n"
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 9 and parts[2] == "gene":
                if current:
                    blocks.append(current)
                gid = attr(parts[8], "ID")
                current = [gid, [line]]
            elif current:
                current[1].append(line)
                if parts[2] == "mRNA":
                    tx_of[current[0]].append(attr(parts[8], "ID"))
                elif parts[2] == "CDS":
                    cds_of[current[0]].append((int(parts[3]), int(parts[4]), parts[6]))
            else:
                header.append(line)
        if current:
            blocks.append(current)
    return header, blocks, tx_of, cds_of


def gene_index(blocks, cds_of):
    spans = {}
    merged = {}
    on_chrom = defaultdict(list)
    for gid, lines in blocks:
        parts = lines[0].rstrip("\n").split("\t")
        spans[gid] = (parts[0], int(parts[3]), int(parts[4]), parts[6])
        iv = merge_intervals([(start, end) for start, end, _strand in cds_of.get(gid, [])])
        merged[gid] = iv
        on_chrom[parts[0]].append((int(parts[3]), int(parts[4]), parts[6], gid))
    for chrom in on_chrom:
        on_chrom[chrom].sort()
    max_len = 1
    for chrom, start, end, _strand in spans.values():
        max_len = max(max_len, end - start + 1)
    return spans, merged, on_chrom, max_len


def other_overlaps(index, chrom, strand, intervals, allowed):
    _spans, merged, on_chrom, max_len = index
    start0 = min(start for start, _end in intervals)
    end0 = max(end for _start, end in intervals)
    rows = on_chrom.get(chrom, [])
    starts = [row[0] for row in rows]
    pos = bisect.bisect_left(starts, start0 - max_len)
    query = merge_intervals(intervals)
    hits = []
    while pos < len(rows) and rows[pos][0] <= end0:
        gstart, gend, gstrand, gid = rows[pos]
        if gstrand == strand and gend >= start0 and gid not in allowed:
            overlap = intersect_len(query, merged.get(gid, []))
            if overlap >= 200:
                hits.append((gid, overlap))
        pos += 1
    return hits


def choose_alignment(gene, by_chrom, zh13_len):
    best = None
    glen = gene["end"] - gene["start"] + 1
    for rec in by_chrom.get(gene["chrom"], []):
        if rec["strand"] != gene["strand"] or not rec["cds"]:
            continue
        if zh13_gene(rec["target"]) != gene["zh13_gene"]:
            continue
        overlap = min(gene["end"], rec["end"]) - max(gene["start"], rec["start"]) + 1
        if overlap < 0.5 * glen:
            continue
        slen = zh13_len.get(rec["target"])
        if not slen:
            continue
        qcov = (rec["t2"] - rec["t1"] + 1) / slen
        if qcov < 0.8:
            continue
        rank = (qcov, rec["ident"])
        if best is None or rank > best[0]:
            best = (rank, rec, qcov)
    return None if best is None else (best[1], best[2])


def cds_sequence(rec, genome):
    chrom_seq = genome.get(rec["chrom"])
    if not chrom_seq:
        return ""
    seq = splice(rec["cds"], chrom_seq)
    if not seq:
        return ""
    if seq[-3:] not in STOP and rec["stops"]:
        strand = rec["strand"]
        stops = sorted(rec["stops"], reverse=(strand == "-"))
        extra = []
        for start, end in stops:
            piece = chrom_seq[start - 1:end]
            extra.append(rc(piece) if strand == "-" else piece)
        tail = "".join(extra).upper()
        if tail in STOP and not seq.endswith(tail):
            seq += tail
    return seq


def render_gene(gid, txid, rec):
    coords = [(start, end) for start, end, _phase, _strand in rec["cds"]]
    coords.extend(rec["stops"])
    start = min(item[0] for item in coords)
    end = max(item[1] for item in coords)
    chrom = rec["chrom"]
    strand = rec["strand"]
    note = f"Note=miniprot:{rec['target']}"
    lines = [
        f"{chrom}\tEVM\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gid}\n",
        f"{chrom}\tEVM\tmRNA\t{start}\t{end}\t.\t{strand}\t.\tID={txid};Parent={gid};{note}\n",
    ]
    for index, (cstart, cend, phase, _strand) in enumerate(sorted(rec["cds"]), start=1):
        lines.append(
            f"{chrom}\tEVM\tCDS\t{cstart}\t{cend}\t.\t{strand}\t0\t"
            f"Parent={txid};ID={txid}.cds{index}\n"
        )
        lines.append(
            f"{chrom}\tEVM\texon\t{cstart}\t{cend}\t.\t{strand}\t.\t"
            f"Parent={txid};ID={txid}.exon{index}\n"
        )
    return lines


def self_test(mrnas, genome, zh13_fasta, zh13_len):
    full = []
    for rec in mrnas.values():
        if not NUC.match(rec["chrom"]) or not rec["cds"] or rec["ident"] < 1:
            continue
        slen = zh13_len.get(rec["target"])
        if slen and rec["t1"] == 1 and rec["t2"] == slen and rec["fs"] == 0 and rec["st"] == 0:
            full.append(rec)
    plus = [rec for rec in full if rec["strand"] == "+"][:20]
    minus = [rec for rec in full if rec["strand"] == "-"][:20]
    sample = plus + minus
    if len(sample) < 5:
        fail(f"miniprot self-test found only {len(sample)} full-length alignments")
    seqs = fasta_seqs(zh13_fasta, [rec["target"] for rec in sample])
    bad = []
    for rec in sample:
        cds = cds_sequence(rec, genome)
        protein = translate(cds) if cds and len(cds) % 3 == 0 else ""
        ref = seqs.get(rec["target"], "").replace("*", "")
        internal = any(codon in STOP for codon in (
            cds[i:i + 3] for i in range(0, max(0, len(cds) - 3), 3)
        ))
        if not protein or protein != ref or internal:
            bad.append((rec["id"], rec["target"], rec["strand"], cds[-3:], len(protein), len(ref)))
            if len(bad) >= 5:
                break
    if bad:
        fail(f"miniprot CDS translation does not match ZH13: {bad}")
    print(f"[INFO] miniprot translation check passed on {len(sample)} full-length proteins", flush=True)


def stream_fasta(src, dst, drop_ids, replace):
    seen = set()
    with open(src) as inp, open(dst, "w") as out:
        mode = "copy"
        for line in inp:
            if line.startswith(">"):
                tid = line[1:].split()[0]
                if tid in drop_ids:
                    mode = "skip"
                elif tid in replace:
                    mode = "skip"
                    seen.add(tid)
                    out.write(f">{tid}\n{wrap(replace[tid])}\n")
                else:
                    mode = "copy"
                    out.write(line if line.endswith("\n") else line + "\n")
            elif mode == "copy":
                out.write(line if line.endswith("\n") else line + "\n")
    missing = set(replace) - seen
    if missing:
        fail(f"FASTA {src} is missing transcript IDs to replace, e.g. {sorted(missing)[:3]}")


def write_text(path, text):
    with open(path, "w") as handle:
        handle.write(text)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--downstream", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--zh13-fasta", required=True)
    parser.add_argument("--miniprot", default="")
    parser.add_argument("--miniprot-gff", default="")
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    final_dir = os.path.join(args.downstream, "13_final")
    renamed = os.path.join(args.downstream, "12_renamed")
    gff_path = os.path.join(final_dir, f"{args.prefix}.final.gff3")
    class_path = os.path.join(final_dir, f"{args.prefix}.final.gene_class.tsv")
    criteria_path = os.path.join(final_dir, f"{args.prefix}.final.criteria.txt")
    for path in (gff_path, class_path, args.genome, args.zh13_fasta):
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            fail(f"missing or empty: {path}")
    if not args.dry_run and not args.force and os.path.isfile(criteria_path):
        with open(criteria_path) as handle:
            if STAMP in handle.read():
                print(f"{args.prefix}\talready_revised\tleft_unchanged", flush=True)
                return

    print(f"[INFO] {args.prefix}: load references", flush=True)
    zh13_len = fasta_lengths(args.zh13_fasta)
    hc = load_hc(class_path)
    summary_hits = load_zh13_hits(
        os.path.join(renamed, f"{args.prefix}.renamed.annotation.summary.tsv")
    )
    banned = load_banned(renamed, args.prefix)
    short = select_short(args.downstream, args.prefix, hc, zh13_len)
    print(f"[INFO] {args.prefix}: short high-confidence genes {len(short)}", flush=True)
    gff_mp = ensure_miniprot(args, short, args.zh13_fasta, zh13_len)
    _mrnas, by_chrom = parse_miniprot(gff_mp)
    print(f"[INFO] {args.prefix}: load genome", flush=True)
    genome = load_genome(args.genome)
    self_test(_mrnas, genome, args.zh13_fasta, zh13_len)

    header, blocks, tx_of, cds_of = load_gff(gff_path)
    index = gene_index(blocks, cds_of)
    in_final = {gid for gid, _lines in blocks}
    if banned & in_final:
        fail(f"dropped genes are already in 13_final, e.g. {sorted(banned & in_final)[:5]}")

    decisions = {}
    chosen = {}
    for gid, gene in short.items():
        if gid not in in_final:
            continue
        if gid in banned:
            decisions[gid] = ("保持原模型", gid, "", "已在删除清单")
            continue
        picked = choose_alignment(gene, by_chrom, zh13_len)
        if picked is None:
            decisions[gid] = ("保持原模型", gid, "", "没有同位点的完整比对")
            continue
        rec, qcov = picked
        chosen[gid] = (rec, qcov)
        if rec["fs"] != 0 or rec["st"] != 0:
            decisions[gid] = ("保持原模型", gid, rec["id"], "移码或框内终止")
            continue
        if rec["ident"] < 0.8:
            decisions[gid] = ("保持原模型", gid, rec["id"], "相似度低于0.8")
            continue
        foreign = []
        for other, overlap in other_overlaps(
            index, rec["chrom"], rec["strand"],
            [(start, end) for start, end, _phase, _strand in rec["cds"]], {gid},
        ):
            other_gene = summary_hits.get(other, "")
            if other_gene != gene["zh13_gene"]:
                foreign.append((other, overlap, other_gene))
        if foreign:
            decisions[gid] = ("保持原模型", gid, rec["id"], "与其他ZH13基因编码区重叠")
            continue
        decisions[gid] = ("clean", gid, rec["id"], "")

    clean = [gid for gid, row in decisions.items() if row[0] == "clean"]
    groups = defaultdict(list)
    for gid in clean:
        groups[short[gid]["zh13_gene"]].append(gid)

    plans = {}
    for zg, members in groups.items():
        if len(members) == 1:
            gid = members[0]
            plans[gid] = ("extend", [gid], chosen[gid][0])
            continue
        by_model = defaultdict(list)
        for gid in members:
            by_model[chosen[gid][0]["id"]].append(gid)
        for mid, cluster in by_model.items():
            if len(cluster) == 1:
                decisions[cluster[0]] = ("保持原模型", cluster[0], mid, "合并证据不足")
                continue
            plans[tuple(sorted(cluster, key=locus_key))] = ("merge", cluster, chosen[cluster[0]][0])

    tx_needed = set()
    for gid in short:
        tx_needed.update(tx_of.get(gid, []))
    proteins = fasta_seqs(
        os.path.join(final_dir, f"{args.prefix}.final.proteins.fa"), tx_needed
    )
    old_aa = {}
    for gid in short:
        lengths = [len(proteins[tid]) for tid in tx_of.get(gid, []) if tid in proteins]
        old_aa[gid] = max(lengths) if lengths else 0

    applied = {}
    reason_counts = defaultdict(int)
    for key, (kind, members, rec) in plans.items():
        seq = cds_sequence(rec, genome)
        reasons = classify_cds(seq)
        if reasons:
            for gid in members:
                decisions[gid] = ("保持原模型", members[0], rec["id"], ";".join(reasons))
            continue
        protein = translate(seq)
        if any(old_aa.get(gid, 0) >= len(protein) for gid in members):
            for gid in members:
                decisions[gid] = ("保持原模型", members[0], rec["id"], "修订后蛋白没有变长")
            continue
        allowed = set(members)
        blocked = other_overlaps(
            index, rec["chrom"], rec["strand"],
            [(start, end) for start, end, _phase, _strand in rec["cds"]], allowed,
        )
        if blocked:
            for gid in members:
                decisions[gid] = ("保持原模型", members[0], rec["id"], "修订后的编码区与保留基因重叠")
            continue
        if any(gid in banned for gid in members):
            for gid in members:
                decisions[gid] = ("保持原模型", members[0], rec["id"], "已在删除清单")
            continue
        kept = sorted(members, key=locus_key)[0]
        kept_tx = sorted(tx_of[kept], key=tx_key)[0]
        applied[kept] = {
            "kind": kind, "members": sorted(members, key=locus_key), "rec": rec,
            "tx": kept_tx, "cds": seq, "protein": protein,
        }
        if kind == "extend":
            decisions[kept] = ("延长", kept, rec["id"], "通过结构QC")
        else:
            decisions[kept] = ("合并保留", kept, rec["id"], "同一条miniprot比对")
            for gid in members:
                if gid != kept:
                    decisions[gid] = ("合并删除", kept, rec["id"], "同一条miniprot比对")

    # Each extension is checked against the original CDS of other genes. Two
    # extensions accepted in the same pass can still overlap each other.
    overlap_reject = set()
    applied_ids = list(applied)
    for left_i, left_id in enumerate(applied_ids):
        left = applied[left_id]
        left_cds = merge_intervals(
            [(start, end) for start, end, _phase, _strand in left["rec"]["cds"]]
        )
        for right_id in applied_ids[left_i + 1:]:
            right = applied[right_id]
            if left["rec"]["chrom"] != right["rec"]["chrom"]:
                continue
            if left["rec"]["strand"] != right["rec"]["strand"]:
                continue
            right_cds = merge_intervals(
                [(start, end) for start, end, _phase, _strand in right["rec"]["cds"]]
            )
            if intersect_len(left_cds, right_cds) >= 200:
                overlap_reject.add(left_id)
                overlap_reject.add(right_id)
    for gid in overlap_reject:
        edit = applied.pop(gid)
        for member in edit["members"]:
            decisions[member] = (
                "保持原模型", member, edit["rec"]["id"], "延长后的编码区相互重叠",
            )

    for row in decisions.values():
        reason_counts[row[0] if row[0] != "保持原模型" else row[3]] += 1

    n_extend = sum(1 for row in applied.values() if row["kind"] == "extend")
    n_merge = sum(1 for row in applied.values() if row["kind"] == "merge")
    n_removed = sum(len(row["members"]) - 1 for row in applied.values() if row["kind"] == "merge")
    print(
        f"{args.prefix}\tshort={len(short)}\tclean={len(clean)}\t"
        f"extend={n_extend}\tmerge_groups={n_merge}\tmerge_removed={n_removed}\t"
        f"kept={sum(1 for row in decisions.values() if row[0] == '保持原模型')}",
        flush=True,
    )
    for key, count in sorted(reason_counts.items(), key=lambda item: -item[1]):
        print(f"[INFO] {args.prefix}\t{key}\t{count}", flush=True)
    if args.dry_run:
        return

    removed_genes = {
        gid for gid, row in decisions.items() if row[0] == "合并删除"
    }
    drop_tx = set()
    replace_tx = {}
    id_rows = []
    for gid, edit in applied.items():
        keep_tx = edit["tx"]
        for member in edit["members"]:
            for tid in tx_of.get(member, []):
                if tid != keep_tx:
                    drop_tx.add(tid)
                    id_rows.append((tid, keep_tx, "transcript", "转录本被修订模型替换"))
        replace_tx[keep_tx] = edit
        if edit["kind"] == "extend":
            id_rows.append((gid, gid, "gene", "延长"))
        else:
            id_rows.append((gid, gid, "gene", "合并保留"))
            for member in edit["members"]:
                if member != gid:
                    id_rows.append((member, gid, "gene", "合并删除"))

    new_blocks = []
    for gid, lines in blocks:
        if gid in removed_genes:
            continue
        if gid in applied:
            new_blocks.append((gid, render_gene(gid, applied[gid]["tx"], applied[gid]["rec"])))
        else:
            new_blocks.append((gid, lines))

    final_genes = [gid for gid, _lines in new_blocks]
    if len(final_genes) != len(blocks) - len(removed_genes):
        fail("gene count after revision does not match the removal list")
    if banned & set(final_genes):
        fail("a previously dropped gene entered the revised set")

    rendered = dict(new_blocks)
    for gid, edit in applied.items():
        check_parts = []
        for line in rendered[gid]:
            parts = line.rstrip("\n").split("\t")
            if parts[2] == "CDS":
                phase = int(parts[7]) if parts[7] in {"0", "1", "2"} else 0
                check_parts.append((int(parts[3]), int(parts[4]), phase, parts[6]))
        if splice(check_parts, genome[edit["rec"]["chrom"]]) != edit["cds"]:
            fail(f"rewritten GFF CDS does not reproduce the sequence: {gid}")
        if classify_cds(edit["cds"]):
            fail(f"applied CDS failed ORF rules: {gid}")

    coord = {}
    for gid, edit in applied.items():
        gene_line = rendered[gid][0].rstrip("\n").split("\t")
        coord[gid] = (int(gene_line[3]), int(gene_line[4]))

    temps = []

    def tmp(path):
        name = path + ".revise.tmp"
        temps.append(name)
        return name

    try:
        gff_tmp = tmp(gff_path)
        with open(gff_tmp, "w") as handle:
            handle.writelines(header)
            for _gid, lines in new_blocks:
                handle.writelines(lines)
        for kind, filename in (
            ("cds", f"{args.prefix}.final.cds.fa"),
            ("tx", f"{args.prefix}.final.transcripts.fa"),
            ("pep", f"{args.prefix}.final.proteins.fa"),
        ):
            src = os.path.join(final_dir, filename)
            replace = {}
            for edit in applied.values():
                if kind == "pep":
                    replace[edit["tx"]] = edit["protein"]
                else:
                    replace[edit["tx"]] = edit["cds"]
            stream_fasta(src, tmp(src), drop_tx, replace)
        list_tmp = tmp(os.path.join(final_dir, f"{args.prefix}.final.genes.list"))
        with open(list_tmp, "w") as handle:
            for gid in final_genes:
                handle.write(gid + "\n")
        class_tmp = tmp(class_path)
        with open(class_path) as inp, open(class_tmp, "w") as out:
            header_line = inp.readline()
            out.write(header_line)
            columns = header_line.rstrip("\n").split("\t")
            for row in csv.DictReader(inp, fieldnames=columns, delimiter="\t"):
                gid = row["gene_id"]
                if gid in removed_genes:
                    continue
                if gid in coord:
                    row["start"] = str(coord[gid][0])
                    row["end"] = str(coord[gid][1])
                out.write("\t".join(row[col] for col in columns) + "\n")
        corr_path = os.path.join(final_dir, f"{args.prefix}.final.corrections.tsv")
        map_path = os.path.join(final_dir, f"{args.prefix}.final.id_map.tsv")
        corr_tmp = tmp(corr_path)
        with open(corr_tmp, "w") as handle:
            handle.write(
                "gene_id\taction\tkept_gene_id\tminiprot_id\tzh13_id\t"
                "old_start\told_end\tnew_start\tnew_end\told_aa\tnew_aa\t"
                "old_protein\tnew_protein\treason\n"
            )
            for gid in sorted(decisions, key=locus_key):
                action, kept, mid, reason = decisions[gid]
                gene = short[gid]
                new_start = new_end = new_aa = ""
                old_protein = new_protein = ""
                if gid in applied:
                    new_start, new_end = coord[gid]
                    new_aa = len(applied[gid]["protein"])
                    old_tx = sorted(tx_of[gid], key=tx_key)[0]
                    old_protein = proteins.get(old_tx, "")
                    new_protein = applied[gid]["protein"]
                elif action == "合并删除":
                    new_start, new_end = coord[kept]
                    new_aa = len(applied[kept]["protein"])
                    old_tx = sorted(tx_of[gid], key=tx_key)[0]
                    old_protein = proteins.get(old_tx, "")
                    new_protein = applied[kept]["protein"]
                handle.write(
                    f"{gid}\t{action}\t{kept}\t{mid}\t{gene['zh13_id']}\t"
                    f"{gene['start']}\t{gene['end']}\t{new_start}\t{new_end}\t"
                    f"{old_aa.get(gid, '')}\t{new_aa}\t{old_protein}\t{new_protein}\t{reason}\n"
                )
        map_tmp = tmp(map_path)
        with open(map_tmp, "w") as handle:
            handle.write("old_id\tnew_id\tid_type\taction\n")
            for old, new, kind, action in id_rows:
                handle.write(f"{old}\t{new}\t{kind}\t{action}\n")
        with open(criteria_path) as handle:
            criteria = handle.read()
        if STAMP in criteria:
            criteria = criteria[:criteria.index(STAMP)].rstrip() + "\n"
        criteria += (
            f"\n{STAMP}\n"
            "只修订 high_confidence 里、相对 ZH13 v3 查询覆盖 >= 0.5 且参考覆盖 < 0.5 的短模型。\n"
            "expressed_zh13_below_50 不修订。这类基因相对 ZH13 v3 的相似度低于 50，"
            "不能用 ZH13 全长蛋白改写。homolog_unexpressed 也不修订。\n"
            "同位点 miniprot 覆盖 ZH13 蛋白至少 80%，无移码、无框内终止、相似度 >= 0.8，"
            "并且不与另一个 ZH13 基因的编码区重叠达到 200 bp，才提出修订。\n"
            "同一 ZH13 基因只有一个这样的模型时，在原基因号上延长。"
            "有多个模型时，只有落在同一条 miniprot 比对上的才合并，保留编号最小的基因号；"
            "不在同一条比对上的保持原模型。\n"
            "新 CDS 必须同时满足：有 ATG、有终止密码子、长度能被 3 整除、没有内部终止、不含非 ACGT。"
            "新蛋白没有变长，或新编码区与仍保留的其他基因重叠达到 200 bp，则不改该基因。\n"
            "同一次修订中，两个新 CDS 在同链上重叠达到 200 bp 时，这两个延长都不写入。\n"
            "已在 12_renamed 删除清单中的基因不会回到注释里。位点编号不重排。\n"
            f"本轮延长 {n_extend} 个基因，合并保留 {n_merge} 个基因，"
            f"合并删除 {n_removed} 个基因。修订后基因数 {len(final_genes)}。\n"
            "覆盖前后的基因号、坐标和蛋白序列记在 final.corrections.tsv；"
            "基因号或转录本号发生变化时记在 final.id_map.tsv。\n"
        )
        crit_tmp = tmp(criteria_path)
        write_text(crit_tmp, criteria)

        verify_outputs(
            args, gff_tmp, final_genes, class_tmp, list_tmp,
            {edit["tx"]: edit for edit in applied.values()}, drop_tx, temps, final_dir,
        )
        replacements = {
            gff_tmp: gff_path,
            list_tmp: os.path.join(final_dir, f"{args.prefix}.final.genes.list"),
            class_tmp: class_path,
            corr_tmp: corr_path,
            map_tmp: map_path,
            crit_tmp: criteria_path,
        }
        for kind, filename in (
            ("cds", f"{args.prefix}.final.cds.fa"),
            ("tx", f"{args.prefix}.final.transcripts.fa"),
            ("pep", f"{args.prefix}.final.proteins.fa"),
        ):
            src = os.path.join(final_dir, filename)
            replacements[src + ".revise.tmp"] = src
        for src, dst in replacements.items():
            os.replace(src, dst)
            temps.remove(src)
    finally:
        for path in temps:
            if os.path.exists(path):
                os.remove(path)
    print(
        f"{args.prefix}\tapplied\textend={n_extend}\tmerge_removed={n_removed}\t"
        f"genes={len(final_genes)}",
        flush=True,
    )


def verify_outputs(args, gff_tmp, final_genes, class_tmp, list_tmp, edits_by_tx, drop_tx, _temps, final_dir):
    genes = []
    tx = []
    with open(gff_tmp) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            if parts[2] == "gene":
                gid = attr(parts[8], "ID")
                genes.append(gid)
                gstart, gend = int(parts[3]), int(parts[4])
                # child bounds are checked after the block; stored on the fly below
                parts = parts
            elif parts[2] == "mRNA":
                tid = attr(parts[8], "ID")
                tx.append(tid)
                if int(parts[3]) < gstart or int(parts[4]) > gend:
                    fail(f"mRNA outside gene: {tid}")
            elif parts[2] in {"CDS", "exon"}:
                if int(parts[3]) < gstart or int(parts[4]) > gend:
                    fail(f"feature outside gene near {tx[-1] if tx else genes[-1]}")
    if genes != final_genes:
        fail("rewritten GFF gene order does not match the planned list")
    with open(list_tmp) as handle:
        listed = [line.strip() for line in handle if line.strip()]
    if listed != genes:
        fail("gene list does not match GFF")
    with open(class_tmp) as handle:
        class_genes = [row["gene_id"] for row in csv.DictReader(handle, delimiter="\t")]
    if class_genes != genes:
        fail("gene class table does not match GFF")
    for filename in (
        f"{args.prefix}.final.cds.fa",
        f"{args.prefix}.final.transcripts.fa",
        f"{args.prefix}.final.proteins.fa",
    ):
        ids = []
        seqs = {}
        name = None
        chunks = []
        with open(os.path.join(final_dir, filename) + ".revise.tmp") as handle:
            for line in handle:
                if line.startswith(">"):
                    if name:
                        seqs[name] = "".join(chunks)
                    name = line[1:].split()[0]
                    ids.append(name)
                    chunks = []
                else:
                    chunks.append(line.strip())
            if name:
                seqs[name] = "".join(chunks)
        if set(ids) != set(tx) or len(ids) != len(tx):
            missing = sorted(set(tx) - set(ids))[:3]
            extra = sorted(set(ids) - set(tx))[:3]
            fail(f"{filename} IDs do not match GFF transcripts, missing {missing}, extra {extra}")
        if set(ids) & drop_tx:
            fail(f"{filename} still contains a removed transcript")
        for tid, edit in edits_by_tx.items():
            if filename.endswith(".proteins.fa"):
                if seqs[tid] != edit["protein"]:
                    fail(f"protein sequence was not overwritten: {tid}")
            elif seqs[tid] != edit["cds"]:
                fail(f"nucleotide sequence was not overwritten: {tid}")
            if filename.endswith(".cds.fa") and classify_cds(seqs[tid]):
                fail(f"overwritten CDS is not a complete ORF: {tid}")


if __name__ == "__main__":
    main()
