#!/usr/bin/env python3
"""Remove unexpressed partial copies from one 13_final gene set.

A high_confidence gene is a candidate when it has no expression support, no
ZH13 v3 or Wm82 a6 protein with identity >= 50 is covered by half or more, and
the same ZH13 v3 gene already has a full-length copy elsewhere in this genome.
Positions are converted to ZH13 v3 coordinates through the nearest anchor genes
on each side, so a segment moved by a translocation is compared at its ZH13
locus. A candidate stays when it maps to the chromosome of that ZH13 gene and
is closer to the ZH13 position than every full-length copy. All candidates are
listed in PREFIX.final.short_fragment_qc.tsv; removed IDs also go to
PREFIX.final.id_map.tsv.
"""

import argparse
import bisect
import os
import subprocess
from collections import Counter, defaultdict

from antisense_qc import ZH13_GFF, attr, gff_blocks, load_fasta, read_table, transcript_ids, write_table

DIAMOND = "/home/zzh/anaconda3/envs/annp/bin/diamond"
REFS = (
    ("ZH13v3", "/tiandata2/zzh/tools/ref_diamond/ZH13v3.0_mRNA_nuclear_20250729_pro.dmnd"),
    ("Wm82a6", "/tiandata2/zzh/tools/ref_diamond/glyma.Wm82.gnm6.ann1.PKSW.protein.dmnd"),
)
STAMP = "短片段筛选已写入"
TMP = ".fragment.tmp"


def zh13_gene(hit):
    return hit.rsplit(".", 1)[0] if hit.startswith("SoyZH13") else hit


def union_len(intervals):
    ordered = sorted(intervals)
    total = 0
    start, end = ordered[0]
    for left, right in ordered[1:]:
        if left > end + 1:
            total += end - start + 1
            start, end = left, right
        else:
            end = max(end, right)
    return total + end - start + 1


def load_zh13_positions(path):
    positions = {}
    with open(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 9 and parts[2] == "gene":
                positions[attr(parts[8])["ID"]] = (parts[0], int(parts[3]))
    return positions


def synteny_anchors(zh13_hits, classes, zh13_pos, min_pident, min_cov):
    """Per chromosome, genes whose top ZH13 hit is near-identical, full length and hit by no other gene."""
    best = {}
    for gene_id, rows in zh13_hits.items():
        _bits, subject, pident, scov = rows[0]
        if gene_id in classes and pident >= min_pident and scov >= min_cov:
            best[gene_id] = zh13_gene(subject)
    uses = Counter(best.values())
    anchors = defaultdict(list)
    for gene_id, locus in best.items():
        if uses[locus] == 1 and locus in zh13_pos:
            anchors[classes[gene_id]["chrom"]].append((int(classes[gene_id]["start"]), *zh13_pos[locus], locus))
    for rows in anchors.values():
        rows.sort()
    return anchors


def project(anchors, chrom, pos, skip):
    """ZH13 chromosome and position of one point, ignoring the anchor of locus `skip`."""
    rows = anchors.get(chrom, [])
    i = bisect.bisect_left(rows, (pos,))
    left = next((rows[j] for j in range(i - 1, -1, -1) if rows[j][3] != skip), None)
    right = next((rows[j] for j in range(i, len(rows)) if rows[j][3] != skip), None)
    if left and right and left[1] == right[1] and right[0] > left[0]:
        return left[1], round(left[2] + (pos - left[0]) / (right[0] - left[0]) * (right[2] - left[2]))
    flank = [row for row in (left, right) if row]
    if not flank:
        return None
    near = min(flank, key=lambda row: abs(row[0] - pos))
    return near[1], near[2]


def write_longest(final, prefix, dest):
    seqs = load_fasta(os.path.join(final, f"{prefix}.final.proteins.fa"))
    best = {}
    for transcript_id, seq in seqs.items():
        gene_id = transcript_id.rsplit(".m", 1)[0]
        seq = seq.rstrip("*")
        if gene_id not in best or len(seq) > len(best[gene_id]):
            best[gene_id] = seq
    with open(dest, "w") as handle:
        for gene_id, seq in best.items():
            handle.write(f">{gene_id}\n{seq}\n")
    return {gene_id: len(seq) for gene_id, seq in best.items()}


def run_diamond(diamond, query, db, out, threads):
    subprocess.run(
        [
            diamond, "blastp", "-q", query, "-d", db, "-o", out,
            "-f", "6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen",
            "--evalue", "1e-5", "--max-target-seqs", "5", "--threads", str(threads), "--quiet",
        ],
        check=True,
    )


def load_hits(path):
    """Per query, one row per subject: HSPs merged on the subject, ordered by summed bitscore."""
    raw = defaultdict(dict)
    with open(path) as handle:
        for line in handle:
            p = line.rstrip("\n").split("\t")
            bits = float(p[11])
            rec = raw[p[0]].get(p[1])
            if rec is None:
                rec = raw[p[0]][p[1]] = {"bits": 0.0, "top": -1.0, "pident": 0.0,
                                         "s": [], "slen": int(p[13])}
            rec["bits"] += bits
            if bits > rec["top"]:
                rec["top"], rec["pident"] = bits, float(p[2])
            s1, s2 = int(p[8]), int(p[9])
            rec["s"].append((min(s1, s2), max(s1, s2)))
    out = {}
    for query, subjects in raw.items():
        rows = [
            (rec["bits"], subject, rec["pident"], union_len(rec["s"]) / rec["slen"])
            for subject, rec in subjects.items()
        ]
        rows.sort(reverse=True)
        out[query] = rows
    return out


def write_fasta(path, records):
    temporary = path + TMP
    with open(temporary, "w") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")
            for pos in range(0, len(seq), 70):
                handle.write(seq[pos:pos + 70] + "\n")
    return temporary


def apply_exclusion(final, prefix, removed, classes):
    for gene_id in removed:
        row = classes[gene_id]
        if row["class"] != "high_confidence" or row["expression_supported"] == "yes":
            raise SystemExit(f"{prefix}: refusing to remove {gene_id} ({row['class']}, "
                             f"expression {row['expression_supported']})")
    gff_path = os.path.join(final, f"{prefix}.final.gff3")
    header, blocks = gff_blocks(gff_path)
    transcripts = {gene_id: transcript_ids(lines) for gene_id, lines in blocks}
    missing = sorted(set(removed) - set(transcripts))
    if missing:
        raise SystemExit(f"{prefix}: exclusion ids missing from GFF {missing[:5]}")

    gff_tmp = gff_path + TMP
    with open(gff_tmp, "w") as handle:
        handle.writelines(header)
        for gene_id, lines in blocks:
            if gene_id not in removed:
                handle.writelines(lines)

    fasta_tmp = {}
    for kind in ("transcripts", "cds", "proteins"):
        path = os.path.join(final, f"{prefix}.final.{kind}.fa")
        seqs = load_fasta(path)
        records = []
        for gene_id, _lines in blocks:
            if gene_id in removed:
                continue
            for transcript_id in transcripts[gene_id]:
                if transcript_id not in seqs:
                    raise SystemExit(f"{prefix}: missing sequence {transcript_id}")
                records.append((transcript_id, seqs[transcript_id]))
        fasta_tmp[path] = write_fasta(path, records)

    list_path = os.path.join(final, f"{prefix}.final.genes.list")
    with open(list_path) as handle:
        gene_list = [line.strip() for line in handle if line.strip() and line.strip() not in removed]
    with open(list_path + TMP, "w") as handle:
        handle.write("\n".join(gene_list) + "\n")

    class_path = os.path.join(final, f"{prefix}.final.gene_class.tsv")
    with open(class_path) as handle:
        class_lines = handle.readlines()
    class_header = class_lines[0].rstrip("\n").split("\t")
    gene_at = class_header.index("gene_id")
    class_at = class_header.index("class")
    counts = defaultdict(int)
    with open(class_path + TMP, "w") as handle:
        handle.write(class_lines[0])
        for line in class_lines[1:]:
            parts = line.rstrip("\n").split("\t")
            if parts[gene_at] in removed:
                continue
            counts[parts[class_at]] += 1
            handle.write(line)

    map_path = os.path.join(final, f"{prefix}.final.id_map.tsv")
    with open(map_path) as handle:
        map_text = handle.read()
    if not map_text.endswith("\n"):
        map_text += "\n"
    with open(map_path + TMP, "w") as handle:
        handle.write(map_text)
        for gene_id in sorted(removed):
            kept = removed[gene_id]
            handle.write(f"{gene_id}\t{kept}\tgene\t短片段排除\n")
            for transcript_id in transcripts[gene_id]:
                handle.write(f"{transcript_id}\t{kept}\ttranscript\t短片段排除\n")

    written = []
    gff_tx = []
    with open(gff_tmp) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 9 and parts[2] == "gene":
                written.append(attr(parts[8]).get("ID", ""))
            elif len(parts) >= 9 and parts[2] == "mRNA":
                gff_tx.append(attr(parts[8]).get("ID", ""))
    if written != gene_list:
        raise SystemExit(f"{prefix}: GFF gene order does not match genes.list")
    for path, temporary in fasta_tmp.items():
        if set(load_fasta(temporary)) != set(gff_tx):
            raise SystemExit(f"{prefix}: {os.path.basename(path)} IDs do not match GFF")

    os.replace(gff_tmp, gff_path)
    for path, temporary in fasta_tmp.items():
        os.replace(temporary, path)
    os.replace(list_path + TMP, list_path)
    os.replace(class_path + TMP, class_path)
    os.replace(map_path + TMP, map_path)
    return counts, len(written)


def run_sample(downstream, prefix, args):
    final = os.path.join(downstream, "13_final")
    criteria_path = os.path.join(final, f"{prefix}.final.criteria.txt")
    with open(criteria_path) as handle:
        criteria = handle.read()
    if STAMP in criteria:
        print(f"{prefix}\talready_filtered", flush=True)
        return
    work = os.path.join(downstream, "05_diamond")
    query = os.path.join(work, f"{prefix}.final_longest.faa")
    aa_len = write_longest(final, prefix, query)
    hits = {}
    for label, db in REFS:
        out = os.path.join(work, f"{prefix}.final_longest.vs_{label}.tsv")
        run_diamond(args.diamond, query, db, out, args.threads)
        hits[label] = load_hits(out)

    full_copies = defaultdict(dict)
    for gene_id, rows in hits["ZH13v3"].items():
        for _bits, subject, pident, scov in rows:
            locus = zh13_gene(subject)
            if pident >= args.min_pident and scov >= args.full_cov:
                full_copies[locus][gene_id] = max(scov, full_copies[locus].get(gene_id, 0.0))

    zh13_pos = load_zh13_positions(args.zh13_gff)
    classes = read_table(os.path.join(final, f"{prefix}.final.gene_class.tsv"))
    anchors = synteny_anchors(hits["ZH13v3"], classes, zh13_pos, args.anchor_pident, args.anchor_cov)

    def mapped(gene_id, locus):
        return project(anchors, classes[gene_id]["chrom"], int(classes[gene_id]["start"]), locus)

    def distance(point, zh13_chrom, zh13_start):
        if point is None or point[0] != zh13_chrom:
            return None
        return abs(point[1] - zh13_start)

    summary = read_table(os.path.join(downstream, "12_renamed", f"{prefix}.renamed.annotation.summary.tsv"))
    repair = read_table(os.path.join(final, f"{prefix}.final.corrections.tsv"))
    _header, blocks = gff_blocks(os.path.join(final, f"{prefix}.final.gff3"))
    exons = {}
    for gene_id, lines in blocks:
        per_tx = defaultdict(int)
        for line in lines:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 9 and parts[2] == "CDS":
                per_tx[attr(parts[8])["Parent"]] += 1
        exons[gene_id] = max(per_tx.values()) if per_tx else 0
    antisense = defaultdict(list)
    pair_path = os.path.join(final, f"{prefix}.final.antisense_pair_qc.tsv")
    if os.path.isfile(pair_path):
        for row in read_table(pair_path).values():
            if row["evidence_tier"] in ("较高", "中等至较高"):
                antisense[row["gene_a"]].append(row["gene_b"])
                antisense[row["gene_b"]].append(row["gene_a"])

    candidates = []
    removed = {}
    for gene_id, row in classes.items():
        if row["class"] != "high_confidence" or row["expression_supported"] == "yes":
            continue
        best = ("NA", "NA", 0.0, 0.0)
        for label, _db in REFS:
            for _bits, subject, pident, scov in hits[label].get(gene_id, []):
                if pident >= args.min_pident and scov > best[3]:
                    best = (label, subject, pident, scov)
        if best[3] >= args.max_cov:
            continue
        top = hits["ZH13v3"].get(gene_id)
        if not top:
            continue
        locus = zh13_gene(top[0][1])
        copies = full_copies.get(locus)
        if not copies:
            continue
        zh13_chrom, zh13_start = zh13_pos[locus]
        point = mapped(gene_id, locus)
        own = distance(point, zh13_chrom, zh13_start)
        near = []
        for copy in copies:
            copy_distance = distance(mapped(copy, locus), zh13_chrom, zh13_start)
            if copy_distance is not None:
                near.append((copy_distance, copy))
        nearest = min(near) if near else None
        if nearest is not None:
            kept_gene = nearest[1]
        else:
            kept_gene = max(copies, key=lambda copy: (copies[copy], copy))
        if own is not None and (nearest is None or own < nearest[0]):
            action, reason = "保留", "同位点短模型：换算到ZH13坐标后位于该基因所在染色体，且比所有完整拷贝更靠近ZH13位置"
        else:
            action, reason = "排除", "已有完整拷贝的短片段"
            removed[gene_id] = kept_gene
        candidates.append({
            "gene_id": gene_id,
            "action": action,
            "reason": reason,
            "full_copy_gene": kept_gene,
            "chrom": row["chrom"],
            "start": row["start"],
            "end": row["end"],
            "strand": row["strand"],
            "protein_aa": aa_len[gene_id],
            "cds_exons": exons.get(gene_id, 0),
            "zh13_gene": locus,
            "zh13_position": f"{zh13_chrom}:{zh13_start}",
            "mapped_zh13_position": "NA" if point is None else f"{point[0]}:{point[1]}",
            "distance_to_zh13": "NA" if own is None else own,
            "full_copy_distance_to_zh13": "NA" if nearest is None else nearest[0],
            "zh13_top_hit": top[0][1],
            "zh13_top_pident": f"{top[0][2]:.1f}",
            "zh13_top_scov": f"{top[0][3]:.4f}",
            "best_ref": best[0],
            "best_ref_hit": best[1],
            "best_ref_pident": f"{best[2]:.1f}",
            "best_ref_scov": f"{best[3]:.4f}",
            "full_copy_scov": f"{copies[kept_gene]:.4f}",
            "max_tpm": summary.get(gene_id, {}).get("max_TPM", "NA"),
            "short_model_repair": repair.get(gene_id, {}).get("reason", "NA"),
            "kept_antisense_partner": ",".join(antisense.get(gene_id, [])) or "NA",
        })
    header = [
        "gene_id", "action", "reason", "full_copy_gene", "chrom", "start", "end", "strand",
        "protein_aa", "cds_exons", "zh13_gene", "zh13_position", "mapped_zh13_position", "distance_to_zh13",
        "full_copy_distance_to_zh13", "zh13_top_hit", "zh13_top_pident", "zh13_top_scov",
        "best_ref", "best_ref_hit", "best_ref_pident", "best_ref_scov", "full_copy_scov",
        "max_tpm", "short_model_repair", "kept_antisense_partner",
    ]
    candidates.sort(key=lambda item: (item["chrom"], int(item["start"]), item["gene_id"]))
    write_table(os.path.join(final, f"{prefix}.final.short_fragment_qc.tsv"), header, candidates)
    class_counts, n_genes = apply_exclusion(final, prefix, removed, classes)
    n_kept = len(candidates) - len(removed)
    pairs_hit = {
        tuple(sorted((gene_id, partner)))
        for gene_id in removed for partner in antisense.get(gene_id, [])
    }
    kept_pairs = sum(len(v) for v in antisense.values()) // 2
    note = (
        f"{STAMP}\n"
        "短片段筛选没有另建一套注释，位点编号不重排。\n"
        f"候选同时满足：high_confidence；没有表达支持；最终最长蛋白与 ZH13 v3、Wm82 a6 中相似度 >= {args.min_pident:g} 的蛋白比对时，"
        f"同一参考蛋白上合并后的覆盖都低于 {args.max_cov:g}；"
        f"最高 bitscore 命中的 ZH13 v3 基因在本基因组已有覆盖 >= {args.full_cov:g} 的完整拷贝。\n"
        "候选和完整拷贝的位置先按两侧最近的锚定基因换算到 ZH13 v3 坐标（锚定基因的最高命中相似度 >= "
        f"{args.anchor_pident:g}、参考蛋白覆盖 >= {args.anchor_cov:g}，且该 ZH13 基因只被它一个命中），"
        "易位片段因此按其在 ZH13 中的位点比较。"
        "候选换算后落在该 ZH13 基因所在染色体，且比所有完整拷贝更靠近 ZH13 基因位置时，视为同位点模型保留；其余候选作为已有完整拷贝的基因片段或假基因排除。"
        "ZH13 基因在本基因组只有短拷贝的模型不是候选，保留。\n"
        f"比对文件在 05_diamond/{prefix}.final_longest.vs_*.tsv。全部候选及处理结果在 final.short_fragment_qc.tsv，排除的基因同时记在 final.id_map.tsv。\n"
        f"候选 {len(candidates)} 个，同位点保留 {n_kept} 个。\n"
        f"其中 {len(pairs_hit)} 对原保留的反义对有一侧被排除，保留的反义对由 {kept_pairs} 对减为 {kept_pairs - len(pairs_hit)} 对。\n"
        f"本轮排除 {len(removed)} 个基因。处理后基因数 {n_genes}，其中 high_confidence {class_counts['high_confidence']}，"
        f"expressed_zh13_below_50 {class_counts['expressed_zh13_below_50']}，"
        f"homolog_unexpressed {class_counts['homolog_unexpressed']}。\n"
    )
    if not criteria.endswith("\n"):
        criteria += "\n"
    with open(criteria_path, "w") as handle:
        handle.write(criteria + "\n" + note)
    print(
        f"{prefix}\tcandidates={len(candidates)}\tkept={n_kept}\tremoved={len(removed)}\tantisense_pairs_hit={len(pairs_hit)}\tgenes={n_genes}\t"
        f"hc={class_counts['high_confidence']}\tbelow={class_counts['expressed_zh13_below_50']}\t"
        f"homolog={class_counts['homolog_unexpressed']}",
        flush=True,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--downstream")
    parser.add_argument("--prefix")
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--diamond", default=DIAMOND)
    parser.add_argument("--zh13-gff", default=ZH13_GFF)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--min-pident", type=float, default=50)
    parser.add_argument("--max-cov", type=float, default=0.5)
    parser.add_argument("--full-cov", type=float, default=0.8)
    parser.add_argument("--anchor-pident", type=float, default=98)
    parser.add_argument("--anchor-cov", type=float, default=0.9)
    args = parser.parse_args()
    if args.all:
        root = "/tiandata2/zzh/sjy_assembly"
        samples = [
            ("KD35", f"{root}/annot/KD35_annotation_downstream"),
            ("N2L", f"{root}/annot/N2L_annotation_downstream"),
            ("BZ", f"{root}/batch2/annot/BZ_annotation_downstream"),
            ("DY11", f"{root}/batch2/annot/DY11_annotation_downstream"),
            ("TH2", f"{root}/batch2/annot/TH2_annotation_downstream"),
            ("TZX1998", f"{root}/batch2/annot/TZX1998_annotation_downstream"),
            ("TZX253", f"{root}/batch2/annot/TZX253_annotation_downstream"),
        ]
        for prefix, downstream in samples:
            run_sample(downstream, prefix, args)
        return
    if not args.downstream or not args.prefix:
        parser.error("give --downstream and --prefix, or --all")
    run_sample(args.downstream, args.prefix, args)


if __name__ == "__main__":
    main()
