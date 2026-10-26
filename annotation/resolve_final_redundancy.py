#!/usr/bin/env python3
"""Clean one 13_final directory in place.

Three edits, all on the same strand:

- Gene IDs that share one CDS sequence and one protein are merged. The ID with
  more UTR exons is kept. The other ID is removed, not added as a second copy
  of the same transcript.
- An extension is restored from 12_renamed when two new CDS overlap by at least
  200 bp and the original CDS overlapped by less than 200 bp.
- An unexpressed gene is removed when its CDS is fully inside a longer gene of
  the same ZH13 locus and every short protein is a continuous piece of a longer
  protein.

Opposite-strand overlaps are not removed. Locus IDs are not renumbered.
"""

import argparse
import os
import re
from collections import defaultdict

STAMP = "冗余合并与碎片排除已写入"


def attr(field):
    out = {}
    for part in field.split(";"):
        if "=" in part:
            key, value = part.split("=", 1)
            out[key] = value
    return out


def locus(hit):
    if not hit or hit == "NA":
        return "NA"
    return re.sub(r"\.m\d+$", "", hit)


def locus_key(gene_id):
    match = re.search(r"_(\d+)G(\d+)$", gene_id)
    if not match:
        return (999, 10**12)
    return (int(match.group(1)), int(match.group(2)))


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


def interval_len(intervals):
    return sum(end - start + 1 for start, end in merge_intervals(intervals))


def overlap_len(left, right):
    left = merge_intervals(left)
    right = merge_intervals(right)
    total = 0
    for start, end in left:
        for other_start, other_end in right:
            if other_end < start:
                continue
            if other_start > end:
                break
            total += min(end, other_end) - max(start, other_start) + 1
    return total


def fully_contained(short, long):
    short_len = interval_len(short)
    return short_len > 0 and short_len < interval_len(long) and overlap_len(short, long) == short_len


def load_table(path, columns):
    rows = {}
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        index = {name: pos for pos, name in enumerate(header)}
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            gene_id = parts[index["gene_id"]]
            rows[gene_id] = {name: parts[index[name]] if name in index else "" for name in columns}
            rows[gene_id]["_line"] = line
    return header, rows


def load_fasta(path):
    seqs = {}
    name = None
    buf = []

    def flush():
        if name:
            seqs[name] = "".join(buf)

    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                flush()
                name = line[1:].split()[0]
                buf = []
            elif name:
                buf.append(line.strip())
        flush()
    return seqs


def gene_of(transcript_id):
    return transcript_id.rsplit(".", 1)[0]


def load_gff_blocks(path):
    header = []
    blocks = []
    current_id = None
    current = None
    seen = False
    with open(path) as handle:
        for line in handle:
            if not seen and (line.startswith("#") or not line.strip()):
                header.append(line)
                continue
            parts = line.rstrip("\n").split("\t")
            is_gene = len(parts) >= 9 and parts[2] == "gene"
            if is_gene:
                seen = True
                if current is not None:
                    blocks.append((current_id, current))
                current_id = attr(parts[8]).get("ID", "")
                current = [line]
            elif current is not None:
                current.append(line)
            else:
                header.append(line)
    if current is not None:
        blocks.append((current_id, current))
    return header, blocks


def block_models(lines):
    transcripts = {}
    order = []
    gene = None
    for line in lines:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue
        info = attr(parts[8])
        start, end = int(parts[3]), int(parts[4])
        if parts[2] == "gene":
            gene = {
                "id": info.get("ID", ""),
                "chrom": parts[0],
                "start": start,
                "end": end,
                "strand": parts[6],
            }
        elif parts[2] == "mRNA":
            transcript_id = info.get("ID", "")
            transcripts[transcript_id] = {"exon": [], "cds": [], "utr": 0}
            order.append(transcript_id)
        else:
            parent = info.get("Parent", "").split(",")[0]
            if parent not in transcripts:
                continue
            if parts[2] == "exon":
                transcripts[parent]["exon"].append((start, end))
            elif parts[2] == "CDS":
                transcripts[parent]["cds"].append((start, end))
            elif parts[2] in ("five_prime_UTR", "three_prime_UTR"):
                transcripts[parent]["utr"] += 1
    return gene, transcripts, order


def write_fasta(path, records):
    temporary = path + ".resolve.tmp"
    with open(temporary, "w") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")
            for pos in range(0, len(seq), 70):
                handle.write(seq[pos:pos + 70] + "\n")
    return temporary


def proteins_are_fragment(short_proteins, long_proteins):
    if not short_proteins or not long_proteins:
        return False
    if min(map(len, short_proteins)) >= max(map(len, long_proteins)):
        return False
    return all(
        any(protein and protein in longer and len(protein) < len(longer) for longer in long_proteins)
        for protein in short_proteins
    )


def class_counts(path, removed):
    counts = defaultdict(int)
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        class_at = header.index("class")
        gene_at = header.index("gene_id")
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if parts[gene_at] in removed:
                continue
            counts[parts[class_at]] += 1
    return counts


def resolve(downstream, prefix, dry_run):
    final_dir = os.path.join(downstream, "13_final")
    renamed = os.path.join(downstream, "12_renamed")
    gff_path = os.path.join(final_dir, f"{prefix}.final.gff3")
    criteria_path = os.path.join(final_dir, f"{prefix}.final.criteria.txt")
    with open(criteria_path) as handle:
        criteria = handle.read()
    if STAMP in criteria and not dry_run:
        print(f"{prefix}\talready_resolved", flush=True)
        return

    header, blocks = load_gff_blocks(gff_path)
    models = {}
    for gene_id, lines in blocks:
        gene, transcripts, order = block_models(lines)
        if gene:
            models[gene_id] = {"gene": gene, "tx": transcripts, "order": order, "lines": lines}

    cds_fa = load_fasta(os.path.join(final_dir, f"{prefix}.final.cds.fa"))
    prot_fa = load_fasta(os.path.join(final_dir, f"{prefix}.final.proteins.fa"))
    tx_fa = load_fasta(os.path.join(final_dir, f"{prefix}.final.transcripts.fa"))
    summary_cols = [
        "ZH13v3.0_mRNA_nuclear_20250729_pro_hit",
        "expression_supported",
        "max_TPM",
    ]
    _summary_header, summary = load_table(
        os.path.join(renamed, f"{prefix}.renamed.annotation.summary.tsv"),
        summary_cols,
    )
    _corr_header, corrections = load_table(
        os.path.join(final_dir, f"{prefix}.final.corrections.tsv"),
        ["action", "zh13_id", "old_start", "old_end", "new_start", "new_end", "old_aa", "new_aa"],
    )

    def sequences_for(gene_id, store):
        values = []
        for transcript_id in models[gene_id]["order"]:
            if transcript_id in store:
                values.append(store[transcript_id])
        return values

    def cds_intervals(gene_id):
        intervals = []
        for transcript in models[gene_id]["tx"].values():
            intervals.extend(transcript["cds"])
        return merge_intervals(intervals)

    def zh13_of(gene_id):
        return locus(summary.get(gene_id, {}).get("ZH13v3.0_mRNA_nuclear_20250729_pro_hit", "NA"))

    def expressed(gene_id):
        return summary.get(gene_id, {}).get("expression_supported", "no") == "yes"

    def structure_score(gene_id):
        transcripts = models[gene_id]["tx"]
        utr = sum(item["utr"] for item in transcripts.values())
        exons = max((len(item["exon"]) for item in transcripts.values()), default=0)
        chrom, number = locus_key(gene_id)
        return (utr, exons, len(transcripts), 1 if expressed(gene_id) else 0, -chrom, -number)

    buckets = defaultdict(list)
    for gene_id in models:
        cds_set = tuple(sorted(set(sequences_for(gene_id, cds_fa))))
        prot_set = tuple(sorted(set(sequences_for(gene_id, prot_fa))))
        if not cds_set or not prot_set:
            continue
        gene = models[gene_id]["gene"]
        buckets[(gene["chrom"], gene["strand"], cds_set, prot_set, zh13_of(gene_id))].append(gene_id)

    parent = {}

    def find(gene_id):
        parent.setdefault(gene_id, gene_id)
        if parent[gene_id] != gene_id:
            parent[gene_id] = find(parent[gene_id])
        return parent[gene_id]

    def union(left, right):
        left, right = find(left), find(right)
        if left != right:
            parent[right] = left

    for group in buckets.values():
        if len(group) < 2:
            continue
        for left_i, left in enumerate(group):
            for right in group[left_i + 1:]:
                if overlap_len(cds_intervals(left), cds_intervals(right)) > 0:
                    union(left, right)
    clusters = defaultdict(list)
    for group in buckets.values():
        for gene_id in group:
            if gene_id in parent:
                clusters[find(gene_id)].append(gene_id)

    drop_merge = {}
    merge_rows = []
    for members in clusters.values():
        if len(members) < 2:
            continue
        ranked = sorted(members, key=structure_score, reverse=True)
        kept = ranked[0]
        for removed in ranked[1:]:
            drop_merge[removed] = kept
            gene = models[removed]["gene"]
            merge_rows.append((
                removed, kept, gene["chrom"], gene["strand"], gene["start"], gene["end"],
                interval_len(cds_intervals(removed)),
                max(len(seq) for seq in sequences_for(removed, prot_fa)),
                zh13_of(removed),
                "yes" if expressed(removed) else "no",
                structure_score(kept)[0],
                structure_score(removed)[0],
            ))

    renamed_header, renamed_blocks = load_gff_blocks(os.path.join(renamed, f"{prefix}.renamed.gff3"))
    renamed_models = {}
    for gene_id, lines in renamed_blocks:
        gene, transcripts, order = block_models(lines)
        if gene:
            renamed_models[gene_id] = {"gene": gene, "tx": transcripts, "order": order, "lines": lines}
    renamed_cds = None
    renamed_prot = None
    renamed_tx = None

    def original_intervals(gene_id):
        intervals = []
        for transcript in renamed_models[gene_id]["tx"].values():
            intervals.extend(transcript["cds"])
        return merge_intervals(intervals)

    extended = [
        gene_id for gene_id, row in corrections.items()
        if row.get("action") == "延长" and gene_id in models and gene_id in renamed_models
    ]
    by_place = defaultdict(list)
    for gene_id in extended:
        gene = models[gene_id]["gene"]
        by_place[(gene["chrom"], gene["strand"])].append(gene_id)
    revert = set()
    conflict_overlap = {}
    for group in by_place.values():
        group.sort(key=lambda gene_id: models[gene_id]["gene"]["start"])
        for left_i, left in enumerate(group):
            for right in group[left_i + 1:]:
                if models[right]["gene"]["start"] > models[left]["gene"]["end"] + 500000:
                    break
                new_overlap = overlap_len(cds_intervals(left), cds_intervals(right))
                old_overlap = overlap_len(original_intervals(left), original_intervals(right))
                if new_overlap >= 200 and old_overlap < 200:
                    revert.add(left)
                    revert.add(right)
                    conflict_overlap[left] = max(conflict_overlap.get(left, 0), new_overlap)
                    conflict_overlap[right] = max(conflict_overlap.get(right, 0), new_overlap)
    both = set(drop_merge) & revert
    if both:
        raise SystemExit(f"{prefix}: gene is both a duplicate and a reverted extension: {sorted(both)}")

    if revert:
        renamed_cds = load_fasta(os.path.join(renamed, f"{prefix}.renamed.cds.fa"))
        renamed_prot = load_fasta(os.path.join(renamed, f"{prefix}.renamed.proteins.fa"))
        renamed_tx = load_fasta(os.path.join(renamed, f"{prefix}.renamed.transcripts.fa"))

    def effective_intervals(gene_id):
        if gene_id in revert:
            return original_intervals(gene_id)
        return cds_intervals(gene_id)

    def effective_proteins(gene_id):
        if gene_id in revert:
            return [renamed_prot[tid] for tid in renamed_models[gene_id]["order"] if tid in renamed_prot]
        return sequences_for(gene_id, prot_fa)

    remain = [gene_id for gene_id in models if gene_id not in drop_merge]
    by_strand = defaultdict(list)
    for gene_id in remain:
        gene = models[gene_id]["gene"]
        intervals = effective_intervals(gene_id)
        if not intervals:
            continue
        by_strand[(gene["chrom"], gene["strand"])].append(
            (intervals[0][0], intervals[-1][1], gene_id)
        )

    fragment_rows = []
    drop_fragment = {}
    for group in by_strand.values():
        group.sort()
        active = []
        for start, end, gene_id in group:
            active = [row for row in active if row[0] >= start]
            for other_end, other_start, other in active:
                pair = (gene_id, other)
                for short, long in (pair, (other, gene_id)):
                    if short in drop_fragment or long in drop_fragment:
                        continue
                    short_iv = effective_intervals(short)
                    long_iv = effective_intervals(long)
                    if not fully_contained(short_iv, long_iv):
                        continue
                    zh13 = zh13_of(short)
                    if zh13 == "NA" or zh13 != zh13_of(long):
                        continue
                    if not proteins_are_fragment(effective_proteins(short), effective_proteins(long)):
                        continue
                    decision = "保留_有表达" if expressed(short) else "排除"
                    fragment_rows.append((
                        short, long, models[short]["gene"]["chrom"], models[short]["gene"]["strand"],
                        interval_len(short_iv), interval_len(long_iv),
                        max(len(seq) for seq in effective_proteins(short)),
                        max(len(seq) for seq in effective_proteins(long)),
                        zh13,
                        "yes" if expressed(short) else "no",
                        summary.get(short, {}).get("max_TPM", "NA"),
                        decision,
                    ))
                    if decision == "排除":
                        drop_fragment[short] = long
            active.append((end, start, gene_id))

    removed = set(drop_merge) | set(drop_fragment)
    opposite = [
        gene_id for gene_id in removed
        if models[gene_id]["gene"]["strand"] != models[drop_merge.get(gene_id, drop_fragment.get(gene_id))]["gene"]["strand"]
    ]
    if opposite:
        raise SystemExit(f"{prefix}: refused to remove opposite-strand genes: {opposite}")

    counts = class_counts(os.path.join(final_dir, f"{prefix}.final.gene_class.tsv"), removed)
    print(
        f"{prefix}\tmerge_removed={len(drop_merge)}\treverted={len(revert)}\t"
        f"fragment_removed={len(drop_fragment)}\tfragment_kept={sum(1 for row in fragment_rows if row[-1] == '保留_有表达')}\t"
        f"genes={sum(counts.values())}\thigh={counts['high_confidence']}\t"
        f"below50={counts['expressed_zh13_below_50']}\thomolog={counts['homolog_unexpressed']}",
        flush=True,
    )
    for row in merge_rows:
        print(f"[MERGE] {prefix}\t{row[0]}\t{row[1]}\tzh13={row[8]}\tdropped_expr={row[9]}", flush=True)
    for gene_id in sorted(revert):
        print(f"[REVERT] {prefix}\t{gene_id}\toverlap={conflict_overlap.get(gene_id, 0)}", flush=True)
    for row in fragment_rows:
        print(f"[FRAGMENT] {prefix}\t{row[-1]}\t{row[0]}\tinside\t{row[1]}\t{row[4]}bp/{row[6]}aa", flush=True)
    if dry_run:
        return

    id_rows = []
    for removed_id, kept in drop_merge.items():
        id_rows.append((removed_id, kept, "gene", "冗余合并"))
        for transcript_id in models[removed_id]["order"]:
            id_rows.append((transcript_id, kept, "transcript", "冗余模型未保留"))
    for removed_id, kept in drop_fragment.items():
        id_rows.append((removed_id, kept, "gene", "碎片排除"))
        for transcript_id in models[removed_id]["order"]:
            id_rows.append((transcript_id, kept, "transcript", "碎片排除"))
    for gene_id in sorted(revert):
        id_rows.append((gene_id, gene_id, "gene", "撤销延长"))

    out_blocks = []
    for gene_id, lines in blocks:
        if gene_id in removed:
            continue
        if gene_id in revert:
            out_blocks.extend(renamed_models[gene_id]["lines"])
        else:
            out_blocks.extend(lines)
    gff_tmp = gff_path + ".resolve.tmp"
    with open(gff_tmp, "w") as handle:
        handle.writelines(header)
        handle.writelines(out_blocks)

    def fasta_records(store, renamed_store, kind):
        records = []
        seen = set()
        for gene_id, _lines in blocks:
            if gene_id in removed:
                continue
            if gene_id in revert:
                source = renamed_store
                order = renamed_models[gene_id]["order"]
            else:
                source = store
                order = models[gene_id]["order"]
            for transcript_id in order:
                if transcript_id not in source:
                    raise SystemExit(f"{prefix}: missing {kind} sequence {transcript_id}")
                records.append((transcript_id, source[transcript_id]))
                seen.add(transcript_id)
        return records

    transcript_tmp = write_fasta(
        os.path.join(final_dir, f"{prefix}.final.transcripts.fa"),
        fasta_records(tx_fa, renamed_tx, "transcript"),
    )
    cds_tmp = write_fasta(
        os.path.join(final_dir, f"{prefix}.final.cds.fa"),
        fasta_records(cds_fa, renamed_cds, "CDS"),
    )
    protein_tmp = write_fasta(
        os.path.join(final_dir, f"{prefix}.final.proteins.fa"),
        fasta_records(prot_fa, renamed_prot, "protein"),
    )

    list_path = os.path.join(final_dir, f"{prefix}.final.genes.list")
    with open(list_path) as handle:
        gene_list = [line.strip() for line in handle if line.strip() and line.strip() not in removed]
    with open(list_path + ".resolve.tmp", "w") as handle:
        handle.write("\n".join(gene_list) + "\n")

    class_path = os.path.join(final_dir, f"{prefix}.final.gene_class.tsv")
    with open(class_path) as handle:
        class_lines = handle.readlines()
    class_header = class_lines[0].rstrip("\n").split("\t")
    gene_at = class_header.index("gene_id")
    start_at = class_header.index("start")
    end_at = class_header.index("end")
    with open(class_path + ".resolve.tmp", "w") as handle:
        handle.write(class_lines[0])
        for line in class_lines[1:]:
            parts = line.rstrip("\n").split("\t")
            gene_id = parts[gene_at]
            if gene_id in removed:
                continue
            if gene_id in revert:
                gene = renamed_models[gene_id]["gene"]
                parts[start_at] = str(gene["start"])
                parts[end_at] = str(gene["end"])
                line = "\t".join(parts) + "\n"
            handle.write(line)

    corr_path = os.path.join(final_dir, f"{prefix}.final.corrections.tsv")
    with open(corr_path) as handle:
        corr_lines = handle.readlines()
    corr_header = corr_lines[0].rstrip("\n").split("\t")
    action_at = corr_header.index("action")
    reason_at = corr_header.index("reason")
    gene_corr_at = corr_header.index("gene_id")
    with open(corr_path + ".resolve.tmp", "w") as handle:
        handle.write(corr_lines[0])
        for line in corr_lines[1:]:
            parts = line.rstrip("\n").split("\t")
            if parts[gene_corr_at] in revert and parts[action_at] == "延长":
                parts[action_at] = "撤销延长"
                parts[reason_at] = (
                    "延长后的同链CDS相互重叠达到200bp，原始CDS没有这样的重叠，"
                    "恢复原始模型，禁止再次自动延长"
                )
                line = "\t".join(parts) + "\n"
            handle.write(line)

    map_path = os.path.join(final_dir, f"{prefix}.final.id_map.tsv")
    with open(map_path) as handle:
        map_text = handle.read()
    if not map_text.endswith("\n"):
        map_text += "\n"
    with open(map_path + ".resolve.tmp", "w") as handle:
        handle.write(map_text)
        for old, new, kind, action in id_rows:
            handle.write(f"{old}\t{new}\t{kind}\t{action}\n")

    def write_tsv(name, header_row, rows):
        path = os.path.join(final_dir, f"{prefix}.final.{name}")
        with open(path, "w") as handle:
            handle.write("\t".join(header_row) + "\n")
            for row in rows:
                handle.write("\t".join(str(item) for item in row) + "\n")

    write_tsv(
        "redundant_merged.tsv",
        ["removed_gene", "kept_gene", "chrom", "strand", "start", "end", "cds_bp",
         "protein_aa", "zh13_id", "removed_expressed", "kept_utr_count", "removed_utr_count"],
        merge_rows,
    )
    conflict_rows = []
    for gene_id in sorted(revert):
        row = corrections[gene_id]
        conflict_rows.append((
            gene_id, conflict_overlap.get(gene_id, 0),
            row.get("old_start", ""), row.get("old_end", ""),
            row.get("new_start", ""), row.get("new_end", ""),
            row.get("old_aa", ""), row.get("new_aa", ""),
            row.get("zh13_id", ""),
            "恢复原始模型",
        ))
    write_tsv(
        "extension_conflict.tsv",
        ["gene_id", "new_cds_overlap_bp", "original_start", "original_end",
         "extended_start", "extended_end", "original_aa", "extended_aa", "zh13_id", "action"],
        conflict_rows,
    )
    write_tsv(
        "nested_fragment.tsv",
        ["gene_id", "containing_gene", "chrom", "strand", "short_cds_bp", "long_cds_bp",
         "short_aa", "long_aa", "zh13_id", "expression_supported", "max_tpm", "decision"],
        fragment_rows,
    )

    if STAMP in criteria:
        criteria = criteria[:criteria.index(STAMP)].rstrip() + "\n"
    criteria += (
        f"\n{STAMP}\n"
        "同链且 CDS 序列相同、蛋白序列相同的重复基因号已合并，保留 UTR 转录本更多的基因号。"
        f"本轮合并删除 {len(drop_merge)} 个基因。\n"
        "延长后新 CDS 在同链相互重叠达到 200 bp、而原始 CDS 重叠不足 200 bp 的模型已恢复为 12_renamed 的原始 CDS。"
        f"本轮撤销 {len(revert)} 个基因。\n"
        "同链短 CDS 完全落在同一 ZH13 基因的长 CDS 内、短蛋白是长蛋白的连续片段、且短模型没有表达支持的基因已从本套注释排除。"
        f"本轮排除 {len(drop_fragment)} 个基因。有表达支持的同类模型保留。\n"
        "反义链重叠留给后续筛选，这一步不删除。内含子嵌套，以及外显子结构不同且有独立表达证据的同链重叠，不在这一步删除。位点编号不重排。\n"
        "对应关系记在 final.redundant_merged.tsv、final.extension_conflict.tsv、final.nested_fragment.tsv 和 final.id_map.tsv。\n"
        f"处理后基因数 {sum(counts.values())}，其中 high_confidence {counts['high_confidence']}，"
        f"expressed_zh13_below_50 {counts['expressed_zh13_below_50']}，"
        f"homolog_unexpressed {counts['homolog_unexpressed']}。\n"
    )
    with open(criteria_path + ".resolve.tmp", "w") as handle:
        handle.write(criteria)

    written = []
    with open(gff_tmp) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 9 and parts[2] == "gene":
                written.append(attr(parts[8]).get("ID", ""))
    if written != gene_list:
        raise SystemExit(f"{prefix}: GFF gene order does not match genes.list")
    gff_tx = []
    with open(gff_tmp) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 9 and parts[2] == "mRNA":
                gff_tx.append(attr(parts[8]).get("ID", ""))
    fasta_tmp = {
        "transcripts": transcript_tmp,
        "cds": cds_tmp,
        "proteins": protein_tmp,
    }
    for kind, tmp_path in fasta_tmp.items():
        ids = set(load_fasta(tmp_path))
        if set(gff_tx) != ids:
            missing = sorted(set(gff_tx) - ids)[:5]
            extra = sorted(ids - set(gff_tx))[:5]
            raise SystemExit(f"{prefix}: {kind} IDs do not match GFF transcripts missing={missing} extra={extra}")
    restored_cds = load_fasta(cds_tmp)
    for gene_id in revert:
        for transcript_id in renamed_models[gene_id]["order"]:
            if restored_cds[transcript_id] != renamed_cds[transcript_id]:
                raise SystemExit(f"{prefix}: restored CDS does not match 12_renamed for {transcript_id}")

    os.replace(gff_tmp, gff_path)
    os.replace(transcript_tmp, os.path.join(final_dir, f"{prefix}.final.transcripts.fa"))
    os.replace(cds_tmp, os.path.join(final_dir, f"{prefix}.final.cds.fa"))
    os.replace(protein_tmp, os.path.join(final_dir, f"{prefix}.final.proteins.fa"))
    os.replace(list_path + ".resolve.tmp", list_path)
    os.replace(class_path + ".resolve.tmp", class_path)
    os.replace(corr_path + ".resolve.tmp", corr_path)
    os.replace(map_path + ".resolve.tmp", map_path)
    os.replace(criteria_path + ".resolve.tmp", criteria_path)
    print(f"{prefix}\twrote {len(written)} genes", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--downstream", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    resolve(args.downstream, args.prefix, args.dry_run)


if __name__ == "__main__":
    main()
