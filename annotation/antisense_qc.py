#!/usr/bin/env python3
"""Classify antisense CDS overlaps and keep one 13_final gene set.

Pairs stay when both genes have reliable homology. A pair that also has the
same antisense CDS overlap in ZH13 v3 is recorded as higher confidence. Other
antisense genes are removed from the GFF and FASTA. The removed IDs stay in
PREFIX.final.antisense_excluded.tsv and PREFIX.final.id_map.tsv.
"""

import argparse
import os
from collections import defaultdict

ZH13_GFF = "/home/zzh/ref/ZH13_v3/ZH13v3.0_mRNA_nuclear_20250729.gff3"
STAMP = "反义链质量评估已写入"
EXCLUSION_STAMP = "反义链筛选已写入"
IDENTITY_COLUMNS = (
    "ZH13.v2.prot_identity",
    "ZH13v3.0_mRNA_nuclear_20250729_pro_identity",
    "Wm82.a2.v1.prot_identity",
    "glyma.Wm82.gnm4.ann1.T8TQ.protein_identity",
    "glyma.Wm82.gnm6.ann1.PKSW.protein_identity",
)
ZH13_HIT = "ZH13v3.0_mRNA_nuclear_20250729_pro_hit"


def attr(text):
    out = {}
    for part in text.strip().rstrip(";").split(";"):
        if "=" in part:
            key, value = part.split("=", 1)
            out[key] = value
    return out


def locus(hit):
    if not hit or hit == "NA":
        return "NA"
    return hit.split(".", 1)[0] if hit.startswith("SoyZH13") and ".m" in hit else hit.split(".m")[0]


def number(text):
    if text is None or text in ("", "NA", "."):
        return None
    try:
        return float(text)
    except ValueError:
        return None


def merge(intervals):
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


def span_len(intervals):
    return sum(end - start + 1 for start, end in intervals)


def overlap_len(left, right):
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


def read_table(path):
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        rows = {}
        for line in handle:
            cells = line.rstrip("\n").split("\t")
            rows[cells[0]] = dict(zip(header, cells))
        return rows


def load_query(gff_path):
    genes = {}
    parent = {}
    cds = defaultdict(list)
    cds_count = defaultdict(int)
    with open(gff_path) as handle:
        for line in handle:
            if not line or line[0] == "#":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] not in ("gene", "mRNA", "CDS"):
                continue
            values = attr(parts[8])
            if parts[2] == "gene":
                genes[values["ID"]] = {
                    "chrom": parts[0],
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "strand": parts[6],
                }
            elif parts[2] == "mRNA":
                parent[values["ID"]] = values.get("Parent")
            else:
                cds[values["Parent"]].append((int(parts[3]), int(parts[4])))
                cds_count[values["Parent"]] += 1
    flags = defaultdict(list)
    for transcript, intervals in cds.items():
        gene_id = parent.get(transcript)
        if gene_id not in genes:
            continue
        genes[gene_id].setdefault("cds", []).extend(intervals)
        flags[gene_id].append(cds_count[transcript] == 1)
        genes[gene_id]["cds_exons"] = max(
            genes[gene_id].get("cds_exons", 0), cds_count[transcript]
        )
    coding = []
    for gene_id, gene in genes.items():
        if not gene.get("cds"):
            continue
        gene["id"] = gene_id
        gene["cds"] = merge(gene["cds"])
        gene["cds_bp"] = span_len(gene["cds"])
        gene["single_exon"] = bool(flags[gene_id]) and all(flags[gene_id])
        coding.append(gene)
    return coding


def load_zh13(path):
    genes = {}
    parent = {}
    cds = defaultdict(list)
    with open(path) as handle:
        for line in handle:
            if not line or line[0] == "#":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] not in ("gene", "mRNA", "CDS"):
                continue
            values = attr(parts[8])
            if parts[2] == "gene":
                genes[values["ID"]] = {
                    "chrom": parts[0],
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "strand": parts[6],
                }
            elif parts[2] == "mRNA":
                parent[values["ID"]] = values.get("Parent")
            else:
                owner = parent.get(values["Parent"], values["Parent"])
                cds[owner].append((int(parts[3]), int(parts[4])))
    for gene_id, gene in genes.items():
        gene["cds"] = merge(cds.get(gene_id, []))
    return genes


def zh13_structure(left_id, right_id, zh13):
    if left_id == "NA" or right_id == "NA":
        return "not_both_zh13", 0
    if left_id == right_id:
        return "same_locus", 0
    left = zh13.get(left_id)
    right = zh13.get(right_id)
    if left is None or right is None:
        return "not_both_zh13", 0
    if left["chrom"] != right["chrom"]:
        return "different_chrom", 0
    cds_overlap = overlap_len(left["cds"], right["cds"]) if left["cds"] and right["cds"] else 0
    spans = not (left["end"] < right["start"] or right["end"] < left["start"])
    if left["strand"] == right["strand"]:
        if cds_overlap:
            return "same_strand_cds_overlap", cds_overlap
        return "same_strand_no_cds_overlap", 0
    if cds_overlap:
        return "antisense_cds_overlap", cds_overlap
    if spans:
        return "antisense_span_only", 0
    return "no_overlap", 0


def evidence(gene_class, summary, gene_id):
    row = gene_class[gene_id]
    info = summary.get(gene_id, {})
    identity = number(row.get("zh13_v3_identity"))
    identities = [number(info.get(column)) for column in IDENTITY_COLUMNS]
    identities = [value for value in identities if value is not None]
    best = max(identities) if identities else None
    reliable = row["class"] in ("high_confidence", "homolog_unexpressed")
    if reliable and row["class"] == "high_confidence":
        basis = "ZH13_v3"
    elif reliable:
        basis = row.get("support_ref", "NA")
    else:
        basis = "none"
    if reliable:
        coding = "RELIABLE"
    elif best is not None and best >= 50:
        coding = "PARTIAL"
    else:
        coding = "UNVERIFIED"
    tpm = number(info.get("max_TPM")) or 0.0
    quantification = "AMBIGUOUS" if tpm >= 1 else "NO_EXPRESSION"
    return {
        "class": row["class"],
        "zh13_id": locus(info.get(ZH13_HIT, "NA")),
        "zh13_identity": row.get("zh13_v3_identity", "NA"),
        "reliable": reliable,
        "basis": basis,
        "coding": coding,
        "tpm": f"{tpm:.4f}",
        "quantification": quantification,
        "identity_value": identity,
    }


def classify(left, right, structure):
    both = left["reliable"] and right["reliable"]
    different = (
        left["zh13_id"] != "NA"
        and right["zh13_id"] != "NA"
        and left["zh13_id"] != right["zh13_id"]
    )
    if both and different and structure == "antisense_cds_overlap":
        return "较高", "保留"
    if both and different:
        return "中等至较高", "保留，核对共线性和结构"
    if both and left["zh13_id"] != "NA" and left["zh13_id"] == right["zh13_id"]:
        return "其他", "排除"
    if both:
        return "中等至较高", "保留"
    unsupported = left if not left["reliable"] else right
    partner = right if unsupported is left else left
    if (
        partner["class"] == "high_confidence"
        and unsupported["single_exon"]
        and not unsupported["reliable"]
    ):
        return "低或未确定", "排除"
    return "其他", "排除"


def write_table(path, header, rows):
    temporary = path + ".tmp"
    with open(temporary, "w") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[name]) for name in header) + "\n")
    os.replace(temporary, path)


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


def write_fasta(path, records):
    temporary = path + ".antisense.tmp"
    with open(temporary, "w") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")
            for pos in range(0, len(seq), 70):
                handle.write(seq[pos:pos + 70] + "\n")
    return temporary


def gff_blocks(path):
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


def transcript_ids(lines):
    found = []
    for line in lines:
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 9 and parts[2] == "mRNA":
            found.append(attr(parts[8]).get("ID", ""))
    return found


def side_fields(row, side):
    return {
        "gene_id": row[f"gene_{side}"],
        "partner_gene_id": row["gene_b" if side == "a" else "gene_a"],
        "gene_class": row[f"class_{side}"],
        "partner_class": row["class_b" if side == "a" else "class_a"],
        "reliable": row[f"reliable_homology_{side}"] == "yes",
        "single_exon": row[f"single_exon_{side}"],
        "cds_bp": row[f"cds_bp_{side}"],
        "cds_exons": int(row[f"cds_exons_{side}"]),
        "zh13_id": row[f"zh13_id_{side}"],
        "partner_zh13_id": row["zh13_id_b" if side == "a" else "zh13_id_a"],
        "max_tpm": row[f"max_tpm_{side}"],
        "overlap_fraction": row[f"overlap_fraction_{side}"],
    }


def weakness(row, side):
    return (
        row[f"single_exon_{side}"] == "yes",
        -int(row[f"cds_exons_{side}"]),
        -int(row[f"cds_bp_{side}"]),
        row[f"gene_{side}"],
    )


def exclusion_rows(pairs):
    protected = set()
    for row in pairs:
        if row["evidence_tier"] in ("较高", "中等至较高"):
            protected.add(row["gene_a"])
            protected.add(row["gene_b"])
    chosen = defaultdict(list)

    def add(side_row, reason):
        gene_id = side_row["gene_id"]
        if gene_id in protected:
            return
        chosen[gene_id].append((side_row, reason))

    for row in pairs:
        tier = row["evidence_tier"]
        if tier in ("较高", "中等至较高"):
            continue
        left = side_fields(row, "a")
        right = side_fields(row, "b")
        if tier == "低或未确定":
            add(right, "低或未确定：单外显子反义模型缺少可靠同源")
            continue
        if left["reliable"] and right["reliable"]:
            drop = "a" if weakness(row, "a") > weakness(row, "b") else "b"
            add(side_fields(row, drop), "其他：两侧对应同一个ZH13基因，保留外显子更多的一侧")
            continue
        if left["reliable"] and not right["reliable"]:
            add(right, "其他：自身没有可靠同源")
            continue
        if right["reliable"] and not left["reliable"]:
            add(left, "其他：自身没有可靠同源")
            continue
        add(left, "其他：两侧都没有可靠同源")
        add(right, "其他：两侧都没有可靠同源")
    return chosen


def apply_exclusion(final, prefix, removed):
    gff_path = os.path.join(final, f"{prefix}.final.gff3")
    header, blocks = gff_blocks(gff_path)
    transcripts = {}
    for gene_id, lines in blocks:
        transcripts[gene_id] = transcript_ids(lines)
    missing = sorted(set(removed) - set(transcripts))
    if missing:
        raise SystemExit(f"{prefix}: exclusion ids missing from GFF {missing[:5]}")
    for gene_id, records in removed.items():
        classes = {item[0]["gene_class"] for item in records}
        reasons = {item[1] for item in records}
        if "high_confidence" in classes and not any(reason.startswith("其他：两侧对应同一个ZH13基因") for reason in reasons):
            raise SystemExit(f"{prefix}: refusing to remove high_confidence gene {gene_id}")
        if "homolog_unexpressed" in classes:
            raise SystemExit(f"{prefix}: refusing to remove homolog_unexpressed gene {gene_id}")

    gff_tmp = gff_path + ".antisense.tmp"
    with open(gff_tmp, "w") as handle:
        handle.writelines(header)
        for gene_id, lines in blocks:
            if gene_id not in removed:
                handle.writelines(lines)

    def kept_records(fasta_path):
        seqs = load_fasta(fasta_path)
        records = []
        for gene_id, lines in blocks:
            if gene_id in removed:
                continue
            for transcript_id in transcripts[gene_id]:
                if transcript_id not in seqs:
                    raise SystemExit(f"{prefix}: missing sequence {transcript_id}")
                records.append((transcript_id, seqs[transcript_id]))
        return records

    transcript_tmp = write_fasta(
        os.path.join(final, f"{prefix}.final.transcripts.fa"),
        kept_records(os.path.join(final, f"{prefix}.final.transcripts.fa")),
    )
    cds_tmp = write_fasta(
        os.path.join(final, f"{prefix}.final.cds.fa"),
        kept_records(os.path.join(final, f"{prefix}.final.cds.fa")),
    )
    protein_tmp = write_fasta(
        os.path.join(final, f"{prefix}.final.proteins.fa"),
        kept_records(os.path.join(final, f"{prefix}.final.proteins.fa")),
    )

    list_path = os.path.join(final, f"{prefix}.final.genes.list")
    with open(list_path) as handle:
        gene_list = [line.strip() for line in handle if line.strip() and line.strip() not in removed]
    with open(list_path + ".antisense.tmp", "w") as handle:
        handle.write("\n".join(gene_list) + "\n")

    class_path = os.path.join(final, f"{prefix}.final.gene_class.tsv")
    with open(class_path) as handle:
        class_lines = handle.readlines()
    class_header = class_lines[0].rstrip("\n").split("\t")
    gene_at = class_header.index("gene_id")
    class_at = class_header.index("class")
    counts = defaultdict(int)
    with open(class_path + ".antisense.tmp", "w") as handle:
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
    with open(map_path + ".antisense.tmp", "w") as handle:
        handle.write(map_text)
        for gene_id in sorted(removed):
            partners = []
            for side, _reason in removed[gene_id]:
                partner = side["partner_gene_id"]
                if partner not in removed and partner not in partners:
                    partners.append(partner)
            kept = partners[0] if partners else "NA"
            handle.write(f"{gene_id}\t{kept}\tgene\t反义排除\n")
            for transcript_id in transcripts[gene_id]:
                handle.write(f"{transcript_id}\t{kept}\ttranscript\t反义排除\n")

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
    for tmp_path, kind in (
        (transcript_tmp, "transcripts"),
        (cds_tmp, "cds"),
        (protein_tmp, "proteins"),
    ):
        ids = set(load_fasta(tmp_path))
        if ids != set(gff_tx):
            raise SystemExit(f"{prefix}: {kind} IDs do not match GFF")

    os.replace(gff_tmp, gff_path)
    os.replace(transcript_tmp, os.path.join(final, f"{prefix}.final.transcripts.fa"))
    os.replace(cds_tmp, os.path.join(final, f"{prefix}.final.cds.fa"))
    os.replace(protein_tmp, os.path.join(final, f"{prefix}.final.proteins.fa"))
    os.replace(list_path + ".antisense.tmp", list_path)
    os.replace(class_path + ".antisense.tmp", class_path)
    os.replace(map_path + ".antisense.tmp", map_path)
    return counts, len(written)


def run_sample(downstream, prefix, zh13):
    final = os.path.join(downstream, "13_final")
    criteria_path = os.path.join(final, f"{prefix}.final.criteria.txt")
    with open(criteria_path) as handle:
        current_criteria = handle.read()
    if EXCLUSION_STAMP in current_criteria:
        print(f"{prefix}\talready_filtered", flush=True)
        return
    coding = load_query(os.path.join(final, f"{prefix}.final.gff3"))
    classes = read_table(os.path.join(final, f"{prefix}.final.gene_class.tsv"))
    summary = read_table(
        os.path.join(downstream, "12_renamed", f"{prefix}.renamed.annotation.summary.tsv")
    )
    missing = [gene["id"] for gene in coding if gene["id"] not in classes]
    if missing:
        raise SystemExit(f"{prefix}: {len(missing)} genes missing from gene_class")

    by_chrom = defaultdict(list)
    looked = {}
    for gene in coding:
        looked[gene["id"]] = evidence(classes, summary, gene["id"])
        looked[gene["id"]]["single_exon"] = gene["single_exon"]
        by_chrom[gene["chrom"]].append(gene)

    pair_header = [
        "pair_id", "chrom", "overlap_bp",
        "gene_a", "strand_a", "class_a", "cds_bp_a", "cds_exons_a", "single_exon_a",
        "zh13_id_a", "zh13_identity_a", "reliable_homology_a", "homology_basis_a",
        "coding_support_a", "max_tpm_a", "quantification_status_a", "overlap_fraction_a",
        "gene_b", "strand_b", "class_b", "cds_bp_b", "cds_exons_b", "single_exon_b",
        "zh13_id_b", "zh13_identity_b", "reliable_homology_b", "homology_basis_b",
        "coding_support_b", "max_tpm_b", "quantification_status_b", "overlap_fraction_b",
        "zh13_structure", "zh13_cds_overlap_bp", "evidence_tier", "action",
        "independent_strand_evidence",
    ]
    model_header = [
        "gene_id", "partner_gene_id", "chrom", "strand", "gene_class", "partner_class",
        "single_exon", "cds_bp", "cds_overlap_bp", "overlap_fraction_of_cds",
        "zh13_id", "zh13_identity", "max_tpm", "evidence_tier", "action",
        "antisense_overlap", "independent_strand_evidence", "coding_support",
        "quantification_status",
    ]
    pairs = []
    models = []
    counts = defaultdict(int)
    for chrom, genes in by_chrom.items():
        genes.sort(key=lambda item: (item["start"], item["end"], item["id"]))
        active = []
        for gene in genes:
            active = [item for item in active if item["end"] >= gene["start"]]
            for other in active:
                if other["strand"] == gene["strand"]:
                    continue
                overlap = overlap_len(other["cds"], gene["cds"])
                if overlap <= 0:
                    continue
                sides = []
                for item in (other, gene):
                    info = looked[item["id"]]
                    sides.append((item, info))
                sides.sort(key=lambda item: (
                    item[1]["class"] != "high_confidence",
                    not item[1]["reliable"],
                    item[0]["single_exon"],
                    -item[0]["cds_exons"],
                    item[0]["id"],
                ))
                (gene_a, info_a), (gene_b, info_b) = sides
                structure, zh13_overlap = zh13_structure(info_a["zh13_id"], info_b["zh13_id"], zh13)
                tier, action = classify(
                    {**info_a, "single_exon": gene_a["single_exon"]},
                    {**info_b, "single_exon": gene_b["single_exon"]},
                    structure,
                )
                counts[tier] += 1
                if overlap >= 200:
                    counts[tier + ">=200"] += 1
                fraction_a = overlap / gene_a["cds_bp"]
                fraction_b = overlap / gene_b["cds_bp"]
                pairs.append({
                    "pair_id": f"{prefix}:{gene_a['id']}:{gene_b['id']}",
                    "chrom": chrom,
                    "overlap_bp": overlap,
                    "gene_a": gene_a["id"],
                    "strand_a": gene_a["strand"],
                    "class_a": info_a["class"],
                    "cds_bp_a": gene_a["cds_bp"],
                    "cds_exons_a": gene_a["cds_exons"],
                    "single_exon_a": "yes" if gene_a["single_exon"] else "no",
                    "zh13_id_a": info_a["zh13_id"],
                    "zh13_identity_a": info_a["zh13_identity"],
                    "reliable_homology_a": "yes" if info_a["reliable"] else "no",
                    "homology_basis_a": info_a["basis"],
                    "coding_support_a": info_a["coding"],
                    "max_tpm_a": info_a["tpm"],
                    "quantification_status_a": info_a["quantification"],
                    "overlap_fraction_a": f"{fraction_a:.4f}",
                    "gene_b": gene_b["id"],
                    "strand_b": gene_b["strand"],
                    "class_b": info_b["class"],
                    "cds_bp_b": gene_b["cds_bp"],
                    "cds_exons_b": gene_b["cds_exons"],
                    "single_exon_b": "yes" if gene_b["single_exon"] else "no",
                    "zh13_id_b": info_b["zh13_id"],
                    "zh13_identity_b": info_b["zh13_identity"],
                    "reliable_homology_b": "yes" if info_b["reliable"] else "no",
                    "homology_basis_b": info_b["basis"],
                    "coding_support_b": info_b["coding"],
                    "max_tpm_b": info_b["tpm"],
                    "quantification_status_b": info_b["quantification"],
                    "overlap_fraction_b": f"{fraction_b:.4f}",
                    "zh13_structure": structure,
                    "zh13_cds_overlap_bp": zh13_overlap,
                    "evidence_tier": tier,
                    "action": action,
                    "independent_strand_evidence": "NO",
                })
                if tier == "低或未确定":
                    models.append({
                        "gene_id": gene_b["id"],
                        "partner_gene_id": gene_a["id"],
                        "chrom": chrom,
                        "strand": gene_b["strand"],
                        "gene_class": info_b["class"],
                        "partner_class": info_a["class"],
                        "single_exon": "yes",
                        "cds_bp": gene_b["cds_bp"],
                        "cds_overlap_bp": overlap,
                        "overlap_fraction_of_cds": f"{fraction_b:.4f}",
                        "zh13_id": info_b["zh13_id"],
                        "zh13_identity": info_b["zh13_identity"],
                        "max_tpm": info_b["tpm"],
                        "evidence_tier": tier,
                        "action": action,
                        "antisense_overlap": info_a["class"],
                        "independent_strand_evidence": "NO",
                        "coding_support": info_b["coding"],
                        "quantification_status": info_b["quantification"],
                    })
            active.append(gene)

    pairs.sort(key=lambda row: (row["chrom"], row["gene_a"], row["gene_b"]))
    models.sort(key=lambda row: (row["chrom"], row["gene_id"], row["partner_gene_id"]))
    removed = exclusion_rows(pairs)
    excluded_header = [
        "removed_gene", "kept_gene", "chrom", "strand", "gene_class", "partner_class",
        "evidence_tier", "reason", "single_exon", "cds_bp", "cds_overlap_bp",
        "overlap_fraction_of_cds", "zh13_id", "partner_zh13_id", "max_tpm",
    ]
    excluded = []
    for gene_id in sorted(removed):
        for side, reason in removed[gene_id]:
            partner = side["partner_gene_id"]
            kept = partner if partner not in removed else "NA"
            pair = next(
                row for row in pairs
                if {row["gene_a"], row["gene_b"]} == {gene_id, partner}
                and row["evidence_tier"] not in ("较高", "中等至较高")
            )
            excluded.append({
                "removed_gene": gene_id,
                "kept_gene": kept,
                "chrom": pair["chrom"],
                "strand": pair["strand_b"] if pair["gene_b"] == gene_id else pair["strand_a"],
                "gene_class": side["gene_class"],
                "partner_class": side["partner_class"],
                "evidence_tier": pair["evidence_tier"],
                "reason": reason,
                "single_exon": side["single_exon"],
                "cds_bp": side["cds_bp"],
                "cds_overlap_bp": pair["overlap_bp"],
                "overlap_fraction_of_cds": side["overlap_fraction"],
                "zh13_id": side["zh13_id"],
                "partner_zh13_id": side["partner_zh13_id"],
                "max_tpm": side["max_tpm"],
            })
    write_table(os.path.join(final, f"{prefix}.final.antisense_pair_qc.tsv"), pair_header, pairs)
    write_table(os.path.join(final, f"{prefix}.final.antisense_model_qc.tsv"), model_header, models)
    write_table(os.path.join(final, f"{prefix}.final.antisense_excluded.tsv"), excluded_header, excluded)
    class_counts, n_genes = apply_exclusion(final, prefix, removed)
    kept_pairs = counts["较高"] + counts["中等至较高"]
    note = (
        f"{EXCLUSION_STAMP}\n"
        "反义筛选没有另建一套注释，也没有把等级写进 GFF。位点编号不重排。\n"
        "保留两侧都有可靠同源的反义 CDS 重叠。较高还要求两个 ZH13 v3 基因自身有反义 CDS 重叠；"
        "中等至较高是两侧有可靠同源，但 ZH13 v3 没有同样的反义 CDS 重叠。\n"
        "低或未确定的单外显子模型，以及其余没有归入这两档的反义对，已从本套蛋白编码注释排除。"
        "有可靠同源的重叠对象保留。两侧对应同一个 ZH13 基因时，只排除外显子更少的一侧。\n"
        "排除前的成对等级在 final.antisense_pair_qc.tsv，低或未确定模型在 final.antisense_model_qc.tsv，"
        "实际排除的基因在 final.antisense_excluded.tsv 和 final.id_map.tsv。\n"
        "可靠同源：high_confidence 为 ZH13 v3 相似度 >= 50；homolog_unexpressed 为五套大豆蛋白的最高 bitscore 命中同时满足相似度、覆盖和 E-value 阈值。"
        "现有 Salmon 文库类型为非链特异性 IU，independent_strand_evidence=NO，TPM 不能作为反义转录证据。\n"
        f"排除前反义 CDS 重叠 {len(pairs)} 对：较高 {counts['较高']}，中等至较高 {counts['中等至较高']}，"
        f"低或未确定 {counts['低或未确定']}，其他 {counts['其他']}。保留 {kept_pairs} 对。\n"
        f"本轮排除 {len(removed)} 个基因。处理后基因数 {n_genes}，其中 high_confidence {class_counts['high_confidence']}，"
        f"expressed_zh13_below_50 {class_counts['expressed_zh13_below_50']}，"
        f"homolog_unexpressed {class_counts['homolog_unexpressed']}。\n"
    )
    criteria = current_criteria
    if STAMP in criteria:
        criteria = criteria[:criteria.index(STAMP)].rstrip() + "\n"
    if not criteria.endswith("\n"):
        criteria += "\n"
    with open(criteria_path, "w") as handle:
        handle.write(criteria + "\n" + note)
    print(
        f"{prefix}\tpairs={len(pairs)}\thigh={counts['较高']}\t"
        f"mid={counts['中等至较高']}\tlow={counts['低或未确定']}\t"
        f"other={counts['其他']}\tremoved={len(removed)}\tgenes={n_genes}\t"
        f"hc={class_counts['high_confidence']}\tbelow={class_counts['expressed_zh13_below_50']}\t"
        f"homolog={class_counts['homolog_unexpressed']}",
        flush=True,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--downstream")
    parser.add_argument("--prefix")
    parser.add_argument("--zh13-gff", default=ZH13_GFF)
    parser.add_argument("--all", action="store_true")
    args = parser.parse_args()
    zh13 = load_zh13(args.zh13_gff)
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
            run_sample(downstream, prefix, zh13)
        return
    if not args.downstream or not args.prefix:
        parser.error("give --downstream and --prefix, or --all")
    run_sample(args.downstream, args.prefix, zh13)


if __name__ == "__main__":
    main()
