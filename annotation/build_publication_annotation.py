#!/usr/bin/env python3
"""Build publication-ready functional annotation table and annotated GFF3."""

from __future__ import annotations

import csv
import os
import re
import sys
from collections import defaultdict

TX_SUFFIX = re.compile(r"\.[tm]\d+")
NUC_RE = re.compile(r"^Chr(0[1-9]|1[0-9]|20)$", re.I)
ZH13_IDENT_KEYS = (
    "ZH13v3.0_mRNA_nuclear_20250729_pro_identity",
    "ZH13_identity",
)


def esc_gff(val: str, max_len: int = 500) -> str:
    if val is None:
        return ""
    s = str(val).strip()
    if not s or s.upper() == "NA":
        return ""
    s = s.replace("\t", " ").replace("\n", " ").replace(";", ",")
    return s[:max_len]


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


def pick_product(row: dict[str, str]) -> str:
    for key in ("eggNOG_description", "SwissProt_description"):
        v = esc_gff(row.get(key, ""))
        if v:
            return v
    return "hypothetical protein"


def pick_name(row: dict[str, str], gid: str) -> str:
    v = esc_gff(row.get("eggNOG_preferred_name", ""))
    if v:
        return v
    v = esc_gff(row.get("SwissProt_hit", ""))
    if v:
        return v
    return gid


def merged_go_terms(row: dict[str, str]) -> str:
    terms: set[str] = set()
    for key in ("eggNOG_GO", "InterPro_GO"):
        v = row.get(key, "")
        if not v or str(v).upper() == "NA":
            continue
        for item in re.split(r"[|,]", str(v)):
            item = item.strip()
            if item.startswith("GO:"):
                terms.add(item)
    return ",".join(sorted(terms))


def build_dbxref(row: dict[str, str]) -> str:
    parts = []
    pn = esc_gff(row.get("eggNOG_preferred_name", ""))
    if pn:
        parts.append(f"eggNOG:{pn}")
    sp = esc_gff(row.get("SwissProt_hit", ""))
    if sp and sp != "NA":
        parts.append(f"UniProtKB:{sp}")
    ipr = row.get("InterPro_IDs", "")
    if ipr and str(ipr).upper() != "NA":
        for x in str(ipr).split(","):
            x = x.strip()
            if x:
                parts.append(f"InterPro:{x}")
    pfam = row.get("InterPro_Pfam", "")
    if pfam and str(pfam).upper() != "NA":
        for x in str(pfam).split(","):
            x = x.strip()
            if x.startswith("PF"):
                parts.append(f"Pfam:{x}")
    return ",".join(parts)


def build_note(row: dict[str, str]) -> str:
    bits = []
    z13 = esc_gff(row.get("ZH13_hit", ""))
    w82 = esc_gff(row.get("Wm82_hit", ""))
    if z13:
        bits.append(f"ZH13_homolog={z13}")
    if w82:
        bits.append(f"Wm82_homolog={w82}")
    te = row.get("TE_overlap_ratio", "NA")
    if te not in ("", "NA"):
        bits.append(f"TE_overlap_ratio={te}")
    conf = row.get("confidence_level", "")
    if conf:
        bits.append(f"confidence={conf}")
    return "; ".join(bits)


def confidence_pass(row: dict[str, str], min_conf: str, max_te: float, require_expr: bool) -> bool:
    conf = row.get("confidence_level", "low_confidence")
    rank = {
        "high_confidence": 3,
        "medium_confidence": 2,
        "low_confidence": 1,
        "low_confidence_possible_TE": 0,
    }
    if rank.get(conf, 0) < rank.get(min_conf, 2):
        return False
    try:
        te = float(row.get("TE_overlap_ratio", 0) or 0)
    except ValueError:
        te = 0.0
    if te > max_te:
        return False
    if require_expr and row.get("expression_supported", "NA") != "yes":
        return False
    return True


def is_nuclear_chrom(chrom: str) -> bool:
    return bool(NUC_RE.match(chrom or ""))


def zh13_pident(row: dict[str, str]) -> float:
    for key in ZH13_IDENT_KEYS:
        raw = row.get(key)
        if raw in (None, "", "NA", "na", "."):
            continue
        try:
            return float(raw)
        except ValueError:
            continue
    return 0.0


def hc_pass(
    row: dict[str, str],
    min_conf: str,
    max_te: float,
    require_expr: bool,
    nuclear_only: bool,
    min_zh13: float,
) -> bool:
    if nuclear_only and not is_nuclear_chrom(row.get("chrom", "")):
        return False
    if not confidence_pass(row, min_conf, max_te, require_expr):
        return False
    if min_zh13 > 0 and zh13_pident(row) < min_zh13:
        return False
    return True


def load_summary(path: str) -> dict[str, dict[str, str]]:
    data: dict[str, dict[str, str]] = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            data[row["gene_id"]] = row
    return data


def write_gene_function_table(rows: dict[str, dict[str, str]], out_path: str) -> None:
    header = [
        "gene_id", "chrom", "start", "end", "strand",
        "gene_name", "product", "Dbxref", "GO_terms", "KEGG_ko", "KEGG_pathway", "PFAM",
        "InterPro_IDs", "InterPro_Pfam", "InterPro_GO",
        "SwissProt_hit", "ZH13_homolog", "Wm82_homolog",
        "TE_overlap_ratio", "max_TPM", "expression_supported", "annotation_confidence",
    ]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for gid in sorted(
            rows,
            key=lambda x: (rows[x].get("chrom", ""), int(rows[x].get("start", 0) or 0)),
        ):
            r = rows[gid]
            w.writerow(
                {
                    "gene_id": gid,
                    "chrom": r.get("chrom", ""),
                    "start": r.get("start", ""),
                    "end": r.get("end", ""),
                    "strand": r.get("strand", ""),
                    "gene_name": pick_name(r, gid),
                    "product": pick_product(r),
                    "Dbxref": build_dbxref(r),
                    "GO_terms": merged_go_terms(r) or "NA",
                    "KEGG_ko": r.get("eggNOG_KEGG_KO", "NA"),
                    "KEGG_pathway": r.get("eggNOG_KEGG_pathway", "NA"),
                    "PFAM": r.get("eggNOG_PFAM", "NA"),
                    "InterPro_IDs": r.get("InterPro_IDs", "NA"),
                    "InterPro_Pfam": r.get("InterPro_Pfam", "NA"),
                    "InterPro_GO": r.get("InterPro_GO", "NA"),
                    "SwissProt_hit": r.get("SwissProt_hit", "NA"),
                    "ZH13_homolog": r.get("ZH13_hit", "NA"),
                    "Wm82_homolog": r.get("Wm82_hit", "NA"),
                    "TE_overlap_ratio": r.get("TE_overlap_ratio", "NA"),
                    "max_TPM": r.get("max_TPM", "NA"),
                    "expression_supported": r.get("expression_supported", "NA"),
                    "annotation_confidence": r.get("confidence_level", "NA"),
                }
            )


def feature_gene_id(feature: str, attrs: dict[str, str], tx2gene: dict[str, str]) -> str | None:
    if feature == "gene":
        return attrs.get("ID")
    if feature in ("mRNA", "transcript"):
        tid = attrs.get("ID")
        if tid and tid in tx2gene:
            return tx2gene[tid]
        parent = attrs.get("Parent", "")
        return parent.split(",")[0] if parent else None
    parent = attrs.get("Parent", "")
    for pid in parent.split(","):
        if not pid:
            continue
        if pid in tx2gene:
            return tx2gene[pid]
        if TX_SUFFIX.search(pid):
            return tx2gene.get(pid, TX_SUFFIX.split(pid)[0])
    return None


def annotate_gff(
    in_gff: str,
    out_gff: str,
    rows: dict[str, dict[str, str]],
    gene_filter: set[str] | None = None,
) -> int:
    tx2gene: dict[str, str] = {}
    n_genes = 0

    with open(in_gff) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            if p[2] in ("mRNA", "transcript"):
                a = parse_attr(p[8])
                tid, parent = a.get("ID"), a.get("Parent", "")
                if tid and parent:
                    tx2gene[tid] = parent.split(",")[0]

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

            feature = p[2]
            a = parse_attr(p[8])
            gid = feature_gene_id(feature, a, tx2gene)

            if gene_filter is not None:
                if gid is None or gid not in gene_filter:
                    continue

            if gid and gid in rows:
                r = rows[gid]
                if feature == "gene":
                    n_genes += 1
                    a["Name"] = pick_name(r, gid)
                    a["product"] = pick_product(r)
                    dbx = build_dbxref(r)
                    if dbx:
                        a["Dbxref"] = dbx
                    go = merged_go_terms(r)
                    if go:
                        a["Ontology_term"] = go
                    note = build_note(r)
                    if note:
                        a["Note"] = note
                    a["annotation_confidence"] = esc_gff(r.get("confidence_level", ""))
                elif feature in ("mRNA", "transcript"):
                    a["Name"] = pick_name(r, gid)
                    a["product"] = pick_product(r)
                    a["annotation_confidence"] = esc_gff(r.get("confidence_level", ""))

            p[8] = fmt_attr(a)
            fout.write("\t".join(p) + "\n")

    return n_genes


def extract_proteins(protein_fa: str, gene_list: set[str], tx2gene: dict[str, str], out_fa: str) -> int:
    wanted_tx = {tid for tid, g in tx2gene.items() if g in gene_list}
    n = 0
    with open(protein_fa) as fin, open(out_fa, "w") as fout:
        keep = False
        for line in fin:
            if line.startswith(">"):
                sid = line[1:].split()[0]
                gid = tx2gene.get(sid, TX_SUFFIX.split(sid)[0] if TX_SUFFIX.search(sid) else sid)
                keep = sid in wanted_tx or gid in gene_list
                if keep:
                    n += 1
                if keep:
                    fout.write(line)
            elif keep:
                fout.write(line)
    return n


def main() -> int:
    in_gff = os.environ["PUB_STD_GFF"]
    summary = os.environ["PUB_SUMMARY"]
    protein_fa = os.environ["PUB_PROTEIN_FA"]
    out_dir = os.environ["PUB_OUT_DIR"]
    prefix = os.environ["PUB_PREFIX"]
    min_conf = os.environ.get("PUB_MIN_CONFIDENCE", "medium_confidence")
    max_te = float(os.environ.get("PUB_MAX_TE_RATIO", "0.5"))
    require_expr = os.environ.get("PUB_REQUIRE_EXPRESSION", "0") == "1"
    nuclear_only = os.environ.get("PUB_NUCLEAR_ONLY", "1") == "1"
    min_zh13 = float(os.environ.get("PUB_MIN_ZH13_PIDENT", "50"))

    os.makedirs(out_dir, exist_ok=True)

    gene_table = os.path.join(out_dir, f"{prefix}.gene_functional_annotation.tsv")
    annot_gff = os.path.join(out_dir, f"{prefix}.annotated.gff3")
    hc_list = os.path.join(out_dir, f"{prefix}.high_confidence.genes.list")
    hc_gff = os.path.join(out_dir, f"{prefix}.high_confidence.annotated.gff3")
    hc_protein = os.path.join(out_dir, f"{prefix}.high_confidence.proteins.fa")
    readme = os.path.join(out_dir, f"{prefix}.annotation_README.txt")

    rows = load_summary(summary)
    write_gene_function_table(rows, gene_table)

    n_all = annotate_gff(in_gff, annot_gff, rows, gene_filter=None)
    print(f"[INFO] Annotated GFF3: {annot_gff} ({n_all} genes with annotation rows)")

    hc_genes = {
        gid
        for gid, r in rows.items()
        if hc_pass(r, min_conf, max_te, require_expr, nuclear_only, min_zh13)
    }
    with open(hc_list, "w") as f:
        for gid in sorted(hc_genes, key=lambda x: (rows[x].get("chrom", ""), int(rows[x].get("start", 0) or 0))):
            f.write(f"{gid}\n")

    n_hc = annotate_gff(in_gff, hc_gff, rows, gene_filter=hc_genes)
    print(f"[INFO] High-confidence GFF3: {hc_gff} ({n_hc} genes)")

    tx2gene: dict[str, str] = {}
    with open(in_gff) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] not in ("mRNA", "transcript"):
                continue
            a = parse_attr(p[8])
            if a.get("ID") and a.get("Parent"):
                tx2gene[a["ID"]] = a["Parent"].split(",")[0]

    n_prot = extract_proteins(protein_fa, hc_genes, tx2gene, hc_protein)
    print(f"[INFO] High-confidence proteins: {hc_protein} ({n_prot} sequences)")

    with open(readme, "w") as f:
        f.write(
            f"""Publication-ready annotation outputs for {prefix}
================================================================

1. {os.path.basename(gene_table)}
   Gene-level functional annotation table (main supplementary table).

2. {os.path.basename(annot_gff)}
   Structural GFF3 with functional attributes on gene/mRNA features:
   Name, product, Dbxref, Ontology_term (GO), Note, annotation_confidence.

3. {os.path.basename(hc_list)}
   Gene IDs passing quality filters (Chr01–20, confidence >= {min_conf}, TE overlap <= {max_te}, ZH13 v3 pident >= {min_zh13}).

4. {os.path.basename(hc_gff)}
   Annotated GFF3 for high-confidence genes only.

5. {os.path.basename(hc_protein)}
   Protein sequences for high-confidence genes.

Filter settings:
  PUB_MIN_CONFIDENCE={min_conf}
  PUB_MAX_TE_RATIO={max_te}
  PUB_REQUIRE_EXPRESSION={require_expr}
  PUB_NUCLEAR_ONLY={nuclear_only}
  PUB_MIN_ZH13_PIDENT={min_zh13}

Upstream summary: {summary}
"""
        )
    print(f"[INFO] README: {readme}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
