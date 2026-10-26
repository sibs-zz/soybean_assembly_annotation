#!/usr/bin/env bash
# Build the final protein-coding annotation from 12_renamed into 13_final.
#
# Final set, three disjoint parts, all on Chr01–Chr20 with confidence >= medium
# and TE overlap <= 0.5:
#   high_confidence            ZH13 v3 identity >= 50; expression not required
#   expressed_zh13_below_50    expression supported and ZH13 v3 identity < 50
#   homolog_unexpressed        no expression support, ZH13 v3 identity < 50, and
#                              at least one soybean reference whose highest-bitscore
#                              DIAMOND hit has identity >= 50, query coverage >= 50%,
#                              subject coverage >= 50%, E-value <= 1e-5
#
# A reference counts only through its own highest-bitscore hit. A lower-scoring
# alignment on that reference is not used. ZH13 v2/v3 are two annotations of one
# cultivar, and the three Wm82 sets are versions of Williams 82, so five hits are
# not five independent confirmations.
#
# Before the merge, abnormal models are removed again so the step stays valid on
# an unfiltered 12_renamed. A transcript whose CDS length is not divisible by 3
# removes the whole gene. A transcript with no CDS, no ATG, an internal stop, or
# no terminal stop is removed; the gene is removed only when nothing valid remains.
# Gene IDs already listed as fully dropped in 12_renamed are rejected if they are
# still present. Locus IDs are not renumbered.
#
# After that set is written, short high-confidence models are revised in place.
# 13_final keeps one GFF and one FASTA set. The previous sequences of an
# overwritten gene are stored in PREFIX.final.corrections.tsv, and any ID change
# is stored in PREFIX.final.id_map.tsv. expressed_zh13_below_50 is not revised.
#
# resolve_final_redundancy.py then edits that same 13_final directory:
#   - merge gene IDs that share the same CDS and protein, keeping the ID with
#     the UTR transcripts
#   - restore an extension when two new same-strand CDS overlap by >= 200 bp
#     and their original CDS did not
#   - drop an unexpressed same-strand fragment whose CDS and protein sit inside
#     a longer gene of the same ZH13 locus
# antisense_qc.py keeps an opposite-strand CDS overlap only when both genes
# have reliable homology. Pairs that also overlap in ZH13, and pairs that do
# not, both stay. Single-exon models without reliable homology, and the other
# unresolved antisense pairs, are removed from this same 13_final set and
# recorded in PREFIX.final.antisense_excluded.tsv. The GFF is not given a
# separate gene type for them.
#
# Usage:
#   bash 10_final_anno.sh
#   SAMPLE=KD35 bash 10_final_anno.sh
#   RUN_BUSCO_FINAL=0 bash 10_final_anno.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${ROOT:-/tiandata2/zzh/sjy_assembly}"
FILTER_SAMPLE="${SAMPLE:-}"

FINAL_MIN_PIDENT="${FINAL_MIN_PIDENT:-50}"
FINAL_MIN_COV="${FINAL_MIN_COV:-0.5}"
FINAL_MAX_EVALUE="${FINAL_MAX_EVALUE:-1e-5}"
MINIPROT="${MINIPROT:-/home/zzh/anaconda3/envs/aligner/bin/miniprot}"
ZH13_FASTA="${ZH13_FASTA:-/home/zzh/ref/blastdb/prot/ZH13v3.0_mRNA_nuclear_20250729_pro.fasta}"
CURATE_THREADS="${CURATE_THREADS:-16}"
RUN_BUSCO_FINAL="${RUN_BUSCO_FINAL:-1}"
BUSCO_BIN="${BUSCO_BIN:-/home/zzh/anaconda3/envs/qc_busco/bin/busco}"
BUSCO_LINEAGE="${BUSCO_LINEAGE:-/tiandata2/zzh/forother/fc/busco_downloads/lineages/embryophyta_odb10}"

genome_for() {
    local prefix="$1"
    local primary="${ROOT}/07_polish_rename/${prefix}_primary/03_final/${prefix}_primary_Chr.fa"
    local batch2="${ROOT}/batch2/07_polish_rename/${prefix}_primary/03_final/${prefix}_primary_Chr.fa"
    if [[ -s "$primary" ]]; then
        echo "$primary"
    elif [[ -s "$batch2" ]]; then
        echo "$batch2"
    else
        die "missing chromosome FASTA for ${prefix}"
    fi
}

run_final_busco() {
    local down="$1"
    local prefix="$2"
    local busco_dir="${down}/04_busco"
    local longest="${busco_dir}/${prefix}.final.longest.faa"
    mkdir -p "$busco_dir"
    python3 - "$down" "$prefix" "${longest}.tmp" <<'PY'
import sys
down, prefix, dest = sys.argv[1:]
src = f"{down}/13_final/{prefix}.final.proteins.fa"
best = {}
name = None
seq = []
def keep():
    if not name:
        return
    gene = name.rsplit(".m", 1)[0]
    text = "".join(seq)
    old = best.get(gene)
    if old is None or len(text) > len(old[1]):
        best[gene] = (name, text)
with open(src) as handle:
    for line in handle:
        if line.startswith(">"):
            keep()
            name = line[1:].split()[0]
            seq = []
        else:
            seq.append(line.strip())
keep()
with open(dest, "w") as handle:
    for gene in best:
        tid, text = best[gene]
        handle.write(f">{tid}\n")
        for i in range(0, len(text), 60):
            handle.write(text[i:i+60] + "\n")
print(f"[INFO] {prefix} longest proteins {len(best)}", flush=True)
PY
    if [[ -s "$longest" ]] && cmp -s "${longest}.tmp" "$longest"; then
        rm -f "${longest}.tmp"
    else
        mv "${longest}.tmp" "$longest"
    fi
    publish_busco() {
        local summary="$1"
        local run_dir json
        [[ -n "$summary" && -s "$summary" ]] || return 1
        run_dir="$(dirname "$summary")/run_embryophyta_odb10"
        json="${summary%.txt}.json"
        cp -f "$summary" "${down}/13_final/${prefix}.final.busco.short_summary.txt"
        [[ -s "$json" ]] && cp -f "$json" "${down}/13_final/${prefix}.final.busco.short_summary.json"
        [[ -s "${run_dir}/full_table.tsv" ]] && cp -f "${run_dir}/full_table.tsv" "${down}/13_final/${prefix}.final.busco.full_table.tsv"
        [[ -s "${run_dir}/missing_busco_list.tsv" ]] && cp -f "${run_dir}/missing_busco_list.tsv" "${down}/13_final/${prefix}.final.busco.missing_busco_list.tsv"
        log "BUSCO written to ${down}/13_final"
        log "BUSCO $(grep -E '^[[:space:]]+C:' "$summary" | tr -s ' ')"
    }
    local summary
    summary="$(ls -1 "${busco_dir}/${prefix}_final_protein_busco"/short_summary.specific.*.txt 2>/dev/null | head -n 1 || true)"
    if [[ -n "$summary" && -s "$summary" && "$summary" -nt "$longest" ]]; then
        publish_busco "$summary"
        return
    fi
    [[ -x "$BUSCO_BIN" ]] || die "missing busco: ${BUSCO_BIN}"
    [[ -s "${BUSCO_LINEAGE}/dataset.cfg" ]] || die "missing BUSCO lineage: ${BUSCO_LINEAGE}"
    log "BUSCO proteins: ${prefix}"
    (
        cd "$busco_dir"
        export PATH="/home/zzh/anaconda3/envs/qc_busco/bin:${PATH}"
        "$BUSCO_BIN" \
            -m proteins \
            -i "$longest" \
            -l "$BUSCO_LINEAGE" \
            --offline \
            -o "${prefix}_final_protein_busco" \
            -c "$CURATE_THREADS" \
            --force
    )
    summary="$(ls -1 "${busco_dir}/${prefix}_final_protein_busco"/short_summary.specific.*.txt 2>/dev/null | head -n 1 || true)"
    publish_busco "$summary" || die "BUSCO summary missing for ${prefix}"
}

REF_PROT_FASTAS=(
    "/home/zzh/ref/blastdb/prot/ZH13.v2.prot.fasta"
    "/home/zzh/ref/blastdb/prot/ZH13v3.0_mRNA_nuclear_20250729_pro.fasta"
    "/home/zzh/ref/blastdb/prot/Wm82.a2.v1.prot.fasta"
    "/home/zzh/ref/blastdb/prot/glyma.Wm82.gnm4.ann1.T8TQ.protein.faa"
    "/home/zzh/ref/blastdb/prot/glyma.Wm82.gnm6.ann1.PKSW.protein.faa"
)

DOWNSTREAM_DIRS=(
    "${ROOT}/annot/KD35_annotation_downstream"
    "${ROOT}/annot/N2L_annotation_downstream"
    "${ROOT}/batch2/annot/BZ_annotation_downstream"
    "${ROOT}/batch2/annot/DY11_annotation_downstream"
    "${ROOT}/batch2/annot/TH2_annotation_downstream"
    "${ROOT}/batch2/annot/TZX1998_annotation_downstream"
    "${ROOT}/batch2/annot/TZX253_annotation_downstream"
)

timestamp() { date "+%F %T"; }
log() { echo "[$(timestamp)] $*"; }
die() { echo "[$(timestamp)] ERROR: $*" >&2; exit 1; }

export FINAL_MIN_PIDENT FINAL_MAX_EVALUE FINAL_MIN_COV
export FINAL_REFS="$(IFS=:; echo "${REF_PROT_FASTAS[*]}")"

for down in "${DOWNSTREAM_DIRS[@]}"; do
    base="$(basename "$down")"
    prefix="${base%%_annotation_downstream}"
    if [[ -n "$FILTER_SAMPLE" && "$prefix" != "$FILTER_SAMPLE" ]]; then
        continue
    fi
    [[ -d "$down/12_renamed" ]] || die "missing 12_renamed: $down"
    log "final annotation: ${prefix}"
    FINAL_DOWNSTREAM="$down" FINAL_PREFIX="$prefix" python3 - <<'PY'
import csv, os, re, shutil, sys
from collections import defaultdict

down = os.environ["FINAL_DOWNSTREAM"]
prefix = os.environ["FINAL_PREFIX"]
pmin = float(os.environ["FINAL_MIN_PIDENT"])
cmin = float(os.environ["FINAL_MIN_COV"])
emax = float(os.environ["FINAL_MAX_EVALUE"])
zh13_split = 50.0
ref_fastas = [p for p in os.environ["FINAL_REFS"].split(":") if p]

renamed = os.path.join(down, "12_renamed")
stage = os.path.join(down, "13_final.partial")
final_dir = os.path.join(down, "13_final")
gff_path = os.path.join(renamed, f"{prefix}.renamed.gff3")
cds_path = os.path.join(renamed, f"{prefix}.renamed.cds.fa")
tx_path = os.path.join(renamed, f"{prefix}.renamed.transcripts.fa")
pep_path = os.path.join(renamed, f"{prefix}.renamed.proteins.fa")
summary_path = os.path.join(renamed, f"{prefix}.renamed.annotation.summary.tsv")
map_path = os.path.join(renamed, f"{prefix}.gene_id_mapping.tsv")
query_pep = os.path.join(down, "02_sequences", f"{prefix}.protein.fa")
drop_tsv = os.path.join(renamed, f"{prefix}.dropped.tsv")
drop_list = os.path.join(renamed, f"{prefix}.dropped_cds_nonmod3.genes.list")
zkey = "ZH13v3.0_mRNA_nuclear_20250729_pro_identity"
tx_suffix = re.compile(r"\.[tm]\d+$")
nuc_re = re.compile(r"^Chr(0[1-9]|1[0-9]|20)$")
stop_codons = {"TAA", "TAG", "TGA"}
conf_ok = {"high_confidence", "medium_confidence"}

def fail(msg):
    print(f"[ERROR] {prefix}: {msg}", file=sys.stderr)
    shutil.rmtree(stage, ignore_errors=True)
    sys.exit(1)

def need(path):
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        fail(f"missing or empty: {path}")

for path in (gff_path, cds_path, tx_path, pep_path, summary_path, map_path, query_pep):
    need(path)
for path in ref_fastas:
    need(path)

def parse_attr(text):
    out = {}
    for item in text.strip().split(";"):
        if "=" in item:
            key, val = item.split("=", 1)
            out[key] = val
    return out

def fasta_seqs(path):
    seqs = {}
    name = None
    chunks = []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            seqs[name] = "".join(chunks)
    return seqs

def fasta_lens(path):
    lens = {}
    name = None
    n = 0
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    lens[name] = n
                name = line[1:].split()[0]
                n = 0
            else:
                n += len(line.strip())
        if name is not None:
            lens[name] = n
    return lens

def fnum(value):
    if value in (None, "", "NA", "na", "."):
        return None
    try:
        return float(value)
    except ValueError:
        return None

def gene_of_tx(tid):
    return tx_suffix.sub("", tid)

def old_gene_of_query(query, old_to_new):
    if not query.startswith("evm.model."):
        return None
    rest = query[len("evm.model."):]
    gene = "evm.TU." + rest
    if gene in old_to_new:
        return gene
    if "." in rest:
        shorter = "evm.TU." + rest.rsplit(".", 1)[0]
        if shorter in old_to_new:
            return shorter
    return None

def ref_label(path):
    base = os.path.basename(path)
    for suffix in (".gz", ".fasta", ".faa", ".fa", ".pep"):
        if base.endswith(suffix):
            base = base[: -len(suffix)]
    return re.sub(r"[^A-Za-z0-9._-]", "_", base)

def classify_cds(seq):
    seq = (seq or "").upper().replace("U", "T")
    if not seq:
        return ["无CDS"]
    reasons = []
    if len(seq) % 3 != 0:
        reasons.append("CDS长度不能被3整除")
        return reasons
    if not seq.startswith("ATG"):
        reasons.append("缺少起始密码子ATG")
    codons = [seq[i:i + 3] for i in range(0, len(seq), 3)]
    if any(codon in stop_codons for codon in codons[:-1]):
        reasons.append("内部终止密码子")
    if not codons or codons[-1] not in stop_codons:
        reasons.append("缺少终止密码子")
    return reasons

def load_prior_drops():
    banned = set()
    if os.path.isfile(drop_tsv) and os.path.getsize(drop_tsv) > 0:
        with open(drop_tsv) as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row.get("删除范围") == "基因" and row.get("基因ID"):
                    banned.add(row["基因ID"].strip())
    if os.path.isfile(drop_list) and os.path.getsize(drop_list) > 0:
        with open(drop_list) as handle:
            for line in handle:
                gene = line.strip()
                if gene:
                    banned.add(gene)
    return banned

# ----- GFF and CDS -----
cds_seq = fasta_seqs(cds_path)
comments = []
order = []
groups = defaultdict(list)
mrna_gene = {}
mrna_has_cds = defaultdict(bool)
with open(gff_path) as handle:
    for line in handle:
        if line.startswith("#") or not line.strip():
            if not order:
                comments.append(line)
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue
        attr = parse_attr(parts[8])
        feature = parts[2]
        if feature == "gene":
            gid = attr.get("ID", "")
            if not gid:
                continue
            order.append(gid)
            groups[gid].append(parts)
        elif feature in ("mRNA", "transcript"):
            tid = attr.get("ID", "")
            gid = attr.get("Parent", "").split(",")[0]
            if not tid or not gid:
                continue
            mrna_gene[tid] = gid
            groups[gid].append(parts)
        else:
            parent = attr.get("Parent", "").split(",")[0]
            gid = mrna_gene.get(parent)
            if gid is None:
                continue
            if feature == "CDS":
                mrna_has_cds[parent] = True
            groups[gid].append(parts)

# ----- ORF filter, including models already removed from this 12_renamed -----
prior = load_prior_drops()
drop_rows = []
kept_tx = {}
orf_genes = set()
prior_blocked = 0
for gid in order:
    if gid in prior:
        prior_blocked += 1
        for parts in groups[gid]:
            if parts[2] in ("mRNA", "transcript"):
                tid = parse_attr(parts[8]).get("ID", "")
                drop_rows.append(("基因", gid, tid, "已在12_renamed删除清单", "脚本复核"))
        continue
    tx_reasons = []
    for parts in groups[gid]:
        if parts[2] not in ("mRNA", "transcript"):
            continue
        tid = parse_attr(parts[8]).get("ID", "")
        if not mrna_has_cds[tid]:
            reasons = ["无CDS"]
        else:
            reasons = classify_cds(cds_seq.get(tid, ""))
        tx_reasons.append((tid, reasons))
    if any("CDS长度不能被3整除" in reasons for _, reasons in tx_reasons):
        orf_genes.add(gid)
        for tid, reasons in tx_reasons:
            why = ";".join(reasons) if reasons else "同基因存在CDS长度不能被3整除的转录本"
            drop_rows.append(("基因", gid, tid, why, "脚本复核"))
        continue
    good = [tid for tid, reasons in tx_reasons if not reasons]
    bad = [(tid, reasons) for tid, reasons in tx_reasons if reasons]
    if not good:
        orf_genes.add(gid)
        for tid, reasons in bad:
            drop_rows.append(("基因", gid, tid, ";".join(reasons), "脚本复核"))
        continue
    kept_tx[gid] = set(good)
    for tid, reasons in bad:
        drop_rows.append(("转录本", gid, tid, ";".join(reasons), "脚本复核"))

# ----- summary classes -----
summary = {}
with open(summary_path) as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        summary[row["gene_id"]] = row

classes = {}
candidates = []
for gid in kept_tx:
    row = summary.get(gid)
    if row is None or not nuc_re.match(row.get("chrom", "")):
        continue
    te = fnum(row.get("TE_overlap_ratio"))
    if te is None:
        te = 0.0
    if te > 0.5 or row.get("confidence_level") not in conf_ok:
        continue
    ident = fnum(row.get(zkey))
    ident = 0.0 if ident is None else ident
    expressed = row.get("expression_supported") == "yes"
    if ident >= zh13_split:
        classes[gid] = "high_confidence"
    elif expressed:
        classes[gid] = "expressed_zh13_below_50"
    else:
        candidates.append(gid)

# ----- homolog evidence: best bitscore hit of each reference -----
old_to_new = {}
with open(map_path) as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        old_to_new[row["old_gene_id"]] = row["new_gene_id"]
want_old = {old for old, new in old_to_new.items() if new in set(candidates)}
qlen = {}
current = None
length = 0
keep = False
with open(query_pep) as handle:
    for line in handle:
        if line.startswith(">"):
            if current is not None and keep:
                qlen[current] = length
            current = line[1:].split()[0]
            length = 0
            old = old_gene_of_query(current, old_to_new)
            keep = old in want_old
        elif keep:
            length += len(line.strip())
    if current is not None and keep:
        qlen[current] = length

subject_len = {}
for fasta in ref_fastas:
    subject_len.update(fasta_lens(fasta))

best_hit = defaultdict(dict)
for fasta in ref_fastas:
    label = ref_label(fasta)
    blast = os.path.join(down, "05_diamond", f"{prefix}.vs_{label}.tsv")
    need(blast)
    with open(blast) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            old = old_gene_of_query(parts[0], old_to_new)
            if old not in want_old:
                continue
            gid = old_to_new[old]
            query_length = qlen.get(parts[0])
            if not query_length:
                continue
            try:
                ident = float(parts[2])
                qstart, qend = int(parts[6]), int(parts[7])
                sstart, send = int(parts[8]), int(parts[9])
                evalue = float(parts[10])
                bits = float(parts[11])
            except ValueError:
                continue
            if evalue > emax:
                continue
            slen = subject_len.get(parts[1])
            scov = None if not slen else (abs(send - sstart) + 1) / slen
            rec = (bits, ident, (abs(qend - qstart) + 1) / query_length, scov, evalue, parts[1])
            prev = best_hit[gid].get(label)
            if prev is None or rec[0] > prev[0]:
                best_hit[gid][label] = rec

support = {}
for gid in candidates:
    chosen = None
    for label, rec in best_hit.get(gid, {}).items():
        bits, ident, qcov, scov, evalue, hit = rec
        if scov is None or ident < pmin or qcov < cmin or scov < cmin:
            continue
        if chosen is None or bits > chosen[0]:
            chosen = (bits, label, hit, ident, qcov, scov, evalue)
    if chosen is not None:
        classes[gid] = "homolog_unexpressed"
        support[gid] = chosen

selected = [gid for gid in order if gid in classes]
if not selected:
    fail("no genes selected")

# ----- write, then check the files agree -----
if os.path.isdir(stage):
    shutil.rmtree(stage)
os.makedirs(stage)

def na(value, digits=None):
    if value is None:
        return "NA"
    if digits is None:
        return str(value)
    return f"{value:.{digits}f}"

class_path = os.path.join(stage, f"{prefix}.final.gene_class.tsv")
with open(class_path, "w") as handle:
    handle.write("\t".join([
        "gene_id", "class", "chrom", "start", "end", "strand",
        "zh13_v3_identity", "expression_supported",
        "support_ref", "support_hit", "support_pident", "support_qcov",
        "support_scov", "support_evalue", "support_bitscore",
    ]) + "\n")
    for gid in selected:
        row = summary[gid]
        sup = support.get(gid)
        if sup is None:
            fields = ["NA"] * 7
        else:
            bits, label, hit, ident, qcov, scov, evalue = sup
            fields = [label, hit, na(ident, 1), na(qcov, 4), na(scov, 4), f"{evalue:.3g}", na(bits, 1)]
        ident = fnum(row.get(zkey))
        handle.write("\t".join([
            gid, classes[gid], row["chrom"], row["start"], row["end"], row["strand"],
            "0" if ident is None else na(ident, 1), row.get("expression_supported", "NA"),
            *fields,
        ]) + "\n")

with open(os.path.join(stage, f"{prefix}.final.genes.list"), "w") as handle:
    for gid in selected:
        handle.write(gid + "\n")

with open(os.path.join(stage, f"{prefix}.final.orf_dropped.tsv"), "w") as handle:
    handle.write("删除范围\t基因ID\t转录本ID\t原因\t批次\n")
    for row in drop_rows:
        handle.write("\t".join(row) + "\n")

criteria = os.path.join(stage, f"{prefix}.final.criteria.txt")
with open(criteria, "w") as handle:
    handle.write(
        "最终蛋白编码注释由三类互不重叠的核基因组成。\n"
        "共同条件：Chr01–Chr20，置信度至少 medium，转座元件重叠不超过 0.5。\n"
        "high_confidence：相对 ZH13 v3 的蛋白相似度 >= "
        f"{zh13_split:g}，不要求表达。\n"
        "expressed_zh13_below_50：有表达支持，相对 ZH13 v3 的相似度 < "
        f"{zh13_split:g}（没有命中记为 0）。这一类不是品系特异基因的最终判定。\n"
        "homolog_unexpressed：没有表达支持，相对 ZH13 v3 的相似度 < "
        f"{zh13_split:g}，并且 ZH13 v2、ZH13 v3、Wm82 a2.v1、Wm82 gnm4、Wm82 gnm6 "
        "中至少一套的最高 bitscore 命中同时满足："
        f"相似度 >= {pmin:g}，查询蛋白覆盖 >= {cmin:g}，参考蛋白覆盖 >= {cmin:g}，"
        f"E-value <= {emax:g}。\n"
        "每个参考只用自己的最高 bitscore 命中，不用该参考上得分更低的比对。\n"
        "五套蛋白是中黄13和 Williams 82 的不同注释版本，不是五次独立验证。\n"
        "写入前复核 ORF：任一转录本 CDS 长度不能被 3 整除则删除整个基因；"
        "无 CDS、缺少 ATG、内部终止或缺少终止密码子的转录本被删除，"
        "仅当该基因不再剩有效转录本时才删除基因。位点编号不重排。\n"
        "12_renamed 删除清单中的基因如果仍出现在输入里，会被拒绝。\n"
    )

out_gff = os.path.join(stage, f"{prefix}.final.gff3")
gff_tx = []
with open(out_gff, "w") as handle:
    handle.writelines(comments)
    if not comments:
        handle.write("##gff-version 3\n")
    for gid in selected:
        body = []
        starts = []
        ends = []
        for parts in groups[gid]:
            feature = parts[2]
            attr = parse_attr(parts[8])
            if feature == "gene":
                continue
            if feature in ("mRNA", "transcript"):
                tid = attr.get("ID", "")
                if tid not in kept_tx[gid]:
                    continue
                starts.append(int(parts[3]))
                ends.append(int(parts[4]))
                body.append(parts)
                gff_tx.append(tid)
            else:
                parent = attr.get("Parent", "").split(",")[0]
                if parent in kept_tx[gid]:
                    body.append(parts)
        gene_line = next(parts for parts in groups[gid] if parts[2] == "gene")
        if starts:
            gene_line = gene_line[:]
            gene_line[3] = str(min(starts))
            gene_line[4] = str(max(ends))
        handle.write("\t".join(gene_line) + "\n")
        for parts in body:
            handle.write("\t".join(parts) + "\n")

def write_fasta(src, dst, wanted):
    found = []
    with open(src) as fin, open(dst, "w") as fout:
        keep = False
        for line in fin:
            if line.startswith(">"):
                sid = line[1:].split()[0]
                keep = sid in wanted
                if keep:
                    found.append(sid)
                    fout.write(line if line.endswith("\n") else line + "\n")
            elif keep:
                fout.write(line if line.endswith("\n") else line + "\n")
    return found

wanted = set(gff_tx)
cds_ids = write_fasta(cds_path, os.path.join(stage, f"{prefix}.final.cds.fa"), wanted)
tx_ids = write_fasta(tx_path, os.path.join(stage, f"{prefix}.final.transcripts.fa"), wanted)
pep_ids = write_fasta(pep_path, os.path.join(stage, f"{prefix}.final.proteins.fa"), wanted)

if set(cds_ids) != wanted or set(tx_ids) != wanted or set(pep_ids) != wanted:
    missing = {
        "cds": sorted(wanted - set(cds_ids))[:5],
        "tx": sorted(wanted - set(tx_ids))[:5],
        "pep": sorted(wanted - set(pep_ids))[:5],
    }
    fail(f"FASTA transcripts do not match GFF: {missing}")

out_cds = fasta_seqs(os.path.join(stage, f"{prefix}.final.cds.fa"))
bad_out = []
for tid, seq in out_cds.items():
    reasons = classify_cds(seq)
    if reasons:
        bad_out.append((tid, ";".join(reasons)))
        if len(bad_out) >= 5:
            break
if bad_out:
    fail(f"output CDS still abnormal: {bad_out}")

written_genes = []
with open(out_gff) as handle:
    for line in handle:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) >= 9 and parts[2] == "gene":
            written_genes.append(parse_attr(parts[8]).get("ID", ""))
if written_genes != selected:
    fail("GFF gene list does not match the class table")
if any(gid in prior or gid in orf_genes for gid in selected):
    fail("a dropped gene entered the final set")

counts = defaultdict(int)
for gid in selected:
    counts[classes[gid]] += 1
n_tx_only = sum(1 for scope, *_ in drop_rows if scope == "转录本")
print(
    f"{prefix}\t{counts['high_confidence']}\t{counts['expressed_zh13_below_50']}\t"
    f"{counts['homolog_unexpressed']}\t{len(selected)}\t{len(gff_tx)}\t"
    f"orf_genes={len(orf_genes)}\torf_tx={n_tx_only}\tprior_blocked={prior_blocked}",
    flush=True,
)

if os.path.isdir(final_dir):
    shutil.rmtree(final_dir)
os.rename(stage, final_dir)
PY
    genome="$(genome_for "$prefix")"
    [[ -x "$MINIPROT" ]] || die "miniprot not executable: ${MINIPROT}"
    log "revise short models in 13_final: ${prefix}"
    python3 "${SCRIPT_DIR}/curate_final_models.py" \
        --downstream "$down" \
        --prefix "$prefix" \
        --genome "$genome" \
        --zh13-fasta "$ZH13_FASTA" \
        --miniprot "$MINIPROT" \
        --threads "$CURATE_THREADS"
    log "resolve redundant models in 13_final: ${prefix}"
    python3 "${SCRIPT_DIR}/resolve_final_redundancy.py" \
        --downstream "$down" \
        --prefix "$prefix"
    log "write antisense QC tables in 13_final: ${prefix}"
    python3 "${SCRIPT_DIR}/antisense_qc.py" \
        --downstream "$down" \
        --prefix "$prefix"
    # BUSCO stays after ORF filtering, short-model repair, same-strand cleanup,
    # and antisense exclusion. It reads the proteins left in this 13_final.
    if [[ "$RUN_BUSCO_FINAL" -eq 1 ]]; then
        run_final_busco "$down" "$prefix"
    fi
done

log "done"
