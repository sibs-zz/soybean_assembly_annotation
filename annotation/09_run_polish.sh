#!/usr/bin/env bash
set -euo pipefail

############################################################
# 09  功能注释收尾（原 04/07_run_polish.sh）
#
# 上游（短读对标 ZH13 v3，无 Iso-seq）：
#   04 transcripts → 05 Mikado → 06 PASA align → 07 EVM → 08 PASA UTR
#   输入 GFF：${sp}_models/${sp}.models.gff3
#
# 其它输入：
#   基因组  07_polish_rename/${SAMPLE}/03_final/${SAMPLE}_Chr.fa
#           回退 ${sp}_TE/${sp}_Chr.fa.masked
#   RM.out  ${sp}_TE/${sp}_Chr.fa.out
#   RNA     rna_clean/${sp}_RNA_Atlas_{1,2}.fq.gz
#
# Python 辅助脚本（同目录）：
#   filter_dot_protein_transcripts.py
#   correct_nterm_starts.py
#   build_publication_annotation.py
#   rename_gene_loci.py
#   pack_structural_gff.py
#   extract_lineage_specific.py
#   extract_hc_sequences.py
#
# 用法：
#   bash 09_run_polish.sh
#   SAMPLE=N2L_primary bash 09_run_polish.sh
#   REQUIRE_MODELS=0 bash 09_run_polish.sh   # 无 08 时回退 BRAKER 原 GFF
############################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ======================== 配置区 =============================

SAMPLES=(
    "KD35_primary"
    "N2L_primary"
)

ROOT="${ROOT:-/tiandata2/zzh/sjy_assembly}"
RNA_DIR="${RNA_DIR:-${ROOT}/rna_clean}"

THREADS="${THREADS:-64}"
SKIP_IF_DONE="${SKIP_IF_DONE:-1}"
CLEAN_SHARED_WORKDIRS="${CLEAN_SHARED_WORKDIRS:-1}"

# 若设置则只跑该样本（需在 SAMPLES 中）
FILTER_SAMPLE="${SAMPLE:-}"
# 1=必须先有 08 的 models.gff3；0=回退 BRAKER braker.gff3
REQUIRE_MODELS="${REQUIRE_MODELS:-1}"
PUB_MIN_AA="${PUB_MIN_AA:-29}"

############################
# Conda env bin paths
############################

ENV_AGAT="/home/zzh/anaconda3/envs/EDTA/bin"
ENV_ANNP="/home/zzh/anaconda3/envs/annp/bin"
ENV_BUSCO="/home/zzh/anaconda3/envs/qc_busco/bin"
ENV_EGGNOG="/home/zzh/anaconda3/envs/eggnog213/bin"

############################
# BUSCO configuration
############################

RUN_BUSCO=1
BUSCO_LINEAGE="/tiandata2/zzh/forother/fc/busco_downloads/lineages/embryophyta_odb10"
BUSCO_MODE="proteins"

############################
# SwissProt / DIAMOND configuration
############################

RUN_DIAMOND=1
SWISSPROT_DMND="/tiandata2/zzh/tools/uniprot_sprot.dmnd"

############################
# Reference proteome homology (fixed 5 soybean refs)
############################

RUN_REF_PROTEOME=1
REF_PROT_DIR="/home/zzh/ref/blastdb"
REF_PROT_SEARCH="${REF_PROT_DIR}/prot"
REF_DMND_DIR="/tiandata2/zzh/tools/ref_diamond"
REF_DIAMOND_EVALUE="1e-5"
REF_DIAMOND_MAX_TARGET_SEQS=5
# 写死：仅这 5 套（ZH13 v2/v3 + Wm82 a2.v1/gnm4/gnm6）
REF_PROT_FASTAS=(
    "${REF_PROT_SEARCH}/ZH13.v2.prot.fasta"
    "${REF_PROT_SEARCH}/ZH13v3.0_mRNA_nuclear_20250729_pro.fasta"
    "${REF_PROT_SEARCH}/Wm82.a2.v1.prot.fasta"
    "${REF_PROT_SEARCH}/glyma.Wm82.gnm4.ann1.T8TQ.protein.faa"
    "${REF_PROT_SEARCH}/glyma.Wm82.gnm6.ann1.PKSW.protein.faa"
)

############################
# Ortholog-guided CDS start correction (before homolog annotation)
############################

RUN_START_FIX=0
START_FIX_ZH13_FA="/home/zzh/ref/ZH13_v3/ZH13v3.0_mRNA_nuclear_20250729_pro.fasta"
START_FIX_WM82_FA="/home/zzh/ref/glyma.Wm82.gnm6.ann1.PKSW.protein.faa"
START_FIX_MIN_PIDENT=95
START_FIX_MIN_NTERM_IDENT=0.95
START_FIX_MIN_MISSING_AA=15
START_FIX_SCRIPT="${SCRIPT_DIR}/correct_nterm_starts.py"

############################
# TE overlap configuration
############################

RUN_TE_OVERLAP=1
TE_SUPERFAMILIES="DNA,LTR,LINE,MITE,RC,SINE,Penelope"

############################
# eggNOG configuration
############################

RUN_EGGNOG=1
EMAPPER_PY="/tiandata2/zzh/tools/eggnog-mapper-2.1.13/emapper.py"
EGGNOG_DATA_DIR="/tiandata2/zzh/tools/eggnog-mapper-2.1.13/data"
export EGGNOG_DATA_DIR

############################
# InterProScan configuration
############################

RUN_INTERPRO=1
INTERPROSCAN_BIN="/home/zzh/software/interproscan/interproscan-5.77-108.0/interproscan.sh"

############################
# Salmon RNA-seq expression support
############################

RUN_SALMON=1
SALMON_LIBTYPE="A"
TPM_CUTOFF=1

############################
# Publication-ready outputs
############################

PUB_MIN_CONFIDENCE="medium_confidence"
PUB_MAX_TE_RATIO="0.5"
PUB_REQUIRE_EXPRESSION=0
# HC 对齐已发表大豆核基因数（约 5.4–6.0 万）：只要 Chr01–20，且 vs ZH13 v3 蛋白 identity ≥ 50
PUB_NUCLEAR_ONLY="${PUB_NUCLEAR_ONLY:-1}"
PUB_MIN_ZH13_PIDENT="${PUB_MIN_ZH13_PIDENT:-50}"

############################
# Locus renaming
############################

RUN_GENE_RENAME=1
GENE_ID_DIGITS=6
GENE_ID_STEP=100

# ============================================================

timestamp() { date "+%F %T"; }
log()       { echo "[$(timestamp)] $*"; }
die()       { echo "[$(timestamp)] ERROR: $*" >&2; exit 1; }

# 按样本初始化路径（PREFIX = sp，如 N2L_primary → N2L）
init_sample_paths() {
    local sample="$1"
    PREFIX="${sample%%_*}"

    GENOME="${ROOT}/07_polish_rename/${sample}/03_final/${sample}_Chr.fa"
    if [[ ! -s "$GENOME" ]]; then
        GENOME="${SCRIPT_DIR}/${PREFIX}_TE/${PREFIX}_Chr.fa.masked"
    fi

    BRAKER_DIR="${SCRIPT_DIR}/BRAKER_${PREFIX}"
    MODELS_GFF="${SCRIPT_DIR}/${PREFIX}_models/${PREFIX}.models.gff3"
    BRAKER_RAW="${BRAKER_DIR}/braker.gff3"
    if [[ -s "$MODELS_GFF" ]]; then
        BRAKER_GFF="$MODELS_GFF"
    else
        BRAKER_GFF="$BRAKER_RAW"
    fi
    OUTDIR="${SCRIPT_DIR}/${PREFIX}_annotation_downstream"

    RM_OUT="${SCRIPT_DIR}/${PREFIX}_TE/${PREFIX}_Chr.fa.out"

    RNA_R1="${RNA_DIR}/${PREFIX}_RNA_Atlas_1.fq.gz"
    RNA_R2="${RNA_DIR}/${PREFIX}_RNA_Atlas_2.fq.gz"
    RNA_SAMPLE="${PREFIX}_RNA_Atlas"

    PUB_DIR="${OUTDIR}/11_publication"
    ANNOTATED_GFF="${PUB_DIR}/${PREFIX}.annotated.gff3"
    GENE_FUNC_TABLE="${PUB_DIR}/${PREFIX}.gene_functional_annotation.tsv"
    HC_GENE_LIST="${PUB_DIR}/${PREFIX}.high_confidence.genes.list"
    HC_ANNOTATED_GFF="${PUB_DIR}/${PREFIX}.high_confidence.annotated.gff3"
    HC_PROTEIN_FA="${PUB_DIR}/${PREFIX}.high_confidence.proteins.fa"

    GENE_ID_PREFIX="${PREFIX}"

    RENAME_DIR="${OUTDIR}/12_renamed"
    RENAMED_GFF="${RENAME_DIR}/${PREFIX}.renamed.gff3"
    RENAMED_GENE_TABLE="${RENAME_DIR}/${PREFIX}.renamed.gene_functional_annotation.tsv"
    RENAMED_SUMMARY="${RENAME_DIR}/${PREFIX}.renamed.annotation.summary.tsv"
    RENAMED_TRANSCRIPT_FA="${RENAME_DIR}/${PREFIX}.renamed.transcripts.fa"
    RENAMED_CDS_FA="${RENAME_DIR}/${PREFIX}.renamed.cds.fa"
    RENAMED_PROTEIN_FA="${RENAME_DIR}/${PREFIX}.renamed.proteins.fa"
    RENAMED_HC_GFF="${RENAME_DIR}/${PREFIX}.renamed.high_confidence.annotated.gff3"
    RENAMED_HC_PROTEIN_FA="${RENAME_DIR}/${PREFIX}.renamed.high_confidence.proteins.fa"
    RENAMED_HC_TRANSCRIPT_FA="${RENAME_DIR}/${PREFIX}.renamed.high_confidence.transcripts.fa"
    RENAMED_HC_CDS_FA="${RENAME_DIR}/${PREFIX}.renamed.high_confidence.cds.fa"
    RENAMED_HC_GENE_LIST="${RENAME_DIR}/${PREFIX}.renamed.high_confidence.genes.list"
    GENE_ID_MAPPING="${RENAME_DIR}/${PREFIX}.gene_id_mapping.tsv"
    RENAMED_NUCLEAR_GFF="${RENAME_DIR}/${PREFIX}.renamed.nuclear.gff3"
    RENAMED_NUCLEAR_PROTEIN_FA="${RENAME_DIR}/${PREFIX}.renamed.nuclear.proteins.fa"

    LOG="${OUTDIR}/00_logs/${PREFIX}.pipeline.log"
}

check_file() {
    local f="$1"
    local name="$2"
    if [[ ! -s "${f}" ]]; then
        echo "[ERROR] ${name} not found or empty: ${f}"
        exit 1
    fi
}

check_dir() {
    local d="$1"
    local name="$2"
    if [[ ! -d "${d}" ]]; then
        echo "[ERROR] ${name} directory not found: ${d}"
        exit 1
    fi
}

outputs_ready() {
    local f
    [[ "${SKIP_IF_DONE}" -eq 1 ]] || return 1
    for f in "$@"; do
        [[ -s "${f}" ]] || return 1
    done
    return 0
}

outputs_ready_glob() {
    local pattern f
    [[ "${SKIP_IF_DONE}" -eq 1 ]] || return 1
    shopt -s nullglob
    local files=(${1})
    shopt -u nullglob
    for f in "${files[@]}"; do
        [[ -s "${f}" ]] && return 0
    done
    return 1
}

log_skip() {
    echo "[SKIP] $*"
}

run_cmd() {
    echo
    echo "[CMD] $*"
    "$@"
}

with_path() {
    local bin_dir="$1"
    shift
    PATH="${bin_dir}:$PATH" "$@"
}

# ============================================================
# 全局依赖检查（一次）
# ============================================================

check_dir "${ENV_AGAT}" "AGAT env bin"
check_dir "${ENV_ANNP}" "annp env bin"
check_dir "${ENV_BUSCO}" "BUSCO env bin"
check_dir "${ENV_EGGNOG}" "eggNOG env bin"

echo "Tool check:"
echo "AGAT convert: ${ENV_AGAT}/agat_convert_sp_gxf2gxf.pl"
echo "AGAT stats  : ${ENV_AGAT}/agat_sp_statistics.pl"
echo "gffread     : ${ENV_ANNP}/gffread"
echo "diamond     : ${ENV_ANNP}/diamond"
echo "salmon      : ${ENV_ANNP}/salmon"
echo "busco       : ${ENV_BUSCO}/busco"
echo "emapper.py  : ${ENV_EGGNOG}/emapper.py"

[[ -x "${ENV_AGAT}/agat_convert_sp_gxf2gxf.pl" ]] || echo "[WARN] AGAT convert not executable or missing."
[[ -x "${ENV_AGAT}/agat_sp_statistics.pl" ]] || echo "[WARN] AGAT statistics not executable or missing."
[[ -x "${ENV_ANNP}/gffread" ]] || { echo "[ERROR] gffread missing in annp env."; exit 1; }
[[ -x "${ENV_ANNP}/diamond" ]] || echo "[WARN] diamond missing in annp env."
[[ -x "${ENV_ANNP}/bedtools" ]] || { echo "[ERROR] bedtools missing in annp env (required for TE overlap)."; exit 1; }
[[ -x "${ENV_ANNP}/salmon" ]] || echo "[WARN] salmon missing in annp env."
[[ -x "${ENV_BUSCO}/busco" ]] || echo "[WARN] BUSCO missing in qc_busco env."
[[ -x "${ENV_EGGNOG}/emapper.py" ]] || echo "[WARN] emapper.py missing in eggnog213 env."
[[ -x "${EMAPPER_PY}" ]] || echo "[WARN] emapper.py missing: ${EMAPPER_PY}"
[[ -d "${EGGNOG_DATA_DIR}" ]] || echo "[WARN] eggNOG data dir missing: ${EGGNOG_DATA_DIR}"
if [[ "${RUN_INTERPRO}" -eq 1 ]]; then
    [[ -x "${INTERPROSCAN_BIN}" ]] || echo "[ERROR] InterProScan not executable: ${INTERPROSCAN_BIN}"
fi

if [[ "${CLEAN_SHARED_WORKDIRS}" -eq 1 ]]; then
    echo
    echo "============================================================"
    echo "Preflight: clean shared work directories"
    echo "============================================================"
    ANNO_TEMP="${SCRIPT_DIR}/temp"
    ANNO_AGAT_LOG="${SCRIPT_DIR}/agat_log_braker"
    if [[ -d "${ANNO_TEMP}" ]]; then
        echo "[INFO] Removing InterProScan temp: ${ANNO_TEMP}"
        rm -rf "${ANNO_TEMP}"
    fi
    mkdir -p "${ANNO_TEMP}"
    if [[ -d "${ANNO_AGAT_LOG}" ]]; then
        echo "[INFO] Removing AGAT log: ${ANNO_AGAT_LOG}"
        rm -rf "${ANNO_AGAT_LOG}"
    fi
    echo "[INFO] Preflight cleanup done."
fi

# ============================================================
# 按样本运行
# ============================================================

for SAMPLE in "${SAMPLES[@]}"; do
    [[ -z "$FILTER_SAMPLE" || "$SAMPLE" == "$FILTER_SAMPLE" ]] || continue

    init_sample_paths "$SAMPLE"

    if [[ "$REQUIRE_MODELS" == "1" && ! -s "$MODELS_GFF" ]]; then
        log "跳过 $SAMPLE：缺少 08 产物 $MODELS_GFF（先跑 04–08，或 REQUIRE_MODELS=0 回退 BRAKER GFF）"
        continue
    fi
    if [[ ! -s "$BRAKER_GFF" ]]; then
        log "跳过 $SAMPLE：缺少注释 GFF $BRAKER_GFF"
        continue
    fi

    mkdir -p "${OUTDIR}"/{00_logs,01_gff,02_sequences,03_stats,04_busco,05_diamond,06_eggnog,07_interpro,08_TE_overlap,09_salmon,10_summary,11_publication,12_renamed,tmp,02_sequences/start_fix}

    {
    echo "============================================================"
    echo "BRAKER downstream annotation pipeline"
    echo "Start time: $(date)"
    echo "SAMPLE : ${SAMPLE}"
    echo "PREFIX : ${PREFIX}"
    echo "GENOME : ${GENOME}"
    echo "BRAKER : ${BRAKER_DIR}"
    echo "OUTDIR : ${OUTDIR}"
    echo "SKIP_IF_DONE : ${SKIP_IF_DONE}"
    echo "REQUIRE_MODELS : ${REQUIRE_MODELS}"
    echo "GFF_IN : ${BRAKER_GFF}"
    echo "RUN_GENE_RENAME : ${RUN_GENE_RENAME}"
    echo "GENE_ID_PREFIX : ${GENE_ID_PREFIX} (format: {prefix}_{chrom}G${GENE_ID_DIGITS}, step=${GENE_ID_STEP})"
    echo "============================================================"

    check_file "${GENOME}" "Genome FASTA"
    check_file "${BRAKER_GFF}" "BRAKER GFF3"

############################
# 3. Standardize GFF3 by AGAT
############################

echo
echo "============================================================"
echo "Step 1. Standardize BRAKER GFF3 by AGAT"
echo "============================================================"

RAW_GFF="${BRAKER_GFF}"
GFF_STD_RAW="${OUTDIR}/01_gff/${PREFIX}.braker.standardized.gff3"
QC_GFF="${OUTDIR}/01_gff/${PREFIX}.braker.qc.gff3"
QC_MARKER="${OUTDIR}/02_sequences/${PREFIX}.qc_dot_protein.done"
QC_DROPPED="${OUTDIR}/02_sequences/${PREFIX}.dropped_dot_protein_transcripts.tsv"

if outputs_ready "${GFF_STD_RAW}"; then
    log_skip "Step 1: standardized GFF3 exists -> ${GFF_STD_RAW}"
else
if [[ -x "${ENV_AGAT}/agat_convert_sp_gxf2gxf.pl" ]]; then
    PATH="${ENV_AGAT}:$PATH" agat_convert_sp_gxf2gxf.pl \
        -g "${RAW_GFF}" \
        -o "${GFF_STD_RAW}"
else
    echo "[WARN] AGAT convert unavailable. Copy raw GFF3 directly."
    cp "${RAW_GFF}" "${GFF_STD_RAW}"
fi
fi

check_file "${GFF_STD_RAW}" "Standardized GFF3"


############################
# 4. Extract transcript/CDS/protein by gffread
############################

echo
echo "============================================================"
echo "Step 2. Extract transcript, CDS, and protein sequences"
echo "============================================================"

TRANSCRIPT_FA="${OUTDIR}/02_sequences/${PREFIX}.transcripts.fa"
CDS_FA="${OUTDIR}/02_sequences/${PREFIX}.cds.fa"
PROTEIN_FA="${OUTDIR}/02_sequences/${PREFIX}.protein.fa"

if outputs_ready "${TRANSCRIPT_FA}" "${CDS_FA}" "${PROTEIN_FA}" "${QC_GFF}" "${QC_MARKER}"; then
    log_skip "Step 2: sequences + QC GFF3 exist"
else
PATH="${ENV_ANNP}:$PATH" gffread "${GFF_STD_RAW}" \
    -g "${GENOME}" \
    -w "${TRANSCRIPT_FA}" \
    -x "${CDS_FA}" \
    -y "${PROTEIN_FA}"

echo
echo "============================================================"
echo "Step 2b. Remove transcripts with '.' in protein (gffread/DIAMOND)"
echo "============================================================"

export QC_PROTEIN_FA="${PROTEIN_FA}"
export QC_TRANSCRIPT_FA="${TRANSCRIPT_FA}"
export QC_CDS_FA="${CDS_FA}"
export QC_IN_GFF="${GFF_STD_RAW}"
export QC_OUT_GFF="${QC_GFF}"
export QC_DROPPED_LIST="${QC_DROPPED}"
export QC_MARKER="${QC_MARKER}"

python3 "${SCRIPT_DIR}/filter_dot_protein_transcripts.py"
fi

# Downstream steps use QC-filtered annotation
STD_GFF="${QC_GFF}"

check_file "${TRANSCRIPT_FA}" "Transcript fasta"
check_file "${CDS_FA}" "CDS fasta"
check_file "${PROTEIN_FA}" "Protein fasta"
check_file "${STD_GFF}" "QC GFF3"

echo "[INFO] Number of transcripts:"
grep -c "^>" "${TRANSCRIPT_FA}" || true

echo "[INFO] Number of CDS:"
grep -c "^>" "${CDS_FA}" || true

echo "[INFO] Number of proteins:"
grep -c "^>" "${PROTEIN_FA}" || true

echo "[INFO] Proteins containing stop codon '*':"
grep -v "^>" "${PROTEIN_FA}" | grep -c "\*" || true

if [[ -s "${QC_DROPPED}" ]]; then
    n_drop="$(awk 'NR>1' "${QC_DROPPED}" | wc -l)"
    echo "[INFO] Dropped transcripts (dot in protein): ${n_drop} (see ${QC_DROPPED})"
fi


############################
# 2c. Ortholog-guided CDS N-terminal start correction (ZH13v3 + Wm82 gnm6)
############################

echo
echo "============================================================"
echo "Step 2c. CDS start correction vs ZH13v3 + Wm82 gnm6 (before homolog)"
echo "============================================================"

START_FIX_DIR="${OUTDIR}/02_sequences/start_fix"
START_FIX_GFF="${START_FIX_DIR}/${PREFIX}.braker.start_fixed.gff3"
START_FIX_EVIDENCE="${START_FIX_DIR}/${PREFIX}.start_codon_correction.tsv"
START_FIX_MARKER="${START_FIX_DIR}/${PREFIX}.start_fix.done"
START_FIX_ZH13_BLAST="${START_FIX_DIR}/${PREFIX}.vs_ZH13v3_startfix.tsv"
START_FIX_WM82_BLAST="${START_FIX_DIR}/${PREFIX}.vs_Wm82_gnm6_startfix.tsv"
START_FIX_ZH13_DMND="${REF_DMND_DIR}/ZH13v3.0_mRNA_nuclear_20250729_pro.dmnd"
START_FIX_WM82_DMND="${REF_DMND_DIR}/glyma.Wm82.gnm6.ann1.PKSW.protein.dmnd"

if [[ "${RUN_START_FIX}" -eq 1 ]]; then
    mkdir -p "${START_FIX_DIR}"
    check_file "${START_FIX_ZH13_FA}" "ZH13v3 proteins for start fix"
    check_file "${START_FIX_WM82_FA}" "Wm82 gnm6 proteins for start fix"
    [[ -f "${START_FIX_SCRIPT}" ]] || die "start-fix script missing: ${START_FIX_SCRIPT}"

    if outputs_ready "${START_FIX_GFF}" "${START_FIX_EVIDENCE}" "${START_FIX_MARKER}" "${PROTEIN_FA}" "${CDS_FA}" "${TRANSCRIPT_FA}"; then
        log_skip "Step 2c: start-fix outputs exist -> ${START_FIX_EVIDENCE}"
        STD_GFF="${START_FIX_GFF}"
    else
        # Build DBs if needed
        if [[ ! -s "${START_FIX_ZH13_DMND}" ]]; then
            echo "[INFO] Building DIAMOND DB: ${START_FIX_ZH13_DMND}"
            PATH="${ENV_ANNP}:$PATH" diamond makedb --in "${START_FIX_ZH13_FA}" -d "${START_FIX_ZH13_DMND%.dmnd}"
        fi
        if [[ ! -s "${START_FIX_WM82_DMND}" ]]; then
            echo "[INFO] Building DIAMOND DB: ${START_FIX_WM82_DMND}"
            PATH="${ENV_ANNP}:$PATH" diamond makedb --in "${START_FIX_WM82_FA}" -d "${START_FIX_WM82_DMND%.dmnd}"
        fi

        echo "[INFO] DIAMOND vs ZH13v3 (start-fix screening)"
        PATH="${ENV_ANNP}:$PATH" diamond blastp \
            -q "${PROTEIN_FA}" \
            -d "${START_FIX_ZH13_DMND}" \
            -o "${START_FIX_ZH13_BLAST}" \
            -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
            --evalue 1e-5 \
            --max-target-seqs 5 \
            --threads "${THREADS}"

        echo "[INFO] DIAMOND vs Wm82 gnm6 (start-fix screening)"
        PATH="${ENV_ANNP}:$PATH" diamond blastp \
            -q "${PROTEIN_FA}" \
            -d "${START_FIX_WM82_DMND}" \
            -o "${START_FIX_WM82_BLAST}" \
            -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
            --evalue 1e-5 \
            --max-target-seqs 5 \
            --threads "${THREADS}"

        echo "[INFO] Applying start-codon corrections"
        python3 "${START_FIX_SCRIPT}" \
            --genome "${GENOME}" \
            --gff-in "${STD_GFF}" \
            --gff-out "${START_FIX_GFF}" \
            --proteins "${PROTEIN_FA}" \
            --zh13-fa "${START_FIX_ZH13_FA}" \
            --wm82-fa "${START_FIX_WM82_FA}" \
            --zh13-blast "${START_FIX_ZH13_BLAST}" \
            --wm82-blast "${START_FIX_WM82_BLAST}" \
            --evidence-tsv "${START_FIX_EVIDENCE}" \
            --min-pident "${START_FIX_MIN_PIDENT}" \
            --min-nterm-ident "${START_FIX_MIN_NTERM_IDENT}" \
            --min-missing-aa "${START_FIX_MIN_MISSING_AA}"

        check_file "${START_FIX_GFF}" "start-fixed GFF3"
        check_file "${START_FIX_EVIDENCE}" "start-fix evidence TSV"

        n_fix="$(awk -F'\t' 'NR>1 && $5=="corrected"' "${START_FIX_EVIDENCE}" | wc -l)"
        echo "[INFO] Start-fix corrected transcripts: ${n_fix}"
        echo "[INFO] Evidence table: ${START_FIX_EVIDENCE}"

        # Re-extract sequences from corrected GFF (downstream homolog uses these)
        echo "[INFO] Re-extract transcript/CDS/protein after start fix"
        PATH="${ENV_ANNP}:$PATH" gffread "${START_FIX_GFF}" \
            -g "${GENOME}" \
            -w "${TRANSCRIPT_FA}" \
            -x "${CDS_FA}" \
            -y "${PROTEIN_FA}"

        STD_GFF="${START_FIX_GFF}"
        date "+%F %T corrected=${n_fix}" > "${START_FIX_MARKER}"

        # Invalidate stale homolog diamond results if any correction applied
        if [[ "${n_fix}" -gt 0 ]]; then
            echo "[INFO] Corrections applied -> clearing ${OUTDIR}/05_diamond for re-run on fixed proteins"
            rm -f "${OUTDIR}/05_diamond/${PREFIX}".vs_*.tsv \
                  "${OUTDIR}/05_diamond/${PREFIX}".vs_*.best_hit.tsv 2>/dev/null || true
        fi
    fi

    check_file "${STD_GFF}" "GFF after start fix"
    check_file "${PROTEIN_FA}" "Protein fasta after start fix"
else
    echo "[INFO] Step 2c skipped because RUN_START_FIX=0"
fi


############################
# 5. Gene structure statistics
############################

echo
echo "============================================================"
echo "Step 3. Gene structure statistics"
echo "============================================================"

BASIC_STATS="${OUTDIR}/03_stats/${PREFIX}.basic_gene_statistics.tsv"
AGAT_STATS="${OUTDIR}/03_stats/${PREFIX}.agat.statistics.txt"

if outputs_ready "${BASIC_STATS}"; then
    log_skip "Step 3: gene statistics exist -> ${BASIC_STATS}"
else
if [[ -x "${ENV_AGAT}/agat_sp_statistics.pl" ]]; then
    PATH="${ENV_AGAT}:$PATH" agat_sp_statistics.pl \
        --gff "${STD_GFF}" \
        --out "${AGAT_STATS}"
else
    echo "[WARN] AGAT statistics unavailable. Skip AGAT statistics."
fi

python - <<PY
from collections import defaultdict

gff = "${STD_GFF}"
out = "${BASIC_STATS}"

genes = {}
mrnas = {}
exons_by_mrna = defaultdict(list)
cds_by_mrna = defaultdict(list)
mrna_parent = {}

def parse_attr(s):
    d = {}
    for item in s.strip().split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d

with open(gff) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        p = line.rstrip("\\n").split("\\t")
        if len(p) < 9:
            continue
        chrom, source, feature, start, end, score, strand, phase, attr = p
        start, end = int(start), int(end)
        a = parse_attr(attr)
        if feature == "gene":
            gid = a.get("ID", "")
            if gid:
                genes[gid] = (chrom, start, end, strand)
        elif feature in ["mRNA", "transcript"]:
            tid = a.get("ID", "")
            parent = a.get("Parent", "")
            if tid:
                mrnas[tid] = (chrom, start, end, strand)
                mrna_parent[tid] = parent
        elif feature == "exon":
            parent = a.get("Parent", "")
            for tid in parent.split(","):
                if tid:
                    exons_by_mrna[tid].append((start, end))
        elif feature == "CDS":
            parent = a.get("Parent", "")
            for tid in parent.split(","):
                if tid:
                    cds_by_mrna[tid].append((start, end))

gene_lengths = [e - s + 1 for _, s, e, _ in genes.values()]
mrna_lengths = [e - s + 1 for _, s, e, _ in mrnas.values()]
exon_counts = [len(exons_by_mrna[t]) for t in mrnas]
cds_lengths = []
for t in mrnas:
    total = sum(e - s + 1 for s, e in cds_by_mrna[t])
    if total > 0:
        cds_lengths.append(total)

def mean(x):
    return sum(x) / len(x) if x else 0

def median(x):
    if not x:
        return 0
    x = sorted(x)
    n = len(x)
    return x[n // 2] if n % 2 else (x[n // 2 - 1] + x[n // 2]) / 2

with open(out, "w") as w:
    w.write("metric\\tvalue\\n")
    w.write(f"gene_number\\t{len(genes)}\\n")
    w.write(f"mRNA_number\\t{len(mrnas)}\\n")
    w.write(f"mean_gene_length\\t{mean(gene_lengths):.2f}\\n")
    w.write(f"median_gene_length\\t{median(gene_lengths):.2f}\\n")
    w.write(f"mean_mRNA_length\\t{mean(mrna_lengths):.2f}\\n")
    w.write(f"median_mRNA_length\\t{median(mrna_lengths):.2f}\\n")
    w.write(f"mean_CDS_length\\t{mean(cds_lengths):.2f}\\n")
    w.write(f"median_CDS_length\\t{median(cds_lengths):.2f}\\n")
    w.write(f"mean_exon_number_per_mRNA\\t{mean(exon_counts):.2f}\\n")
    w.write(f"single_exon_mRNA_number\\t{sum(1 for x in exon_counts if x == 1)}\\n")

print(f"[INFO] Basic statistics written to {out}")
PY
fi


############################
# 6. BUSCO protein completeness
############################

echo
echo "============================================================"
echo "Step 4. BUSCO protein completeness"
echo "============================================================"

BUSCO_OUT_DIR="${OUTDIR}/04_busco/${PREFIX}_protein_busco"

if [[ "${RUN_BUSCO}" -eq 1 ]]; then
    if outputs_ready_glob "${BUSCO_OUT_DIR}/short_summary"*.txt; then
        log_skip "Step 4: BUSCO results exist -> ${BUSCO_OUT_DIR}"
    else
    check_file "${BUSCO_LINEAGE}/dataset.cfg" "BUSCO lineage dataset.cfg"

    (
        export PATH="${ENV_BUSCO}:$PATH"

        cd "${OUTDIR}/04_busco"

        busco \
            -m "${BUSCO_MODE}" \
            -i "${PROTEIN_FA}" \
            -l "${BUSCO_LINEAGE}" \
            --offline \
            -o "${PREFIX}_protein_busco" \
            -c "${THREADS}" \
            --force
    )
    fi
else
    echo "[INFO] BUSCO skipped because RUN_BUSCO=0"
fi


############################
# 7. DIAMOND SwissProt annotation
############################

echo
echo "============================================================"
echo "Step 5. DIAMOND SwissProt annotation"
echo "============================================================"

SWISSPROT_OUT="${OUTDIR}/05_diamond/${PREFIX}.vs_swissprot.tsv"

if [[ "${RUN_DIAMOND}" -eq 1 ]]; then
    if outputs_ready "${SWISSPROT_OUT}"; then
        log_skip "Step 5: SwissProt DIAMOND exists -> ${SWISSPROT_OUT}"
    else
    check_file "${SWISSPROT_DMND}" "SwissProt DIAMOND database"

    PATH="${ENV_ANNP}:$PATH" diamond blastp \
        -q "${PROTEIN_FA}" \
        -d "${SWISSPROT_DMND}" \
        -o "${SWISSPROT_OUT}" \
        -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
        --evalue 1e-5 \
        --max-target-seqs 5 \
        --threads "${THREADS}"
    fi
else
    echo "[INFO] DIAMOND skipped because RUN_DIAMOND=0"
fi


############################
# 7b. DIAMOND vs all reference proteomes in blastdb
############################

echo
echo "============================================================"
echo "Step 5b. DIAMOND vs reference proteomes (fixed 5: ZH13 v2/v3 + Wm82 a2.v1/gnm4/gnm6)"
echo "============================================================"

REF_PROT_MANIFEST="${OUTDIR}/05_diamond/ref_proteome_manifest.tsv"

proteome_label_from_path() {
    local base
    base="$(basename "$1")"
    base="${base%.gz}"
    base="${base%.fasta}"
    base="${base%.faa}"
    base="${base%.fa}"
    base="${base%.pep}"
    # keep readable but file-safe
    echo "${base}" | sed 's/[^A-Za-z0-9._-]/_/g'
}

discover_ref_proteomes() {
    # 兼容旧调用：直接返回写死的 REF_PROT_FASTAS
    local fa
    if [[ "${#REF_PROT_FASTAS[@]}" -eq 0 ]]; then
        echo "[ERROR] REF_PROT_FASTAS is empty" >&2
        return 1
    fi
    for fa in "${REF_PROT_FASTAS[@]}"; do
        if [[ ! -s "$fa" ]]; then
            echo "[ERROR] Missing reference proteome: $fa" >&2
            return 1
        fi
        printf '%s\n' "$fa"
    done
}

run_ref_proteome_diamond() {
    local label="$1"
    local pep_fa="$2"
    local dmnd_db="$3"
    local blast_out="$4"
    local best_out="$5"

    if outputs_ready "${blast_out}" "${best_out}"; then
        log_skip "Step 5b (${label}): DIAMOND results exist -> ${best_out}"
        return 0
    fi

    check_file "${pep_fa}" "${label} reference proteome"
    mkdir -p "${REF_DMND_DIR}"

    if [[ ! -s "${dmnd_db}" ]]; then
        echo "[INFO] Building DIAMOND DB for ${label}: ${dmnd_db}"
        PATH="${ENV_ANNP}:$PATH" diamond makedb \
            --in "${pep_fa}" \
            -d "${dmnd_db%.dmnd}"
    fi
    check_file "${dmnd_db}" "${label} DIAMOND database"

    echo "[INFO] DIAMOND blastp vs ${label}"
    PATH="${ENV_ANNP}:$PATH" diamond blastp \
        -q "${PROTEIN_FA}" \
        -d "${dmnd_db}" \
        -o "${blast_out}" \
        -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
        --evalue "${REF_DIAMOND_EVALUE}" \
        --max-target-seqs "${REF_DIAMOND_MAX_TARGET_SEQS}" \
        --threads "${THREADS}"

    python - <<PY
import os

blast = "${blast_out}"
best = "${best_out}"
label = "${label}"

best_hit = {}
if os.path.exists(blast) and os.path.getsize(blast) > 0:
    with open(blast) as f:
        for line in f:
            p = line.rstrip("\\n").split("\\t")
            if len(p) < 12:
                continue
            q = p[0]
            if q not in best_hit:
                best_hit[q] = p

with open(best, "w") as w:
    w.write("query_id\\tref_id\\tpident\\tlength\\tevalue\\tbitscore\\tref_label\\n")
    for q in sorted(best_hit):
        p = best_hit[q]
        w.write(f"{p[0]}\\t{p[1]}\\t{p[2]}\\t{p[3]}\\t{p[10]}\\t{p[11]}\\t{label}\\n")

print(f"[INFO] {label} best hits: {len(best_hit)} -> {best}")
PY
}

if [[ "${RUN_REF_PROTEOME}" -eq 1 ]]; then
    mapfile -t REF_PROT_FASTAS < <(discover_ref_proteomes)
    echo "[INFO] Using ${#REF_PROT_FASTAS[@]} fixed reference proteomes (ZH13 v2/v3 + Wm82 a2.v1/gnm4/gnm6)"
    mkdir -p "${OUTDIR}/05_diamond" "${REF_DMND_DIR}"
    printf 'label\tpep_fa\tdmnd\tblast_out\tbest_out\n' > "${REF_PROT_MANIFEST}"

    for pep_fa in "${REF_PROT_FASTAS[@]}"; do
        label="$(proteome_label_from_path "${pep_fa}")"
        dmnd_db="${REF_DMND_DIR}/${label}.dmnd"
        blast_out="${OUTDIR}/05_diamond/${PREFIX}.vs_${label}.tsv"
        best_out="${OUTDIR}/05_diamond/${PREFIX}.vs_${label}.best_hit.tsv"
        printf '%s\t%s\t%s\t%s\t%s\n' \
            "${label}" "${pep_fa}" "${dmnd_db}" "${blast_out}" "${best_out}" \
            >> "${REF_PROT_MANIFEST}"
        echo "[INFO] Ref proteome: ${label} <- ${pep_fa}"
        run_ref_proteome_diamond "${label}" "${pep_fa}" "${dmnd_db}" "${blast_out}" "${best_out}"
    done
    echo "[INFO] Reference proteome manifest -> ${REF_PROT_MANIFEST}"
else
    echo "[INFO] Reference proteome DIAMOND skipped because RUN_REF_PROTEOME=0"
    mkdir -p "${OUTDIR}/05_diamond"
    printf 'label\tpep_fa\tdmnd\tblast_out\tbest_out\n' > "${REF_PROT_MANIFEST}"
fi


############################
# 8. Convert RepeatMasker .out to BED and calculate TE overlap
############################

echo
echo "============================================================"
echo "Step 6. RepeatMasker-based TE overlap"
echo "============================================================"

GENE_BED="${OUTDIR}/08_TE_overlap/${PREFIX}.genes.bed"
CDS_BED="${OUTDIR}/08_TE_overlap/${PREFIX}.CDS.bed"
RM_BED="${OUTDIR}/08_TE_overlap/${PREFIX}.RepeatMasker.TE.bed"
GENE_TE_OVERLAP="${OUTDIR}/08_TE_overlap/${PREFIX}.gene_TE_overlap.raw.tsv"
CDS_TE_OVERLAP="${OUTDIR}/08_TE_overlap/${PREFIX}.CDS_TE_overlap.raw.tsv"
GENE_TE_SUMMARY="${OUTDIR}/08_TE_overlap/${PREFIX}.gene_TE_overlap_summary.tsv"

if [[ "${RUN_TE_OVERLAP}" -eq 1 ]]; then
    if outputs_ready "${GENE_TE_SUMMARY}" "${RM_BED}"; then
        log_skip "Step 6: TE overlap results exist -> ${GENE_TE_SUMMARY}"
    else
    check_file "${RM_OUT}" "RepeatMasker .out file"

    python - <<PY
gff = "${STD_GFF}"
gene_bed = "${GENE_BED}"
cds_bed = "${CDS_BED}"

def parse_attr(s):
    d = {}
    for item in s.strip().split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d

with open(gff) as f, open(gene_bed, "w") as wg, open(cds_bed, "w") as wc:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        p = line.rstrip("\\n").split("\\t")
        if len(p) < 9:
            continue
        chrom, source, feature, start, end, score, strand, phase, attr = p
        start0 = int(start) - 1
        end1 = int(end)
        a = parse_attr(attr)

        if feature == "gene":
            gid = a.get("ID", "NA")
            wg.write(f"{chrom}\\t{start0}\\t{end1}\\t{gid}\\n")

        elif feature == "CDS":
            parent = a.get("Parent", "NA")
            wc.write(f"{chrom}\\t{start0}\\t{end1}\\t{parent}\\n")
PY

    # Convert RepeatMasker .out to BED (TE superfamilies only).
    # Standard RepeatMasker .out columns:
    # SW perc_div perc_del perc_ins query begin end left strand repeat class/family r_begin r_end r_left ID
    python - <<PY
rm_out = "${RM_OUT}"
rm_bed = "${RM_BED}"
te_super = {x.strip() for x in "${TE_SUPERFAMILIES}".split(",") if x.strip()}

def is_te(repeat_class: str) -> bool:
    return repeat_class.split("/")[0] in te_super

n_in = 0
n_te = 0
n_skip = 0
with open(rm_out) as f, open(rm_bed, "w") as w:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith("SW") or line.startswith("score") or line.startswith("There were") or line.startswith("No repetitive"):
            continue
        parts = line.split()
        if len(parts) < 11:
            continue
        n_in += 1
        try:
            chrom = parts[4]
            start = int(parts[5]) - 1
            end = int(parts[6])
            strand = parts[8]
            repeat_name = parts[9]
            repeat_class = parts[10]
        except Exception:
            n_skip += 1
            continue

        if not is_te(repeat_class):
            n_skip += 1
            continue

        if end <= start:
            n_skip += 1
            continue

        w.write(f"{chrom}\\t{start}\\t{end}\\t{repeat_name}\\t{repeat_class}\\t{strand}\\n")
        n_te += 1

print(f"[INFO] RepeatMasker .out records parsed: {n_in}")
print(f"[INFO] TE-only BED records kept: {n_te} (superfamilies: {', '.join(sorted(te_super))})")
print(f"[INFO] Skipped (non-TE or invalid): {n_skip}")
print(f"[INFO] RepeatMasker TE BED: {rm_bed}")
PY

    check_file "${RM_BED}" "RepeatMasker BED"

    PATH="${ENV_ANNP}:$PATH" bedtools intersect \
        -a "${GENE_BED}" \
        -b "${RM_BED}" \
        -wo \
        > "${GENE_TE_OVERLAP}" || true

    PATH="${ENV_ANNP}:$PATH" bedtools intersect \
        -a "${CDS_BED}" \
        -b "${RM_BED}" \
        -wo \
        > "${CDS_TE_OVERLAP}" || true

    python - <<PY
from collections import defaultdict

gene_bed = "${GENE_BED}"
overlap = "${GENE_TE_OVERLAP}"
out = "${GENE_TE_SUMMARY}"

gene_len = {}
with open(gene_bed) as f:
    for line in f:
        p = line.rstrip("\\n").split("\\t")
        if len(p) < 4:
            continue
        chrom, start, end, gene = p[:4]
        gene_len[gene] = int(end) - int(start)

overlap_bp = defaultdict(int)
repeat_names = defaultdict(set)
repeat_classes = defaultdict(set)

try:
    with open(overlap) as f:
        for line in f:
            p = line.rstrip("\\n").split("\\t")
            if len(p) < 11:
                continue
            gene = p[3]
            rep_name = p[7]
            rep_class = p[8]
            bp = int(p[-1])
            overlap_bp[gene] += bp
            repeat_names[gene].add(rep_name)
            repeat_classes[gene].add(rep_class)
except FileNotFoundError:
    pass

with open(out, "w") as w:
    w.write("gene_id\\tgene_length\\tTE_overlap_bp\\tTE_overlap_ratio\\tTE_repeat_names\\tTE_repeat_classes\\n")
    for gene in sorted(gene_len):
        length = gene_len[gene]
        bp = overlap_bp.get(gene, 0)
        ratio = bp / length if length > 0 else 0
        names = ",".join(sorted(repeat_names[gene])) if gene in repeat_names else "NA"
        classes = ",".join(sorted(repeat_classes[gene])) if gene in repeat_classes else "NA"
        w.write(f"{gene}\\t{length}\\t{bp}\\t{ratio:.4f}\\t{names}\\t{classes}\\n")

print(f"[INFO] TE overlap summary written to {out}")
PY
    fi

else
    echo "[INFO] TE overlap skipped because RUN_TE_OVERLAP=0"
fi


############################
# 9. eggNOG-mapper annotation
############################

echo
echo "============================================================"
echo "Step 7. eggNOG-mapper annotation"
echo "============================================================"

EGGNOG_DIR="${OUTDIR}/06_eggnog"
EGGNOG_ANNOT="${EGGNOG_DIR}/${PREFIX}_eggnog.emapper.annotations"

if [[ "${RUN_EGGNOG}" -eq 1 ]]; then
    if outputs_ready "${EGGNOG_ANNOT}"; then
        log_skip "Step 7: eggNOG annotations exist -> ${EGGNOG_ANNOT}"
    else
    check_file "${EMAPPER_PY}" "emapper.py"
    check_dir "${EGGNOG_DATA_DIR}" "eggNOG database directory"

    echo "[INFO] eggNOG data_dir: ${EGGNOG_DATA_DIR}"
    # eggnog213 没有 diamond；系统 PATH 上的 0.9.19 不认 --iterate。
    # python 仍走 eggnog213，diamond 走 annp 2.x（与 Step 5 相同）。
    EGGNOG_PATH="${ENV_EGGNOG}:${ENV_ANNP}:$PATH"
    echo "[INFO] eggNOG diamond: $(PATH="${EGGNOG_PATH}" command -v diamond)"

    PATH="${EGGNOG_PATH}" "${EMAPPER_PY}" \
        -i "${PROTEIN_FA}" \
        --itype proteins \
        --data_dir "${EGGNOG_DATA_DIR}" \
        -o "${PREFIX}_eggnog" \
        --output_dir "${EGGNOG_DIR}" \
        --cpu "${THREADS}" \
        --override
    fi
else
    echo "[INFO] eggNOG skipped because RUN_EGGNOG=0"
fi


############################
# 10. InterProScan annotation, optional
############################

echo
echo "============================================================"
echo "Step 8. InterProScan annotation"
echo "============================================================"

# -b is output file base (InterProScan appends .tsv / .gff3; do NOT use -o with multiple -f)
INTERPRO_BASE="${OUTDIR}/07_interpro/${PREFIX}.interproscan"
INTERPRO_TSV="${INTERPRO_BASE}.tsv"
INTERPRO_GFF="${INTERPRO_BASE}.gff3"

if [[ "${RUN_INTERPRO}" -eq 1 ]]; then
    if outputs_ready "${INTERPRO_TSV}"; then
        log_skip "Step 8: InterProScan TSV exists -> ${INTERPRO_TSV}"
    else
    if [[ ! -x "${INTERPROSCAN_BIN}" ]]; then
        echo "[ERROR] InterProScan not executable: ${INTERPROSCAN_BIN}" >&2
        exit 1
    fi

    mkdir -p "${OUTDIR}/07_interpro"
    echo "[INFO] InterProScan: ${INTERPROSCAN_BIN}"
    echo "[INFO] Input proteins: ${PROTEIN_FA}"
    echo "[INFO] Output base: ${INTERPRO_BASE}"
    echo "[INFO] InterProScan may take several hours for a large proteome."

    "${INTERPROSCAN_BIN}" \
        -i "${PROTEIN_FA}" \
        -f tsv,gff3 \
        -dp \
        -goterms \
        -pa \
        -cpu "${THREADS}" \
        -b "${INTERPRO_BASE}"

    check_file "${INTERPRO_TSV}" "InterProScan TSV (${INTERPRO_TSV})"
    if [[ -s "${INTERPRO_GFF}" ]]; then
        echo "[INFO] InterProScan GFF3: ${INTERPRO_GFF}"
    else
        echo "[WARN] InterProScan GFF3 not found: ${INTERPRO_GFF}"
    fi
    fi
else
    echo "[INFO] InterProScan skipped because RUN_INTERPRO=0"
fi


############################
# 11. Salmon expression support, optional
############################

echo
echo "============================================================"
echo "Step 9. Salmon RNA-seq expression support"
echo "============================================================"

SALMON_INDEX="${OUTDIR}/09_salmon/${PREFIX}.salmon_index"
SALMON_QUANT_DIR="${OUTDIR}/09_salmon/quants"
EXPR_SUPPORT="${OUTDIR}/09_salmon/${PREFIX}.gene_expression_support.tsv"

if [[ "${RUN_SALMON}" -eq 1 ]]; then
    if outputs_ready "${EXPR_SUPPORT}" "${SALMON_QUANT_DIR}/${RNA_SAMPLE}/quant.sf"; then
        log_skip "Step 9: Salmon results exist -> ${EXPR_SUPPORT}"
    else
    check_file "${RNA_R1}" "RNA-seq R1"
    check_file "${RNA_R2}" "RNA-seq R2"

    mkdir -p "${SALMON_QUANT_DIR}"

    echo "[INFO] Salmon RNA-seq: ${RNA_SAMPLE}"
    echo "[INFO]   R1: ${RNA_R1}"
    echo "[INFO]   R2: ${RNA_R2}"

    PATH="${ENV_ANNP}:$PATH" salmon index \
        -t "${TRANSCRIPT_FA}" \
        -i "${SALMON_INDEX}" \
        -k 31 \
        -p "${THREADS}"

    PATH="${ENV_ANNP}:$PATH" salmon quant \
        -i "${SALMON_INDEX}" \
        -l "${SALMON_LIBTYPE}" \
        -1 "${RNA_R1}" \
        -2 "${RNA_R2}" \
        -p "${THREADS}" \
        --validateMappings \
        -o "${SALMON_QUANT_DIR}/${RNA_SAMPLE}"

    echo "[INFO] Salmon quantification finished for ${RNA_SAMPLE}."

    python - <<PY
import os
import glob
import csv
from collections import defaultdict

quant_dir = "${SALMON_QUANT_DIR}"
out_gene = "${EXPR_SUPPORT}"
tpm_cutoff = float("${TPM_CUTOFF}")
gff = "${STD_GFF}"

def parse_attr(s):
    d = {}
    for item in s.strip().split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d

tx2gene = {}
genes = set()

with open(gff) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        p = line.rstrip("\\n").split("\\t")
        if len(p) < 9:
            continue
        feature = p[2]
        a = parse_attr(p[8])
        if feature == "gene":
            gid = a.get("ID")
            if gid:
                genes.add(gid)
        elif feature in ["mRNA", "transcript"]:
            tid = a.get("ID")
            parent = a.get("Parent")
            if tid and parent:
                tx2gene[tid] = parent.split(",")[0]

gene_max_tpm = defaultdict(float)
gene_expressed_samples = defaultdict(int)
samples = []

for q in sorted(glob.glob(os.path.join(quant_dir, "*", "quant.sf"))):
    sample = os.path.basename(os.path.dirname(q))
    samples.append(sample)
    sample_gene_tpm = defaultdict(float)

    with open(q) as f:
        r = csv.DictReader(f, delimiter="\\t")
        for row in r:
            tid = row["Name"]
            tpm = float(row["TPM"])
            gid = tx2gene.get(tid)
            if gid:
                sample_gene_tpm[gid] += tpm

    for gid, tpm in sample_gene_tpm.items():
        if tpm > gene_max_tpm[gid]:
            gene_max_tpm[gid] = tpm
        if tpm >= tpm_cutoff:
            gene_expressed_samples[gid] += 1

with open(out_gene, "w") as w:
    w.write("gene_id\\tmax_TPM\\texpressed_sample_count\\texpression_supported\\n")
    for gid in sorted(genes):
        max_tpm = gene_max_tpm.get(gid, 0.0)
        count = gene_expressed_samples.get(gid, 0)
        supported = "yes" if count > 0 else "no"
        w.write(f"{gid}\\t{max_tpm:.4f}\\t{count}\\t{supported}\\n")

print(f"[INFO] Gene expression support written to {out_gene}")
print(f"[INFO] Salmon samples parsed: {len(samples)}")
PY
    fi

else
    echo "[INFO] Salmon skipped because RUN_SALMON=0"
fi


############################
# 12. Final summary table
############################

echo
echo "============================================================"
echo "Step 10. Generate final annotation summary table"
echo "============================================================"

FINAL_SUMMARY="${OUTDIR}/10_summary/${PREFIX}.annotation.summary.tsv"

if outputs_ready "${FINAL_SUMMARY}"; then
    log_skip "Step 10: annotation summary exists -> ${FINAL_SUMMARY}"
else
python - <<PY
import os
import csv
from collections import defaultdict

gff = "${STD_GFF}"
out = "${FINAL_SUMMARY}"

swiss = "${SWISSPROT_OUT}"
ref_manifest = "${REF_PROT_MANIFEST}"
eggnog = "${EGGNOG_ANNOT}"
interpro = "${INTERPRO_TSV}"
te_summary = "${GENE_TE_SUMMARY}"
expr = "${EXPR_SUPPORT}"

def parse_attr(s):
    d = {}
    for item in s.strip().split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d

genes = {}
tx2gene = {}
gene_tx = defaultdict(list)
exon_count = defaultdict(int)
cds_len = defaultdict(int)

with open(gff) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        p = line.rstrip("\\n").split("\\t")
        if len(p) < 9:
            continue
        chrom, source, feature, start, end, score, strand, phase, attr = p
        a = parse_attr(attr)

        if feature == "gene":
            gid = a.get("ID")
            if gid:
                genes[gid] = {
                    "chrom": chrom,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                }

        elif feature in ["mRNA", "transcript"]:
            tid = a.get("ID")
            parent = a.get("Parent")
            if tid and parent:
                gid = parent.split(",")[0]
                tx2gene[tid] = gid
                gene_tx[gid].append(tid)

        elif feature == "exon":
            parent = a.get("Parent", "")
            for tid in parent.split(","):
                gid = tx2gene.get(tid)
                if gid:
                    exon_count[gid] += 1

        elif feature == "CDS":
            parent = a.get("Parent", "")
            length = int(end) - int(start) + 1
            for tid in parent.split(","):
                gid = tx2gene.get(tid)
                if gid:
                    cds_len[gid] += length

def load_ref_best(path, tx2gene_map):
    gene_hit = {}
    if not (os.path.exists(path) and os.path.getsize(path) > 0):
        return gene_hit
    with open(path) as f:
        f.readline()
        for line in f:
            p = line.rstrip("\\n").split("\\t")
            if len(p) < 6:
                continue
            q, ref, pident, length, evalue, bitscore = p[:6]
            gid = tx2gene_map.get(q, q.split(".t")[0] if ".t" in q else q)
            if gid not in gene_hit:
                gene_hit[gid] = {
                    "hit": ref,
                    "pident": pident,
                    "length": length,
                    "evalue": evalue,
                    "bitscore": bitscore,
                }
    return gene_hit

# Reference proteome best hits (all proteomes in blastdb manifest)
ref_hits = {}  # label -> {gid -> hit dict}
ref_labels = []
if os.path.exists(ref_manifest) and os.path.getsize(ref_manifest) > 0:
    with open(ref_manifest) as mf:
        mf.readline()
        for line in mf:
            p = line.rstrip("\\n").split("\\t")
            if len(p) < 5:
                continue
            label, _pep, _dmnd, _blast, best_path = p[:5]
            ref_labels.append(label)
            ref_hits[label] = load_ref_best(best_path, tx2gene)

# SwissProt top hit
swiss_hit = {}
if os.path.exists(swiss) and os.path.getsize(swiss) > 0:
    with open(swiss) as f:
        for line in f:
            p = line.rstrip("\\n").split("\\t")
            if len(p) >= 13:
                q = p[0]
                gid = tx2gene.get(q, q.split(".t")[0] if ".t" in q else q)
                if gid not in swiss_hit:
                    swiss_hit[gid] = {
                        "hit": p[1],
                        "pident": p[2],
                        "evalue": p[10],
                        "bitscore": p[11],
                        "desc": p[12],
                    }

# eggNOG
eggnog_anno = {}
if os.path.exists(eggnog) and os.path.getsize(eggnog) > 0:
    with open(eggnog) as f:
        header = None
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#query"):
                header = line.lstrip("#").rstrip("\\n").split("\\t")
                continue
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\\n").split("\\t")
            if header and len(p) == len(header):
                row = dict(zip(header, p))
                q = row.get("query", p[0])
                gid = tx2gene.get(q, q.split(".t")[0] if ".t" in q else q)
                eggnog_anno[gid] = {
                    "preferred_name": row.get("Preferred_name", ""),
                    "go": row.get("GOs", ""),
                    "ko": row.get("KEGG_ko", ""),
                    "pathway": row.get("KEGG_Pathway", ""),
                    "pfam": row.get("PFAMs", ""),
                    "desc": row.get("Description", ""),
                }

# InterPro
interpro_anno = defaultdict(lambda: {"ipr": set(), "pfam": set(), "go": set()})
if os.path.exists(interpro) and os.path.getsize(interpro) > 0:
    with open(interpro) as f:
        for line in f:
            p = line.rstrip("\\n").split("\\t")
            if len(p) < 12:
                continue
            q = p[0]
            gid = tx2gene.get(q, q.split(".t")[0] if ".t" in q else q)
            analysis = p[3]
            signature = p[4]
            ipr = p[11] if len(p) > 11 else ""
            gos = p[13] if len(p) > 13 else ""

            if ipr and ipr != "-":
                interpro_anno[gid]["ipr"].add(ipr)
            if analysis == "Pfam" and signature and signature != "-":
                interpro_anno[gid]["pfam"].add(signature)
            if gos and gos != "-":
                for go in gos.split("|"):
                    if go:
                        interpro_anno[gid]["go"].add(go)

# TE overlap
te_anno = {}
if os.path.exists(te_summary) and os.path.getsize(te_summary) > 0:
    with open(te_summary) as f:
        r = csv.DictReader(f, delimiter="\\t")
        for row in r:
            te_anno[row["gene_id"]] = row

# Expression
expr_anno = {}
if os.path.exists(expr) and os.path.getsize(expr) > 0:
    with open(expr) as f:
        r = csv.DictReader(f, delimiter="\\t")
        for row in r:
            expr_anno[row["gene_id"]] = row

with open(out, "w") as w:
    header = [
        "gene_id", "chrom", "start", "end", "strand",
        "transcript_count", "total_exon_count", "total_CDS_length",
        "SwissProt_hit", "SwissProt_identity", "SwissProt_evalue", "SwissProt_description",
    ]
    for label in ref_labels:
        header.extend([
            f"{label}_hit",
            f"{label}_identity",
            f"{label}_evalue",
            f"{label}_bitscore",
        ])
    header.extend([
        "eggNOG_preferred_name", "eggNOG_description", "eggNOG_GO", "eggNOG_KEGG_KO", "eggNOG_KEGG_pathway", "eggNOG_PFAM",
        "InterPro_IDs", "InterPro_Pfam", "InterPro_GO",
        "TE_overlap_bp", "TE_overlap_ratio", "TE_repeat_classes",
        "max_TPM", "expressed_sample_count", "expression_supported",
        "confidence_level"
    ])
    w.write("\\t".join(header) + "\\n")

    for gid in sorted(genes, key=lambda x: (genes[x]["chrom"], genes[x]["start"])):
        g = genes[gid]
        sw = swiss_hit.get(gid, {})
        eg = eggnog_anno.get(gid, {})
        ip = interpro_anno.get(gid, {"ipr": set(), "pfam": set(), "go": set()})
        te = te_anno.get(gid, {})
        ex = expr_anno.get(gid, {})

        te_ratio = float(te.get("TE_overlap_ratio", 0) or 0)
        expr_support = ex.get("expression_supported", "NA")

        any_ref_hit = any(
            ref_hits.get(label, {}).get(gid, {}).get("hit") not in [None, "", "NA"]
            for label in ref_labels
        )

        has_function = bool(
            sw.get("hit") not in [None, "", "NA"] or
            any_ref_hit or
            eg.get("desc") not in [None, "", "NA"] or
            ip["ipr"] or
            ip["pfam"] or
            ip["go"]
        )
        has_expr = expr_support == "yes"
        high_te = te_ratio >= 0.5

        if high_te and not has_expr and not has_function:
            conf = "low_confidence_possible_TE"
        elif has_expr and has_function:
            conf = "high_confidence"
        elif has_expr or has_function:
            conf = "medium_confidence"
        else:
            conf = "low_confidence"

        row = [
            gid, g["chrom"], g["start"], g["end"], g["strand"],
            len(gene_tx.get(gid, [])),
            exon_count.get(gid, 0),
            cds_len.get(gid, 0),

            sw.get("hit", "NA"),
            sw.get("pident", "NA"),
            sw.get("evalue", "NA"),
            sw.get("desc", "NA"),
        ]
        for label in ref_labels:
            rh = ref_hits.get(label, {}).get(gid, {})
            row.extend([
                rh.get("hit", "NA"),
                rh.get("pident", "NA"),
                rh.get("evalue", "NA"),
                rh.get("bitscore", "NA"),
            ])
        row.extend([
            eg.get("preferred_name", "NA"),
            eg.get("desc", "NA"),
            eg.get("go", "NA"),
            eg.get("ko", "NA"),
            eg.get("pathway", "NA"),
            eg.get("pfam", "NA"),

            ",".join(sorted(ip["ipr"])) if ip["ipr"] else "NA",
            ",".join(sorted(ip["pfam"])) if ip["pfam"] else "NA",
            ",".join(sorted(ip["go"])) if ip["go"] else "NA",

            te.get("TE_overlap_bp", "NA"),
            te.get("TE_overlap_ratio", "NA"),
            te.get("TE_repeat_classes", "NA"),

            ex.get("max_TPM", "NA"),
            ex.get("expressed_sample_count", "NA"),
            ex.get("expression_supported", "NA"),

            conf
        ])

        w.write("\\t".join(str(x).replace("\\t", " ").replace("\\n", " ") for x in row) + "\\n")

print(f"[INFO] Final annotation summary written to {out}")
PY
fi


############################
# 12b. Publication-ready annotation deliverables
############################

echo
echo "============================================================"
echo "Step 11. Publication-ready functional annotation outputs"
echo "============================================================"

check_file "${FINAL_SUMMARY}" "annotation summary table"

if outputs_ready "${GENE_FUNC_TABLE}" "${ANNOTATED_GFF}" "${HC_ANNOTATED_GFF}" "${HC_PROTEIN_FA}"; then
    log_skip "Step 11: publication outputs exist under ${PUB_DIR}"
else
export PUB_STD_GFF="${STD_GFF}"
export PUB_SUMMARY="${FINAL_SUMMARY}"
export PUB_PROTEIN_FA="${PROTEIN_FA}"
export PUB_OUT_DIR="${PUB_DIR}"
export PUB_PREFIX="${PREFIX}"
export PUB_MIN_CONFIDENCE
export PUB_MAX_TE_RATIO
export PUB_REQUIRE_EXPRESSION
export PUB_NUCLEAR_ONLY
export PUB_MIN_ZH13_PIDENT

python3 "${SCRIPT_DIR}/build_publication_annotation.py"

check_file "${GENE_FUNC_TABLE}" "gene functional annotation table"
check_file "${ANNOTATED_GFF}" "annotated GFF3"
check_file "${HC_ANNOTATED_GFF}" "high-confidence annotated GFF3"
check_file "${HC_PROTEIN_FA}" "high-confidence protein FASTA"
fi


############################
# 12c. Rename gene loci (${PREFIX}_01G000001)
############################

if [[ "${RUN_GENE_RENAME}" -eq 1 ]]; then
    echo
    echo "============================================================"
    echo "Step 12. Rename gene loci (${GENE_ID_PREFIX}_{{chrom}}G${GENE_ID_DIGITS} digits)"
    echo "============================================================"

    check_file "${ANNOTATED_GFF}" "annotated GFF3 (run Step 11 first)"

    if outputs_ready "${RENAMED_GFF}" "${GENE_ID_MAPPING}" "${RENAMED_GENE_TABLE}"; then
        log_skip "Step 12: renamed outputs exist under ${RENAME_DIR}"
    else
        mkdir -p "${RENAME_DIR}"

        export RENAME_IN_GFF="${ANNOTATED_GFF}"
        export RENAME_OUT_GFF="${RENAMED_GFF}"
        export RENAME_PREFIX="${PREFIX}"
        export RENAME_ID_PREFIX="${GENE_ID_PREFIX}"
        export RENAME_GENE_DIGITS="${GENE_ID_DIGITS}"
        export RENAME_GENE_ID_STEP="${GENE_ID_STEP}"
        export RENAME_MAPPING_TSV="${GENE_ID_MAPPING}"

        export RENAME_IN_GENE_TABLE="${GENE_FUNC_TABLE}"
        export RENAME_OUT_GENE_TABLE="${RENAMED_GENE_TABLE}"
        export RENAME_IN_SUMMARY="${FINAL_SUMMARY}"
        export RENAME_OUT_SUMMARY="${RENAMED_SUMMARY}"

        export RENAME_IN_TRANSCRIPT_FA="${TRANSCRIPT_FA}"
        export RENAME_OUT_TRANSCRIPT_FA="${RENAMED_TRANSCRIPT_FA}"
        export RENAME_IN_CDS_FA="${CDS_FA}"
        export RENAME_OUT_CDS_FA="${RENAMED_CDS_FA}"
        export RENAME_IN_PROTEIN_FA="${PROTEIN_FA}"
        export RENAME_OUT_PROTEIN_FA="${RENAMED_PROTEIN_FA}"

        export RENAME_IN_HC_GFF="${HC_ANNOTATED_GFF}"
        export RENAME_OUT_HC_GFF="${RENAMED_HC_GFF}"
        export RENAME_IN_HC_PROTEIN="${HC_PROTEIN_FA}"
        export RENAME_OUT_HC_PROTEIN="${RENAMED_HC_PROTEIN_FA}"
        export RENAME_IN_HC_LIST="${HC_GENE_LIST}"
        export RENAME_OUT_HC_LIST="${RENAMED_HC_GENE_LIST}"

        python3 "${SCRIPT_DIR}/rename_gene_loci.py"

        check_file "${RENAMED_GFF}" "renamed annotated GFF3"
        check_file "${GENE_ID_MAPPING}" "gene ID mapping table"
        check_file "${RENAMED_GENE_TABLE}" "renamed gene functional table"
        check_file "${RENAMED_PROTEIN_FA}" "renamed protein FASTA"

        echo
        echo "============================================================"
        echo "Step 12b. ZH13-style structural GFF (gene/mRNA/exon/CDS/UTR only)"
        echo "============================================================"
        python3 "${SCRIPT_DIR}/pack_structural_gff.py" \
            -i "${RENAMED_GFF}" \
            -o "${RENAMED_GFF}.structural.tmp" \
            --nuclear-out "${RENAMED_NUCLEAR_GFF}" \
            --min-aa "${PUB_MIN_AA}" \
            --genome "${GENOME}"
        mv "${RENAMED_GFF}.structural.tmp" "${RENAMED_GFF}"
        if [[ -s "${RENAMED_NUCLEAR_GFF}" ]]; then
            PATH="${ENV_ANNP}:$PATH" gffread "${RENAMED_NUCLEAR_GFF}" \
                -g "${GENOME}" -y "${RENAMED_NUCLEAR_PROTEIN_FA}" || true
        fi
        if [[ -s "${RENAMED_HC_GFF}" ]]; then
            python3 "${SCRIPT_DIR}/pack_structural_gff.py" \
                -i "${RENAMED_HC_GFF}" \
                -o "${RENAMED_HC_GFF}.structural.tmp"
            mv "${RENAMED_HC_GFF}.structural.tmp" "${RENAMED_HC_GFF}"
        fi
    fi
else
    echo "[INFO] Gene locus renaming skipped (RUN_GENE_RENAME=0)"
fi

if [[ "${RUN_GENE_RENAME}" -eq 1 && -s "${RENAMED_SUMMARY}" && -s "${RENAMED_GFF}" ]]; then
    echo
    echo "============================================================"
    echo "Step 12c. Lineage-specific set (nuclear + medium+ + TE<=0.5 + expressed + ZH13 pident < ${PUB_MIN_ZH13_PIDENT})"
    echo "============================================================"
    python3 "${SCRIPT_DIR}/extract_lineage_specific.py" \
        --dir "${RENAME_DIR}" \
        --prefix "${PREFIX}" \
        --max-zh13 "${PUB_MIN_ZH13_PIDENT}"
fi

if [[ "${RUN_GENE_RENAME}" -eq 1 && -s "${RENAMED_HC_GENE_LIST}" && -s "${RENAMED_TRANSCRIPT_FA}" ]]; then
    echo
    echo "============================================================"
    echo "Step 12d. High-confidence transcripts (mRNA) and CDS FASTA"
    echo "============================================================"
    python3 "${SCRIPT_DIR}/extract_hc_sequences.py" \
        --dir "${RENAME_DIR}" \
        --prefix "${PREFIX}"
fi


############################
# 13. Final message
############################

echo
echo "============================================================"
echo "Pipeline finished"
echo "End time: $(date)"
echo "============================================================"

echo
echo "Main outputs:"
echo
echo "1. Standardized GFF3:"
echo "   ${GFF_STD_RAW}"
echo "   ${STD_GFF} (QC, used downstream)"
echo
echo "2. Transcript/CDS/protein:"
echo "   ${TRANSCRIPT_FA}"
echo "   ${CDS_FA}"
echo "   ${PROTEIN_FA}"
if [[ "${RUN_START_FIX}" -eq 1 ]]; then
echo "2c. Start-codon correction evidence:"
echo "   ${OUTDIR}/02_sequences/start_fix/${PREFIX}.start_codon_correction.tsv"
echo "   ${OUTDIR}/02_sequences/start_fix/${PREFIX}.braker.start_fixed.gff3"
fi
echo
echo "3. Basic gene statistics:"
echo "   ${OUTDIR}/03_stats/${PREFIX}.basic_gene_statistics.tsv"
echo
echo "4. AGAT statistics:"
echo "   ${OUTDIR}/03_stats/${PREFIX}.agat.statistics.txt"
echo
echo "5. BUSCO result:"
echo "   ${OUTDIR}/04_busco/${PREFIX}_protein_busco"
echo
echo "6. SwissProt Diamond:"
echo "   ${SWISSPROT_OUT}"
echo
echo "6b. Reference proteome Diamond (blastdb/prot):"
echo "   manifest: ${REF_PROT_MANIFEST}"
if [[ -s "${REF_PROT_MANIFEST}" ]]; then
    awk -F'\t' 'NR>1 {print "   "$5}' "${REF_PROT_MANIFEST}"
fi
echo
echo "7. RepeatMasker TE BED:"
echo "   ${RM_BED}"
echo
echo "8. TE overlap summary:"
echo "   ${GENE_TE_SUMMARY}"
echo
echo "9. eggNOG annotation:"
echo "   ${EGGNOG_ANNOT}"
echo
echo "10. InterProScan:"
echo "   ${INTERPRO_TSV}"
echo "   ${INTERPRO_GFF}"
echo
echo "11. Salmon expression support:"
echo "   ${EXPR_SUPPORT}"
echo
echo "12. Final annotation summary:"
echo "   ${FINAL_SUMMARY}"
echo
echo "13. Publication-ready outputs:"
echo "   ${GENE_FUNC_TABLE}"
echo "   ${ANNOTATED_GFF}"
echo "   ${HC_GENE_LIST}"
echo "   ${HC_ANNOTATED_GFF}"
echo "   ${HC_PROTEIN_FA}"
echo "   ${PUB_DIR}/${PREFIX}.annotation_README.txt"
echo
if [[ "${RUN_GENE_RENAME}" -eq 1 ]]; then
echo "14. Renamed locus outputs (${GENE_ID_PREFIX}_XXG${GENE_ID_DIGITS}):"
echo "   ${RENAMED_GFF}"
echo "   ${GENE_ID_MAPPING}"
echo "   ${RENAMED_GENE_TABLE}"
echo "   ${RENAMED_SUMMARY}"
echo "   ${RENAMED_TRANSCRIPT_FA}"
echo "   ${RENAMED_CDS_FA}"
echo "   ${RENAMED_PROTEIN_FA}"
echo "   ${RENAMED_HC_GFF}"
echo "   ${RENAMED_HC_PROTEIN_FA}"
echo "   ${RENAMED_HC_TRANSCRIPT_FA}"
echo "   ${RENAMED_HC_CDS_FA}"
echo "   ${RENAMED_HC_GENE_LIST}"
echo "   ${RENAMED_NUCLEAR_GFF}  (Chr01–20, protein >= ${PUB_MIN_AA} aa)"
echo "   ${RENAMED_NUCLEAR_PROTEIN_FA}"
echo "   ${RENAME_DIR}/${PREFIX}.renamed.lineage_specific.genes.list"
echo "   ${RENAME_DIR}/${PREFIX}.renamed.lineage_specific.annotated.gff3"
echo "   ${RENAME_DIR}/${PREFIX}.renamed.lineage_specific.transcripts.fa"
echo "   ${RENAME_DIR}/${PREFIX}.renamed.lineage_specific.cds.fa"
echo "   ${RENAME_DIR}/${PREFIX}.renamed.lineage_specific.proteins.fa"
echo
fi
echo "Log:"
echo "   ${LOG}"

    } 2>&1 | tee -a "${LOG}"

    log "样本 $SAMPLE 完成"
done

log "========== 全部结束 =========="