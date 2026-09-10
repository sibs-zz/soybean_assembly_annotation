#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# HISAT2：软屏蔽基因组 → Atlas BAM（供 BRAKER3 使用）
#
# 基因组（01_EDTA 产物）：
#   ${sp}_TE/${sp}_Chr.fa.masked
#
# RNA 输入（手动指定，合并后放入 rna_clean/）：
#   ${sp}_RNA_Atlas_1.fq.gz / ${sp}_RNA_Atlas_2.fq.gz
#
# 输出：
#   ${sp}_RNA/${sp}_idx.*          HISAT2 索引
#   ${sp}_RNA/${sp}_RNA_sorted.bam Atlas BAM（含 .bai）
#   HISAT2_align_summary.tsv       比对汇总
#
# 用法：
#   bash 02_run_HISAT2_align.sh                # 建索引 + 比对
#   STAGE=index bash 02_run_HISAT2_align.sh    # 只建索引
#   STAGE=align bash 02_run_HISAT2_align.sh    # 只比对（已有索引）
#   STAGE=summary bash 02_run_HISAT2_align.sh  # 只汇总
#   FORCE_RERUN=1 bash 02_run_HISAT2_align.sh  # 强制重跑
# ============================================================

# ======================== 配置区 =============================

SAMPLES=(
    "KD35_primary"
    "N2L_primary"
)

ROOT="/tiandata2/zzh/sjy_assembly"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RNA_DIR="${ROOT}/rna_clean"

THREADS="${THREADS:-64}"
SORT_THREADS="${SORT_THREADS:-64}"
FORCE_RERUN="${FORCE_RERUN:-0}"
# index | align | all | summary
STAGE="${STAGE:-all}"

CONDA_ALIGNER="${CONDA_ALIGNER:-/home/zzh/anaconda3/envs/aligner}"
HISAT2="${HISAT2:-${CONDA_ALIGNER}/bin/hisat2}"
HISAT2_BUILD="${HISAT2_BUILD:-${CONDA_ALIGNER}/bin/hisat2-build}"
SAMTOOLS="${SAMTOOLS:-${CONDA_ALIGNER}/bin/samtools}"
SUMMARY_TSV="${SCRIPT_DIR}/HISAT2_align_summary.tsv"

# ============================================================

timestamp() { date "+%F %T"; }
log() { echo "[$(timestamp)] $*"; }
die() { echo "[$(timestamp)] ERROR: $*" >&2; exit 1; }

run_index() {
    local sample="$1"
    local sp="${sample%%_*}"
    local outdir masked idx_prefix idx_done

    outdir="${SCRIPT_DIR}/${sp}_RNA"
    masked="${SCRIPT_DIR}/${sp}_TE/${sp}_Chr.fa.masked"
    idx_prefix="${outdir}/${sp}_idx"
    idx_done="${idx_prefix}.1.ht2"

    [[ -s "$masked" ]] || die "缺少软屏蔽基因组：$masked（请先完成 01_EDTA.sh RepeatMasker）"
    mkdir -p "$outdir"

    if [[ -s "$idx_done" && "$FORCE_RERUN" != "1" ]]; then
        log "已有 HISAT2 索引，跳过：$idx_prefix"
        return 0
    fi
    [[ "$FORCE_RERUN" == "1" ]] && rm -f "${idx_prefix}".*.ht2 "${outdir}/hisat2_build.log"

    log "HISAT2 建索引：$sample"
    log "  genome : $masked"
    log "  index  : $idx_prefix"

    "$HISAT2_BUILD" -p "$THREADS" "$masked" "$idx_prefix" \
        >"${outdir}/hisat2_build.log" 2>&1

    [[ -s "$idx_done" ]] || die "HISAT2 索引失败：见 ${outdir}/hisat2_build.log"
    log "索引完成：$idx_prefix"
}

run_align() {
    local sample="$1"
    local sp="${sample%%_*}"
    local outdir idx_prefix r1 r2 bam logf

    outdir="${SCRIPT_DIR}/${sp}_RNA"
    idx_prefix="${outdir}/${sp}_idx"
    r1="${RNA_DIR}/${sp}_RNA_Atlas_1.fq.gz"
    r2="${RNA_DIR}/${sp}_RNA_Atlas_2.fq.gz"
    bam="${outdir}/${sp}_RNA_sorted.bam"
    logf="${outdir}/${sp}_RNA_Atlas.hisat2.log"

    [[ -s "$r1" ]] || die "缺少 Atlas R1：$r1"
    [[ -s "$r2" ]] || die "缺少 Atlas R2：$r2"
    [[ -s "${idx_prefix}.1.ht2" ]] || die "缺少索引：${idx_prefix}（先跑 STAGE=index）"

    if [[ -s "$bam" && -s "${bam}.bai" && "$FORCE_RERUN" != "1" ]]; then
        log "已有 BAM，跳过：$bam"
        return 0
    fi
    [[ "$FORCE_RERUN" == "1" ]] && rm -f "$bam" "${bam}.bai" "$logf"

    log "HISAT2 比对：$sample  Atlas"
    log "  R1     : $r1"
    log "  R2     : $r2"
    log "  output : $bam"
    log "  log    : $logf"

    "$HISAT2" -p "$THREADS" -x "$idx_prefix" -1 "$r1" -2 "$r2" 2>"$logf" \
        | "$SAMTOOLS" sort -@ "$SORT_THREADS" -o "$bam" -
    "$SAMTOOLS" index "$bam"
    log "比对完成：$bam"
}

write_summary() {
    local tmp
    tmp="$(mktemp)"
    {
        printf 'sample\tsp\tstatus\tbam\tmapped_reads\tgenome\n'
        for SAMPLE in "${SAMPLES[@]}"; do
            local sp="${SAMPLE%%_*}"
            local bam="${SCRIPT_DIR}/${sp}_RNA/${sp}_RNA_sorted.bam"
            local masked="${SCRIPT_DIR}/${sp}_TE/${sp}_Chr.fa.masked"
            local status mapped
            if [[ ! -s "$masked" ]]; then
                status="skipped_no_genome"
                printf '%s\t%s\t%s\t\t\t%s\n' "$SAMPLE" "$sp" "$status" "$masked"
            elif [[ ! -s "$bam" ]]; then
                status="missing_bam"
                printf '%s\t%s\t%s\t%s\t\t%s\n' "$SAMPLE" "$sp" "$status" "$bam" "$masked"
            else
                status="ok"
                mapped="$("$SAMTOOLS" flagstat "$bam" 2>/dev/null \
                    | grep ' mapped ' | grep -v primary | head -1 | awk '{print $1}' || true)"
                printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
                    "$SAMPLE" "$sp" "$status" "$bam" "${mapped:-}" "$masked"
            fi
        done
    } >"$tmp"
    mv "$tmp" "$SUMMARY_TSV"
    log "HISAT2 summary 已写入：$SUMMARY_TSV"
}

print_summary() {
    [[ -s "$SUMMARY_TSV" ]] || die "缺少 summary 文件：$SUMMARY_TSV"
    echo "------------------------------------------------"
    log "HISAT2 align summary"
    echo "------------------------------------------------"
    column -t -s $'\t' "$SUMMARY_TSV" 2>/dev/null || cat "$SUMMARY_TSV"
    echo "------------------------------------------------"
    while IFS=$'\t' read -r sample sp status bam mapped genome; do
        [[ "$sample" == "sample" ]] && continue
        if [[ "$status" == "ok" ]]; then
            printf '  %s: mapped=%s  %s\n' "$sample" "$mapped" "$bam"
        else
            printf '  %s: %s\n' "$sample" "$status"
        fi
    done < "$SUMMARY_TSV"
}

[[ -x "$HISAT2" || -f "$HISAT2" ]] || die "找不到 hisat2：$HISAT2"
[[ -x "$HISAT2_BUILD" || -f "$HISAT2_BUILD" ]] || die "找不到 hisat2-build：$HISAT2_BUILD"
[[ -x "$SAMTOOLS" || -f "$SAMTOOLS" ]] || die "找不到 samtools：$SAMTOOLS"
[[ -d "$RNA_DIR" ]] || die "缺少 RNA 目录：$RNA_DIR"

log "SCRIPT_DIR=$SCRIPT_DIR"
log "RNA_DIR=$RNA_DIR"
log "STAGE=$STAGE THREADS=$THREADS SORT_THREADS=$SORT_THREADS FORCE_RERUN=$FORCE_RERUN"
log "CONDA_ALIGNER=$CONDA_ALIGNER"
log "SUMMARY_TSV=$SUMMARY_TSV"

case "$STAGE" in
    index|align|all)
        for SAMPLE in "${SAMPLES[@]}"; do
            log "============================================================"
            log "样本：$SAMPLE"
            log "============================================================"
            case "$STAGE" in
                index) run_index "$SAMPLE" ;;
                align) run_align "$SAMPLE" ;;
                all)   run_index "$SAMPLE"; run_align "$SAMPLE" ;;
            esac
        done
        ;;
    summary) ;;
    *) die "未知 STAGE=$STAGE（index|align|all|summary）" ;;
esac

write_summary
print_summary
log "========== 结束 =========="
