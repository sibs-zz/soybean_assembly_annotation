#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# 04  转录本组装（短读，对标 ZH13 的 RNA evidence）
#
#   StringTie  : HISAT2 BAM（aligner env）
#   Trinity-GG : genome-guided Trinity（annot_trinity env）
#
# 输入：${sp}_RNA/${sp}_RNA_sorted.bam
# 输出：
#   ${sp}_models/04_transcripts/${sp}.stringtie.gtf
#   ${sp}_models/04_transcripts/${sp}.stringtie.fa
#   ${sp}_models/04_transcripts/trinity_gg/Trinity-GG.fasta
#   ${sp}_models/04_transcripts/${sp}.transcripts.fa   # 合并去冗余，给 PASA
#
# 不依赖 03 BRAKER。BAM 齐了就能跑。
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SAMPLES=(
    "KD35_primary"
    "N2L_primary"
)
ROOT="/tiandata2/zzh/sjy_assembly"
# shellcheck source=/dev/null
source "${SCRIPT_DIR}/annot_common.sh"

THREADS="${THREADS:-16}"
FORCE_RERUN="${FORCE_RERUN:-0}"
FILTER_SAMPLE="${SAMPLE:-}"
RUN_STRINGTIE="${RUN_STRINGTIE:-1}"
RUN_TRINITY="${RUN_TRINITY:-1}"
TRINITY_MAX_MEM="${TRINITY_MAX_MEM:-64G}"
TRINITY_MAX_INTRON="${TRINITY_MAX_INTRON:-20000}"
STRINGTIE_MIN_FRAC="${STRINGTIE_MIN_FRAC:-0.01}"
STRINGTIE_MIN_COV="${STRINGTIE_MIN_COV:-1.0}"
STRINGTIE_MIN_LEN="${STRINGTIE_MIN_LEN:-150}"

[[ -x "$STRINGTIE" ]] || die "找不到 stringtie：$STRINGTIE（aligner env）"
[[ -x "$GFFREAD" ]] || die "找不到 gffread：$GFFREAD"

log "THREADS=$THREADS RUN_STRINGTIE=$RUN_STRINGTIE RUN_TRINITY=$RUN_TRINITY"

for SAMPLE in "${SAMPLES[@]}"; do
    [[ -z "$FILTER_SAMPLE" || "$SAMPLE" == "$FILTER_SAMPLE" ]] || continue
    sp="${SAMPLE%%_*}"
    bam="$(resolve_bam "$sp")"
    outdir="${SCRIPT_DIR}/${sp}_models/04_transcripts"
    st_gtf="${outdir}/${sp}.stringtie.gtf"
    st_fa="${outdir}/${sp}.stringtie.fa"
    tr_dir="${outdir}/trinity_gg"
    tr_fa="${tr_dir}/Trinity-GG.fasta"
    all_fa="${outdir}/${sp}.transcripts.fa"
    logf_st="${outdir}/${sp}.stringtie.log"
    logf_tr="${outdir}/${sp}.trinity.log"

    log "============================================================"
    log "样本：$SAMPLE"
    log "============================================================"
    [[ -s "$bam" ]] || { log "跳过 $SAMPLE：缺少 BAM $bam"; continue; }
    mkdir -p "$outdir"

    if [[ "$RUN_STRINGTIE" == "1" ]]; then
        if [[ -s "$st_gtf" && "$FORCE_RERUN" != "1" ]]; then
            log "已有 StringTie GTF，跳过：$st_gtf"
        else
            use_env aligner
            log "StringTie -> $st_gtf"
            "$STRINGTIE" "$bam" -o "$st_gtf" -p "$THREADS" \
                -m "$STRINGTIE_MIN_LEN" -c "$STRINGTIE_MIN_COV" -f "$STRINGTIE_MIN_FRAC" \
                >"$logf_st" 2>&1
            [[ -s "$st_gtf" ]] || die "StringTie 失败：见 $logf_st"
        fi
        if [[ -s "$st_gtf" ]]; then
            genome="$(resolve_genome "$SAMPLE")"
            "$GFFREAD" "$st_gtf" -g "$genome" -w "$st_fa"
        fi
    fi

    if [[ "$RUN_TRINITY" == "1" ]]; then
        link_trinity_gg_fasta "$tr_dir" || true
        if [[ -s "$tr_fa" && "$FORCE_RERUN" != "1" ]]; then
            log "已有 Trinity-GG，跳过：$tr_fa"
        else
            use_env annot_trinity
            require_bin Trinity "annot_trinity env"
            require_salmon_1x
            mkdir -p "$tr_dir"
            log "Trinity genome-guided -> $tr_dir  (max_memory=$TRINITY_MAX_MEM)"
            Trinity --genome_guided_bam "$bam" \
                --genome_guided_max_intron "$TRINITY_MAX_INTRON" \
                --max_memory "$TRINITY_MAX_MEM" \
                --CPU "$THREADS" \
                --output "$tr_dir" \
                --full_cleanup \
                >"$logf_tr" 2>&1
            link_trinity_gg_fasta "$tr_dir" || true
            if [[ ! -s "$tr_fa" ]]; then
                alt=$(find "$tr_dir" -name 'Trinity-GG.fasta' -o -name 'Trinity.fasta' | head -1 || true)
                [[ -n "${alt:-}" && -s "$alt" ]] && ln -sfn "$alt" "$tr_fa"
            fi
            [[ -s "$tr_fa" ]] || die "Trinity 失败：见 $logf_tr"
        fi
    fi

    # 合并转录本给 PASA
    tmp="${outdir}/.cat.fa"
    : >"$tmp"
    [[ -s "$st_fa" ]] && cat "$st_fa" >>"$tmp"
    [[ -s "$tr_fa" ]] && cat "$tr_fa" >>"$tmp"
    if [[ ! -s "$tmp" ]]; then
        log "跳过合并：还没有 StringTie/Trinity fasta"
        continue
    fi
    if [[ -x "$SEQKIT" ]]; then
        "$SEQKIT" seq -i -w 0 "$tmp" >"${tmp}.clean"
        mv "${tmp}.clean" "$tmp"
    fi
    if [[ -x "$CDHIT" ]]; then
        log "cd-hit-est 95% -> $all_fa"
        use_env EDTA
        "$CDHIT" -i "$tmp" -o "$all_fa" -c 0.95 -n 10 -T "$THREADS" -M 0 \
            >"${outdir}/${sp}.cdhit.log" 2>&1 || cp "$tmp" "$all_fa"
    else
        cp "$tmp" "$all_fa"
    fi
    rm -f "$tmp"
    n=$(grep -c '^>' "$all_fa" || true)
    log "完成 $SAMPLE  transcripts=$n  $all_fa"
done

log "========== 04 transcripts 结束 =========="
log "下一步：bash 05_run_mikado.sh"
