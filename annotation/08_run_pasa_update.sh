#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# 08  PASA annotation compare：在 EVM 模型上加 UTR / 异构体
#     （ZH13 发表注释的典型后处理；无 Iso-seq，TRUST_FL_STATUS=0）
#
# 环境：annot_pasa    PASAHOME 必须已 export
# 输入：07 evm.gff3、06 sqlite + transcripts.fa
# 输出：
#   ${sp}_models/${sp}.models.gff3              # 给 09 polish
#   ${sp}_models/${sp}.models.nuclear.gff3      # Chr01–20
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

[[ -d "$PASAHOME" ]] || die "PASAHOME 不存在：$PASAHOME"
use_env annot_pasa
export PASAHOME
require_bin Launch_PASA_pipeline.pl "annot_pasa"

log "PASAHOME=$PASAHOME"

for SAMPLE in "${SAMPLES[@]}"; do
    [[ -z "$FILTER_SAMPLE" || "$SAMPLE" == "$FILTER_SAMPLE" ]] || continue
    sp="${SAMPLE%%_*}"
    genome="$(resolve_genome "$SAMPLE")"
    tx_fa="${SCRIPT_DIR}/${sp}_models/04_transcripts/${sp}.transcripts.fa"
    sqlite="${SCRIPT_DIR}/${sp}_models/06_pasa/${sp}.pasa.sqlite"
    evm_gff="${SCRIPT_DIR}/${sp}_models/07_evm/${sp}.evm.gff3"
    outdir="${SCRIPT_DIR}/${sp}_models/08_pasa_update"
    models="${SCRIPT_DIR}/${sp}_models/${sp}.models.gff3"
    nuclear="${SCRIPT_DIR}/${sp}_models/${sp}.models.nuclear.gff3"
    cmp_conf="${outdir}/annotCompare.config"

    log "============================================================"
    log "样本：$SAMPLE"
    log "============================================================"
    [[ -s "$evm_gff" ]] || { log "跳过 $SAMPLE：缺少 $evm_gff"; continue; }
    [[ -s "$tx_fa" ]] || { log "跳过 $SAMPLE：缺少 $tx_fa"; continue; }
    [[ -s "$sqlite" ]] || { log "跳过 $SAMPLE：缺少 PASA sqlite $sqlite（先跑 06）"; continue; }

    if [[ -s "$models" && "$FORCE_RERUN" != "1" ]]; then
        log "已有 models GFF3，跳过：$models"
        continue
    fi
    mkdir -p "$outdir"

    cat >"$cmp_conf" <<EOF
DATABASE=${sqlite}
cDNA_annotation_comparer.dbi:--MIN_PERCENT_OVERLAP=50
cDNA_annotation_comparer.dbi:--MIN_PERCENT_PROT_CODING=40
cDNA_annotation_comparer.dbi:--MIN_PERID_PROT_COMPARE=70
cDNA_annotation_comparer.dbi:--MIN_PERCENT_LENGTH_FL_COMPARE=70
cDNA_annotation_comparer.dbi:--MIN_PERCENT_LENGTH_NONFL_COMPARE=70
cDNA_annotation_comparer.dbi:--MIN_FL_ORF_SIZE=90
cDNA_annotation_comparer.dbi:--MIN_PERCENT_ALIGN_LENGTH=70
cDNA_annotation_comparer.dbi:--MIN_PERCENT_OVERLAP_GENE_REPLACE=80
cDNA_annotation_comparer.dbi:--STOMP_HIGH_PERCENTAGE_OVERLAPPING_GENE=1
cDNA_annotation_comparer.dbi:--TRUST_FL_STATUS=0
cDNA_annotation_comparer.dbi:--MAX_UTR_EXONS=3
EOF

    log "PASA annot_compare -> $outdir"
    (
        cd "$outdir"
        Launch_PASA_pipeline.pl \
            -c "$cmp_conf" \
            -A -L \
            -g "$genome" \
            -t "$tx_fa" \
            --annots "$evm_gff" \
            --CPU "$THREADS"
    ) >"${outdir}/${sp}.pasa_update.log" 2>&1

    found=$(find "$outdir" "$SCRIPT_DIR/${sp}_models/06_pasa" -name '*gene_structures_post_PASA_updates*.gff3' 2>/dev/null | sort | tail -1 || true)
    if [[ -z "${found:-}" ]]; then
        found=$(find "$outdir" -name '*gff3' | sort | tail -1 || true)
    fi
    [[ -s "${found:-}" ]] || die "PASA update 未产出 GFF3：见 ${outdir}/${sp}.pasa_update.log"

    python3 "${SCRIPT_DIR}/pack_structural_gff.py" \
        -i "$found" \
        -o "$models" \
        --nuclear-out "$nuclear"

    [[ -s "$models" ]] || die "未写出 $models"
    log "完成 $SAMPLE"
    log "  models  : $models"
    log "  nuclear : $nuclear"
done

log "========== 08 PASA update 结束 =========="
log "下一步收尾：bash 09_run_polish.sh"
