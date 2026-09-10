#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# 06  PASA alignment assembly（SQLite，ALIGNERS=gmap）
#
# 环境：annot_pasa（自带 GMAP/BLAT/SQLite）
#   export PASAHOME=.../opt/pasa-2.5.3
#
# 输入：04 ${sp}.transcripts.fa（已含 StringTie+Trinity）、基因组
# 不要再传 --trans_gtf：StringTie ID 已在 -t 里，再导会 UNIQUE cdna_acc 炸库
# 输出：
#   ${sp}_models/06_pasa/${sp}.pasa.sqlite
#   ${sp}_models/06_pasa/*pasa_assemblies*.gff3
#   ${sp}_models/06_pasa/*.assemblies.fasta
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
MAX_INTRON="${MAX_INTRON:-20000}"

[[ -d "$PASAHOME" ]] || die "PASAHOME 不存在：$PASAHOME"
use_env annot_pasa
export PASAHOME
require_bin Launch_PASA_pipeline.pl "annot_pasa；并确认 PASAHOME=$PASAHOME"

log "PASAHOME=$PASAHOME THREADS=$THREADS"

for SAMPLE in "${SAMPLES[@]}"; do
    [[ -z "$FILTER_SAMPLE" || "$SAMPLE" == "$FILTER_SAMPLE" ]] || continue
    sp="${SAMPLE%%_*}"
    genome="$(resolve_genome "$SAMPLE")"
    d04="${SCRIPT_DIR}/${sp}_models/04_transcripts"
    tx_fa="${d04}/${sp}.transcripts.fa"
    outdir="${SCRIPT_DIR}/${sp}_models/06_pasa"
    sqlite="${outdir}/${sp}.pasa.sqlite"
    align_conf="${outdir}/alignAssembly.config"
    done_gff="${outdir}/${sp}.pasa_assemblies.gff3"

    log "============================================================"
    log "样本：$SAMPLE"
    log "============================================================"
    [[ -s "$tx_fa" ]] || { log "跳过 $SAMPLE：缺少 $tx_fa（先跑 04）"; continue; }

    if [[ -s "$done_gff" && "$FORCE_RERUN" != "1" ]]; then
        log "已有 PASA assemblies GFF，跳过：$done_gff"
        continue
    fi
    mkdir -p "$outdir"
    if [[ "$FORCE_RERUN" == "1" ]]; then
        rm -f "$sqlite" "${sqlite}-"* "$done_gff"
        rm -rf "${outdir}/__pasa_${sp}.pasa.sqlite_SQLite_chkpts"
        rm -f "${outdir}/__pasa_${sp}.pasa.sqlite_SQLite_chkpts.cmds_log"
    fi

    cat >"$align_conf" <<EOF
DATABASE=${sqlite}
validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=75
validate_alignments_in_db.dbi:--MIN_AVG_PER_ID=80
subcluster_builder.dbi:-m=50
EOF

    log "PASA align/assemble -> $outdir  （-t 已含 StringTie，不用 --trans_gtf）"
    pasa_cmd=(
        Launch_PASA_pipeline.pl
        -c "$align_conf"
        -C -R
        -g "$genome"
        -t "$tx_fa"
        --ALIGNERS gmap
        --MAX_INTRON_LENGTH "$MAX_INTRON"
        --CPU "$THREADS"
        --TRANSDECODER
    )
    (
        cd "$outdir"
        "${pasa_cmd[@]}"
    ) >>"${outdir}/${sp}.pasa_align.log" 2>&1

    # 收集 assemblies gff
    found=$(find "$outdir" -name '*pasa_assemblies*.gff3' | head -1 || true)
    if [[ -z "${found:-}" ]]; then
        found=$(find "$outdir" -name '*assemblies*.gff3' | head -1 || true)
    fi
    [[ -s "${found:-}" ]] || die "PASA 未产出 assemblies GFF3：见 ${outdir}/${sp}.pasa_align.log"
    cp "$found" "$done_gff"
    log "完成 $SAMPLE  $done_gff"
done

log "========== 06 PASA align 结束 =========="
log "下一步（等 03 BRAKER 出预测）：bash 07_run_evm.sh"
