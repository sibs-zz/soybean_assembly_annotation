#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# TIR_fixed / EDTA：最终 Chr.fa → 软链 → 容器 EDTA → conda RepeatMasker
#
# 输入优先：
#   07_polish_rename/${SAMPLE}/03_final/${SAMPLE}_Chr.fa
#
# 步骤：
#   1) 将 Chr.fa 软链到本脚本目录（annot/）
#   2) 在 ${sp}_TE/ 下用 apptainer 跑 EDTA（必须容器）
#   3) EDTA 完成后，在 EDTA conda 环境跑 RepeatMasker 软屏蔽
#
# 用法：
#   bash 01_EDTA.sh                   # 默认 STAGE=all：EDTA 前台完成后立刻 RM
#   STAGE=edta bash 01_EDTA.sh        # 仅后台启动 EDTA
#   RUN_RM=1 bash 01_EDTA.sh          # 仅跑 RepeatMasker（需 TElib 已生成）
#   FORCE_RERUN=1 bash 01_EDTA.sh     # 清理旧 TIR 中间结果后重跑 EDTA
# ============================================================

# ======================== 配置区 =============================

SAMPLES=(
    "KD35_primary"
#    "N2L_primary"
)

ROOT="/tiandata2/zzh/sjy_assembly"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

THREADS="${THREADS:-32}"
RM_THREADS="${RM_THREADS:-64}"
FORCE_RERUN="${FORCE_RERUN:-0}"
# edta | rm | all（默认全跑：EDTA 前台完成后立刻 RM）
STAGE="${STAGE:-all}"
# 兼容旧习惯：RUN_RM=1 等价于 STAGE=rm
if [[ "${RUN_RM:-0}" == "1" ]]; then
    STAGE="rm"
fi

SIF="${SIF:-/tiandata2/zzh/tools/EDTA_rebuild/EDTA_2.3.0_tirfix.sif}"
APPTAINER="${APPTAINER:-apptainer}"
CONDA_EDTA="${CONDA_EDTA:-/home/zzh/anaconda3/envs/EDTA}"

# ============================================================

timestamp() { date "+%F %T"; }
log() { echo "[$(timestamp)] $*"; }
die() { echo "[$(timestamp)] ERROR: $*" >&2; exit 1; }

resolve_chr_fa() {
    local sample="$1"
    local fa07="${ROOT}/07_polish_rename/${sample}/03_final/${sample}_Chr.fa"
    if [[ -s "$fa07" ]]; then
        echo "$fa07"
        return 0
    fi
    echo ""
}

link_chr_fa() {
    local sample="$1" src="$2"
    local sp="${sample%%_*}"
    local link_annot="${SCRIPT_DIR}/${sample}_Chr.fa"
    local te_dir="${SCRIPT_DIR}/${sp}_TE"
    local link_te="${te_dir}/${sp}_Chr.fa"

    mkdir -p "$te_dir"
    ln -sfn "$src" "$link_annot"
    ln -sfn "$src" "$link_te"
    log "软链：$link_annot -> $src"
    log "软链：$link_te -> $src"
}

clean_tir_intermediates() {
    local sp="$1"
    local te_dir="${SCRIPT_DIR}/${sp}_TE"
    local genome="${sp}_Chr.fa"
    (
        cd "$te_dir"
        rm -rf TIR-Learner-Result
        rm -rf "${genome}.mod.EDTA.raw/"*TIR* 2>/dev/null || true
        rm -f ./*TIR* 2>/dev/null || true
    )
}

run_edta() {
    local sample="$1"
    local sp="${sample%%_*}"
    local te_dir="${SCRIPT_DIR}/${sp}_TE"
    local genome="${sp}_Chr.fa"
    local logf="${te_dir}/edta_${sp}.log"
    local bg="${2:-1}"

    [[ -L "${te_dir}/${genome}" || -s "${te_dir}/${genome}" ]] \
        || die "缺少基因组软链：${te_dir}/${genome}"

    if [[ "$FORCE_RERUN" == "1" ]]; then
        log "FORCE_RERUN=1：清理 ${sp} 旧 TIR 中间结果"
        clean_tir_intermediates "$sp"
    fi

    if [[ -s "${te_dir}/${genome}.mod.EDTA.TElib.fa" && "$FORCE_RERUN" != "1" ]]; then
        log "已有 TElib，跳过 EDTA：${te_dir}/${genome}.mod.EDTA.TElib.fa"
        return 0
    fi

    log "EDTA（容器）启动：$sample  日志：$logf"
    if [[ "$bg" == "1" ]]; then
        (
            cd "$te_dir"
            nohup "$APPTAINER" exec --cleanenv -B /tiandata2:/tiandata2 \
                "$SIF" \
                EDTA.pl --genome "$genome" --species others --step all --threads "$THREADS" \
                >"$logf" 2>&1
        ) &
        log "EDTA 已推入后台 PID=$!  sample=$sample"
    else
        (
            cd "$te_dir"
            "$APPTAINER" exec --cleanenv -B /tiandata2:/tiandata2 \
                "$SIF" \
                EDTA.pl --genome "$genome" --species others --step all --threads "$THREADS" \
                >"$logf" 2>&1
        )
        log "EDTA 前台完成：$sample"
    fi
}

run_repeatmasker() {
    local sample="$1"
    local sp="${sample%%_*}"
    local te_dir="${SCRIPT_DIR}/${sp}_TE"
    local genome="${sp}_Chr.fa"
    local lib="${te_dir}/${genome}.mod.EDTA.TElib.fa"
    local logf="${te_dir}/repeatmasker_${sp}.log"
    local rm_bin="${CONDA_EDTA}/bin/RepeatMasker"

    [[ -s "$lib" ]] || die "缺少 TElib（EDTA 未完成？）：$lib"
    [[ -x "$rm_bin" || -f "$rm_bin" ]] || die "EDTA 环境中无 RepeatMasker：$rm_bin"
    [[ -L "${te_dir}/${genome}" || -s "${te_dir}/${genome}" ]] \
        || die "缺少基因组：$te_dir/$genome"

    if [[ -s "${te_dir}/${genome}.masked" && "$FORCE_RERUN" != "1" ]]; then
        log "已有 soft-masked，跳过 RepeatMasker：${te_dir}/${genome}.masked"
        return 0
    fi

    log "RepeatMasker（EDTA conda）启动：$sample  日志：$logf"
    (
        cd "$te_dir"
        env PATH="${CONDA_EDTA}/bin:${PATH}" \
            CONDA_PREFIX="$CONDA_EDTA" \
            "$rm_bin" -pa "$RM_THREADS" -q -xsmall \
            -lib "$lib" \
            "$genome" \
            >"$logf" 2>&1
    )
    log "RepeatMasker 完成：$sample -> ${te_dir}/${genome}.masked / ${genome}.out"
}

command -v "$APPTAINER" >/dev/null 2>&1 || die "找不到 apptainer：$APPTAINER"
[[ -s "$SIF" ]] || die "缺少 EDTA sif：$SIF"
[[ -d "$CONDA_EDTA" ]] || die "缺少 EDTA conda：$CONDA_EDTA"

log "SCRIPT_DIR=$SCRIPT_DIR"
log "STAGE=$STAGE THREADS=$THREADS RM_THREADS=$RM_THREADS FORCE_RERUN=$FORCE_RERUN"
log "SIF=$SIF"
log "CONDA_EDTA=$CONDA_EDTA"

for SAMPLE in "${SAMPLES[@]}"; do
    sp="${SAMPLE%%_*}"
    SRC_FA="$(resolve_chr_fa "$SAMPLE")"

    log "============================================================"
    log "样本：$SAMPLE"
    log "============================================================"

    if [[ -z "$SRC_FA" ]]; then
        log "跳过 $SAMPLE：找不到 Chr.fa（需先完成 07_polish_rename）"
        continue
    fi
    log "Chr.fa=$SRC_FA"

    link_chr_fa "$SAMPLE" "$SRC_FA"

    case "$STAGE" in
        edta)
            run_edta "$SAMPLE" 1
            ;;
        rm)
            run_repeatmasker "$SAMPLE"
            ;;
        all)
            run_edta "$SAMPLE" 0
            run_repeatmasker "$SAMPLE"
            ;;
        *)
            die "未知 STAGE=$STAGE（edta|rm|all）"
            ;;
    esac
done

echo "------------------------------------------------"
case "$STAGE" in
    edta)
        log "EDTA 已按样本推入后台（约 1–3 天）。监控："
        for SAMPLE in "${SAMPLES[@]}"; do
            sp="${SAMPLE%%_*}"
            echo "  tail -f ${SCRIPT_DIR}/${sp}_TE/edta_${sp}.log"
        done
        log "EDTA 全部完成后跑 RepeatMasker："
        echo "  RUN_RM=1 bash ${SCRIPT_DIR}/01_EDTA.sh"
        ;;
    rm)
        log "RepeatMasker 阶段结束。"
        ;;
    all)
        log "EDTA + RepeatMasker 全部完成。"
        ;;
esac
log "========== 结束 =========="
