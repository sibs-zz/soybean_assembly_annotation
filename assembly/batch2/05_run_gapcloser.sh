#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# LR_Gapcloser 填缝（ONT）
#
# 输入：
#   04_t2t_output/${sample}_T2T.FINAL.fa
#   ont/${sp}_ONT_clean.fa
#
# 输出：
#   05_LRgapcloser_out/${sample}/${sample}_gapclosed.fa   # 仅保留最终轮
#
# 说明：
#   LR_Gapcloser 默认 -r 3，会生成 iteration-1/2/3。
#   跑完后删除 iteration-1、iteration-2 以及 tmp/db 等中间文件，
#   只留下最终 gapclosed 序列（对应第 3 轮，或提前收敛时的最后一轮）。
# ============================================================

# ======================== 配置区 =============================

LR_GAPCLOSER="/tiandata2/zzh/tools/LR_Gapcloser/src/LR_Gapcloser.sh"

ROOT="${ROOT:-/tiandata2/zzh/sjy_assembly/batch2}"

SAMPLES=(
    "BZ_primary"
    "DY11_primary"
    "TH2_primary"
    "TZX1998_primary"
    "TZX253_primary"
)

SCAFF_DIR="${ROOT}/04_t2t_output"
ONT_DIR="${ROOT}/ont"
OUT_DIR="${ROOT}/05_LRgapcloser_out"

THREADS=64

# nanopore
PLATFORM="n"

# 迭代轮数（正式结果用 3；只要最终一轮）
ITERATIONS=3

# 0：已有最终 fa 则跳过
# 1：删旧目录重跑
# 用法：FORCE_RERUN=1 bash 05_run_gapcloser.sh
FORCE_RERUN="${FORCE_RERUN:-0}"

# ============================================================

timestamp() { date "+%F %T"; }
log() { echo "[$(timestamp)] $*"; }
die() { echo "[$(timestamp)] ERROR: $*" >&2; exit 1; }

[[ -x "$LR_GAPCLOSER" || -f "$LR_GAPCLOSER" ]] ||
    die "找不到 LR_Gapcloser：$LR_GAPCLOSER"

mkdir -p "$OUT_DIR"

for sample in "${SAMPLES[@]}"; do
    log "============================================================"
    log "开始处理：$sample"
    log "============================================================"

    sp="${sample%%_*}"
    scaff="${SCAFF_DIR}/${sample}_T2T.FINAL.fa"
    reads="${ONT_DIR}/${sp}_ONT_clean.fa"

    sample_out="${OUT_DIR}/${sample}"
    final_fa="${sample_out}/${sample}_gapclosed.fa"
    work_dir="${sample_out}/work"   # LR 要求 -o 目录事先不存在
    pipe_log="${sample_out}/pipe.log"

    [[ -s "$scaff" ]] || die "缺少 scaffold：$scaff"
    [[ -s "$reads" ]] || die "缺少 ONT：$reads"

    if [[ -s "$final_fa" && "$FORCE_RERUN" -eq 0 ]]; then
        log "已有最终结果，跳过：$final_fa"
        ls -lh "$final_fa"
        continue
    fi

    if [[ "$FORCE_RERUN" -eq 1 ]]; then
        log "FORCE_RERUN=1，清理旧目录：$sample_out"
        rm -rf "$sample_out"
    fi

    # LR_Gapcloser：若 -o 目录已存在会直接退出
    mkdir -p "$sample_out"
    rm -rf "$work_dir"

    log "scaffold : $scaff"
    log "ONT      : $reads"
    log "work dir : $work_dir"
    log "threads  : $THREADS ; platform=$PLATFORM ; -r $ITERATIONS"

    (
        bash "$LR_GAPCLOSER" \
            -i "$scaff" \
            -l "$reads" \
            -s "$PLATFORM" \
            -t "$THREADS" \
            -r "$ITERATIONS" \
            -o "$work_dir"
    ) >"$pipe_log" 2>&1

    # 取最高一轮 iteration-*/gapclosed.fasta（通常是 iteration-3）
    last_gapclosed="$(
        find "$work_dir" -maxdepth 2 -type f -path '*/iteration-*/gapclosed.fasta' \
            | sort -V \
            | tail -n 1
    )"

    [[ -n "$last_gapclosed" && -s "$last_gapclosed" ]] ||
        die "$sample 未找到 gapclosed.fasta，见 $pipe_log"

    iter_name="$(basename "$(dirname "$last_gapclosed")")"
    log "最终轮次：$iter_name -> $last_gapclosed"

    cp -f "$last_gapclosed" "$final_fa"

    # 只保留最终 fa + 日志；删掉 iteration-1/2/3、tmp、LR 中间文件
    log "清理中间文件（仅保留最终 fa 与 pipe.log）..."
    rm -rf "$work_dir"

    # 顺便丢掉可能散落的中间产物
    find "$sample_out" -mindepth 1 -maxdepth 1 ! -name "$(basename "$final_fa")" ! -name 'pipe.log' -exec rm -rf {} +

    log "$sample 完成："
    ls -lh "$final_fa" "$pipe_log"
    echo "------------------------------------------------"
done

log "========== 全部 LR_Gapcloser 任务结束 =========="
