#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# BRAKER3：软屏蔽基因组 + Atlas BAM + 蛋白数据库 → 基因注释
#
# 依赖（均来自前序脚本）：
#   基因组  01_EDTA.sh  →  ${sp}_TE/${sp}_Chr.fa.masked
#   BAM     02_HISAT2   →  ${sp}_RNA/${sp}_RNA_sorted.bam
#
# 输出：
#   BRAKER_${sp}/braker.gff3   BRAKER3 注释结果
#   braker_${sp}.log           运行日志
#
# 用法：
#   bash 03_run_BRAKER3.sh                # 后台并行跑所有样本
#   FORCE_RERUN=1 bash 03_run_BRAKER3.sh  # 清理旧结果后重跑
# ============================================================

# ======================== 配置区 =============================

SAMPLES=(
    "KD35_primary"
    "N2L_primary"
)

ROOT="/tiandata2/zzh/sjy_assembly"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

THREADS="${THREADS:-32}"
FORCE_RERUN="${FORCE_RERUN:-0}"

CONDA_BRAKER="${CONDA_BRAKER:-/home/zzh/anaconda3/envs/braker3}"
BRAKER="${BRAKER:-${CONDA_BRAKER}/bin/braker.pl}"
GENEMARK_PATH="${GENEMARK_PATH:-/tiandata2/zzh/tools/GeneMark-ETP-main/bin}"
PROTHINT_PATH="${PROTHINT_PATH:-/tiandata2/zzh/tools/ProtHint/bin}"
PROT_DB="${PROT_DB:-${ROOT}/anno/Viridiplantae.fa}"

# ============================================================

timestamp() { date "+%F %T"; }
log()       { echo "[$(timestamp)] $*"; }
die()       { echo "[$(timestamp)] ERROR: $*" >&2; exit 1; }

# 激活 braker3 conda 并清理 Perl 环境变量（避免污染）
source /home/zzh/anaconda3/bin/activate braker3
unset PERL5LIB PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT
export PATH="${CONDA_BRAKER}/bin:${PATH}"
hash -r

# 依赖检查
log "依赖检查..."
log "  perl         : $(which perl)"
log "  braker.pl    : $(which braker.pl)"
log "  gffread      : $(which gffread)"
log "  GENEMARK_PATH: $GENEMARK_PATH"
log "  PROTHINT_PATH: $PROTHINT_PATH"

[[ -f "${GENEMARK_PATH}/gmetp.pl" ]]   || die "缺少 gmetp.pl：${GENEMARK_PATH}/gmetp.pl"
[[ -f "${PROTHINT_PATH}/prothint.py" ]] || die "缺少 prothint.py：${PROTHINT_PATH}/prothint.py"
[[ -s "$PROT_DB" ]]                     || die "缺少蛋白数据库：$PROT_DB"
which gffread >/dev/null                || die "找不到 gffread"

declare -a PIDS=()

for SAMPLE in "${SAMPLES[@]}"; do
    sp="${SAMPLE%%_*}"
    sp_lower="${sp,,}"

    genome="${SCRIPT_DIR}/${sp}_TE/${sp}_Chr.fa.masked"
    bam="${SCRIPT_DIR}/${sp}_RNA/${sp}_RNA_sorted.bam"
    workdir="${SCRIPT_DIR}/BRAKER_${sp}"
    logf="${SCRIPT_DIR}/braker_${sp}.log"

    log "============================================================"
    log "样本：$SAMPLE"
    log "============================================================"

    [[ -s "$genome" ]] || { log "跳过 $SAMPLE：缺少软屏蔽基因组 $genome"; continue; }
    [[ -s "$bam" ]]    || { log "跳过 $SAMPLE：缺少 Atlas BAM $bam"; continue; }

    species_dir="${CONDA_BRAKER}/config/species/soybean_${sp_lower}"

    if [[ -s "${workdir}/braker.gff3" && "$FORCE_RERUN" != "1" ]]; then
        log "已有 braker.gff3，跳过（FORCE_RERUN=1 可重跑）：$workdir"
        continue
    fi
    # 空目录/失败残留会挡住下次启动；FORCE_RERUN 或未完成结果都清掉
    if [[ -d "$workdir" ]]; then
        log "清理未完成工作目录：$workdir"
        rm -rf "$workdir"
    fi
    # AUGUSTUS 物种参数目录一旦存在，braker.pl 会直接退出（除非 --useexisting）
    if [[ -d "$species_dir" ]]; then
        log "清理旧 AUGUSTUS species：$species_dir"
        rm -rf "$species_dir"
    fi

    log "BRAKER3 后台启动：$SAMPLE"
    log "  genome  : $genome"
    log "  bam     : $bam"
    log "  workdir : $workdir"
    log "  log     : $logf"

    nohup "$BRAKER" \
        --genome="$genome" \
        --bam="$bam" \
        --prot_seq="$PROT_DB" \
        --softmasking \
        --species="soybean_${sp_lower}" \
        --gff3 \
        --threads="$THREADS" \
        --GENEMARK_PATH="$GENEMARK_PATH" \
        --PROTHINT_PATH="$PROTHINT_PATH" \
        --workingdir="$workdir" \
        >"$logf" 2>&1 &

    PIDS+=("$!")
    log "BRAKER3 已推入后台 PID=${PIDS[-1]}  sample=$SAMPLE"
done

echo "------------------------------------------------"
if [[ ${#PIDS[@]} -gt 0 ]]; then
    log "BRAKER3 已按样本推入后台（约数小时）。监控："
    for SAMPLE in "${SAMPLES[@]}"; do
        sp="${SAMPLE%%_*}"
        echo "  tail -f ${SCRIPT_DIR}/braker_${sp}.log"
    done
else
    log "无样本启动（所有样本已跳过或缺少输入）。"
fi
log "========== 结束 =========="
