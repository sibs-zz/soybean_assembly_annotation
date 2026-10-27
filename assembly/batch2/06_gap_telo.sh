#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# 06_gap_telo：端粒/gap 审计 + ONT 自动补端粒 + 局部填 gap
#
# 方案 A（全自动）：高置信 ONT 命中直接写入最终 fa，并保留记录。
#
# 用法：
#   bash 06_gap_telo.sh
#   FORCE_RERUN=1 bash 06_gap_telo.sh   # 清除中间结果重跑
#
# 输出：
#   06_gap_telo/${SAMPLE}/
#     01_audit/          端粒+gap 统计
#     02_telo_fill/      端粒补全中间文件与记录
#     03_gap_fill/       gap 填补中间文件与记录
#     04_final/          最终 fa + 汇总记录
# ============================================================

# ======================== 配置区 =============================

SAMPLES=(
    "BZ_primary"
    "DY11_primary"
    "TH2_primary"
    "TZX1998_primary"
    "TZX253_primary"
)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${ROOT:-/tiandata2/zzh/sjy_assembly/batch2}"
TOP_N=20
THREADS="${THREADS:-64}"
FORCE_RERUN="${FORCE_RERUN:-1}"

# 放宽参数：尽量多补端粒、多填 gap
BAIT_LEN=100000
GAP_FLANK=80000
USE_FULL_ONT=1

# apply-telo
TELO_MIN_ALEN=500
TELO_MIN_IDENT=0.65
TELO_MIN_COPIES=8
TELO_EDGE_BP=8000
TELO_MIN_OVERHANG=50

# apply-gap (unified bait)
GAP_MIN_ALEN=500
GAP_MIN_IDENT=0.60
GAP_MAX_FILL=5000

QUARTET="/tiandata2/zzh/tools/quarTeT/quartet.py"
CHECK_GAP="${SCRIPT_DIR}/../check_gap.sh"
HELPER="${SCRIPT_DIR}/../06_gap_telo_lib.py"
SEQKIT="${SEQKIT:-/home/zzh/anaconda3/envs/py37/bin/seqkit}"
MINIMAP2="${MINIMAP2:-minimap2}"

# quarTeT 依赖环境（TeloExplorer 等）
CONDA_ENV="quarTeTdependencies"

# ============================================================

timestamp() { date "+%F %T"; }
log() { echo "[$(timestamp)] $*"; }
die() { echo "[$(timestamp)] ERROR: $*" >&2; exit 1; }

run_telo_explorer() {
    local fa="$1" prefix="$2" outdir="$3"
    (
        cd "$outdir"
        $PYTHON3 "$QUARTET" TeloExplorer -i "$fa" -c plant -p "$prefix"
    )
}

run_telo_fill_round() {
    local round="$1" in_fa="$2" telo_info="$3" missing_tsv="$4" out_fa="$5" record="$6" workdir="$7"
    local bait_fa="${workdir}/telo_baits_r${round}.fa"
    local paf="${workdir}/telo_baits_r${round}.paf"
    local hit_ids="${workdir}/telo_hit_ids_r${round}.txt"
    local hit_fa="${workdir}/telo_hit_reads_r${round}.fa"

    if [[ ! -s "$missing_tsv" ]]; then
        log "端粒 round${round}：无缺失端，跳过"
        cp -f "$in_fa" "$out_fa"
        printf "scaffold\tside\tstatus\tread\tadded_bp\ttelo_copies\tident\tnote\n" >"$record"
        return 0
    fi

    log "端粒 round${round}：$(wc -l <"$missing_tsv") 个缺失端"
    $PYTHON3 "$HELPER" make-telo-baits \
        --fasta "$in_fa" \
        --missing-telo "$missing_tsv" \
        --bait-len "$BAIT_LEN" \
        --out-fa "$bait_fa"

    log "minimap2 map-ont：端粒 bait round${round} <- ONT..."
    "$MINIMAP2" -x map-ont -t "$THREADS" "$bait_fa" "$ONT_FA" >"$paf"

    awk '{print $1}' "$paf" | sort -u >"$hit_ids"
    "$SEQKIT" grep -f "$hit_ids" "$ONT_FA" >"$hit_fa"

    $PYTHON3 "$HELPER" apply-telo \
        --fasta "$in_fa" \
        --telo-info "$telo_info" \
        --top-list "${AUD}/top${TOP_N}.list" \
        --missing-telo "$missing_tsv" \
        --paf "$paf" \
        --reads-fa "$hit_fa" \
        --out-fa "$out_fa" \
        --record "$record" \
        --min-alen "$TELO_MIN_ALEN" \
        --min-ident "$TELO_MIN_IDENT" \
        --min-telo-copies "$TELO_MIN_COPIES" \
        --edge-bp "$TELO_EDGE_BP" \
        --min-overhang "$TELO_MIN_OVERHANG"
}

run_gap_fill_round() {
    local round="$1" mode="$2" in_fa="$3" gaps_tsv="$4" out_fa="$5" record="$6" workdir="$7"
    local bait_fa="${workdir}/gap_baits_r${round}.fa"
    local paf="${workdir}/gap_baits_r${round}.paf"
    local hit_ids="${workdir}/gap_hit_ids_r${round}.txt"
    local hit_fa="${workdir}/gap_hit_reads_r${round}.fa"
    local unified_flag=()
    local min_alen="$GAP_MIN_ALEN"
    local min_ident="$GAP_MIN_IDENT"

    local gap_n
    gap_n=$(($(wc -l <"$gaps_tsv") - 1))
    if [[ "$gap_n" -le 0 ]]; then
        log "gap round${round}：无 gap，跳过"
        cp -f "$in_fa" "$out_fa"
        printf "scaffold\tgap_index\told_start\told_end\told_length\tstatus\tread\tfill_bp\tident\tnote\n" >"$record"
        return 0
    fi

    if [[ "$mode" == "unified" ]]; then
        unified_flag=(--unified)
    fi
    if [[ "$round" -ge 2 ]]; then
        min_alen=300
        min_ident=0.55
    fi

    log "gap round${round} (${mode})：${gap_n} 个 gap"
    $PYTHON3 "$HELPER" make-gap-baits \
        --fasta "$in_fa" \
        --gaps-tsv "$gaps_tsv" \
        --flank "$GAP_FLANK" \
        --out-fa "$bait_fa" \
        "${unified_flag[@]}"

    log "minimap2 map-ont：gap bait round${round} <- ONT..."
    "$MINIMAP2" -x map-ont -t "$THREADS" "$bait_fa" "$ONT_FA" >"$paf"

    awk '{print $1}' "$paf" | sort -u >"$hit_ids"
    "$SEQKIT" grep -f "$hit_ids" "$ONT_FA" >"$hit_fa"

    $PYTHON3 "$HELPER" apply-gap \
        --fasta "$in_fa" \
        --gaps-tsv "$gaps_tsv" \
        --paf "$paf" \
        --reads-fa "$hit_fa" \
        --out-fa "$out_fa" \
        --record "$record" \
        --flank "$GAP_FLANK" \
        --min-alen "$min_alen" \
        --min-ident "$min_ident" \
        --max-fill "$GAP_MAX_FILL" \
        "${unified_flag[@]}"
}

# -------------------- 激活 conda 环境 --------------------
if [[ -z "${CONDA_DEFAULT_ENV:-}" || "${CONDA_DEFAULT_ENV}" != "$CONDA_ENV" ]]; then
    command -v conda >/dev/null 2>&1 || die "未找到 conda"
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$CONDA_ENV"
fi
log "conda env: ${CONDA_DEFAULT_ENV} ($(command -v python3))"

# 强制使用该环境下的 python，避免 PATH 里其他 python 抢先
PYTHON3="${CONDA_PREFIX}/bin/python3"
[[ -x "$PYTHON3" ]] || die "环境 ${CONDA_ENV} 中找不到 python3：$PYTHON3"

[[ -f "$HELPER" ]] || die "缺少 helper：$HELPER"
[[ -f "$CHECK_GAP" ]] || die "缺少 check_gap.sh"
[[ -x "$SEQKIT" || -f "$SEQKIT" ]] || die "缺少 seqkit：$SEQKIT"
command -v "$MINIMAP2" >/dev/null 2>&1 || die "缺少 minimap2"
[[ -f "$QUARTET" ]] || die "缺少 quarTeT：$QUARTET"

log "参数：BAIT_LEN=$BAIT_LEN GAP_FLANK=$GAP_FLANK FORCE_RERUN=$FORCE_RERUN USE_FULL_ONT=$USE_FULL_ONT"
log "端粒阈值：alen>=$TELO_MIN_ALEN ident>=$TELO_MIN_IDENT copies>=$TELO_MIN_COPIES"
log "gap阈值：alen>=$GAP_MIN_ALEN ident>=$GAP_MIN_IDENT max_fill=$GAP_MAX_FILL"

for SAMPLE in "${SAMPLES[@]}"; do
    sp="${SAMPLE%%_*}"
    IN_FA="${ROOT}/05_LRgapcloser_out/${SAMPLE}/${SAMPLE}_gapclosed.fa"
    ONT_FA="${ROOT}/ont/${sp}_ONT_clean.fa"
    if [[ "$USE_FULL_ONT" != "1" && -s "${ROOT}/ont/${sp}_ONT_min20000.fa" ]]; then
        ONT_FA="${ROOT}/ont/${sp}_ONT_min20000.fa"
    fi
    OUT_ROOT="${ROOT}/06_gap_telo/${SAMPLE}"

    log "============================================================"
    log "开始处理：$SAMPLE"
    log "============================================================"

    if [[ ! -s "$IN_FA" ]]; then
        log "跳过 $SAMPLE：输入 fa 尚不存在：$IN_FA（等 LR_Gapcloser 完成后再跑）"
        continue
    fi
    [[ -s "$ONT_FA" ]] || die "$SAMPLE ONT 不存在：$ONT_FA"

    AUD="${OUT_ROOT}/01_audit"
    TELO="${OUT_ROOT}/02_telo_fill"
    GAPD="${OUT_ROOT}/03_gap_fill"
    FIN="${OUT_ROOT}/04_final"
    mkdir -p "$AUD" "$TELO" "$GAPD" "$FIN"

    if [[ "$FORCE_RERUN" == "1" ]]; then
        log "FORCE_RERUN=1：清除 02/03/04 中间结果"
        rm -rf "${TELO:?}/"* "${GAPD:?}/"* "${FIN:?}/"*
    fi

    log "IN_FA=$IN_FA"
    log "ONT_FA=$ONT_FA"
    log "OUT=$OUT_ROOT"

    # ----------------------------------------------------------
    # 0. 工作拷贝
    # ----------------------------------------------------------
    WORK_FA="${OUT_ROOT}/00_input.fa"
    cp -f "$IN_FA" "$WORK_FA"
    samtools faidx "$WORK_FA"

    # ----------------------------------------------------------
    # 1. 审计：TeloExplorer + gap
    # ----------------------------------------------------------
    log "[1/4] 审计前${TOP_N} scaffold 端粒与 gap..."

    TELO_PREFIX="${AUD}/${SAMPLE}_Telo"
    run_telo_explorer "$WORK_FA" "${SAMPLE}_Telo" "$AUD"

    GAPS_TSV="${AUD}/${SAMPLE}_gaps.tsv"
    bash "$CHECK_GAP" "$WORK_FA" >"$GAPS_TSV" 2>"${AUD}/${SAMPLE}_gaps.summary.txt"

    $PYTHON3 "$HELPER" audit \
        --fasta "$WORK_FA" \
        --telo-info "${TELO_PREFIX}.telo.info" \
        --gaps-tsv "$GAPS_TSV" \
        --top-n "$TOP_N" \
        --top-list "${AUD}/top${TOP_N}.list" \
        --missing-telo "${AUD}/missing_telo.tsv" \
        --summary "${AUD}/audit_summary.txt"

    log "审计完成：${AUD}/audit_summary.txt"
    sed -n '1,5p; /top20_both/,$p' "${AUD}/audit_summary.txt"

    # ----------------------------------------------------------
    # 2. 自动补端粒（两轮：初补 + 复检后再补）
    # ----------------------------------------------------------
    log "[2/4] ONT 自动补端粒（最多 2 轮）..."

    TELO_R1_FA="${TELO}/after_telo_r1.fa"
    TELO_R1_REC="${TELO}/telo_fill_r1.tsv"
    run_telo_fill_round 1 "$WORK_FA" "${TELO_PREFIX}.telo.info" \
        "${AUD}/missing_telo.tsv" "$TELO_R1_FA" "$TELO_R1_REC" "$TELO"
    samtools faidx "$TELO_R1_FA"

    TELO_R2_PREFIX="${TELO}/${SAMPLE}_Telo_r2"
    run_telo_explorer "$TELO_R1_FA" "${SAMPLE}_Telo_r2" "$TELO"
    $PYTHON3 "$HELPER" audit \
        --fasta "$TELO_R1_FA" \
        --telo-info "${TELO_R2_PREFIX}.telo.info" \
        --gaps-tsv "$GAPS_TSV" \
        --top-n "$TOP_N" \
        --top-list "${TELO}/top${TOP_N}_r2.list" \
        --missing-telo "${TELO}/missing_telo_r2.tsv" \
        --summary "${TELO}/audit_r2_summary.txt"

    TELO_AFTER_FA="${TELO}/after_telo.fa"
    TELO_RECORD="${TELO}/telo_fill_record.tsv"
    TELO_R2_REC="${TELO}/telo_fill_r2.tsv"
    run_telo_fill_round 2 "$TELO_R1_FA" "${TELO_R2_PREFIX}.telo.info" \
        "${TELO}/missing_telo_r2.tsv" "$TELO_AFTER_FA" "$TELO_R2_REC" "$TELO"

    {
        echo -e "scaffold\tside\tstatus\tread\tadded_bp\ttelo_copies\tident\tnote"
        tail -n +2 "$TELO_R1_REC" 2>/dev/null || true
        tail -n +2 "$TELO_R2_REC" 2>/dev/null || true
    } >"$TELO_RECORD"

    samtools faidx "$TELO_AFTER_FA"

    GAPS_AFTER_TELO="${TELO}/gaps_after_telo.tsv"
    bash "$CHECK_GAP" "$TELO_AFTER_FA" >"$GAPS_AFTER_TELO" 2>"${TELO}/gaps_after_telo.summary.txt"

    # ----------------------------------------------------------
    # 3. 剩余 gap：统一 bait + 双翼 bait 两轮
    # ----------------------------------------------------------
    log "[3/4] ONT 局部比对填剩余 gap（unified + split 两轮）..."

    GAP_R1_FA="${GAPD}/after_gap_r1.fa"
    GAP_R1_REC="${GAPD}/gap_fill_r1.tsv"
    run_gap_fill_round 1 unified "$TELO_AFTER_FA" "$GAPS_AFTER_TELO" \
        "$GAP_R1_FA" "$GAP_R1_REC" "$GAPD"
    samtools faidx "$GAP_R1_FA"

    GAPS_AFTER_R1="${GAPD}/gaps_after_r1.tsv"
    bash "$CHECK_GAP" "$GAP_R1_FA" >"$GAPS_AFTER_R1" 2>"${GAPD}/gaps_after_r1.summary.txt"

    GAP_AFTER_FA="${GAPD}/after_gap.fa"
    GAP_RECORD="${GAPD}/gap_fill_record.tsv"
    GAP_R2_REC="${GAPD}/gap_fill_r2.tsv"
    run_gap_fill_round 2 split "$GAP_R1_FA" "$GAPS_AFTER_R1" \
        "$GAP_AFTER_FA" "$GAP_R2_REC" "$GAPD"

    {
        echo -e "scaffold\tgap_index\told_start\told_end\told_length\tstatus\tread\tfill_bp\tident\tnote"
        tail -n +2 "$GAP_R1_REC" 2>/dev/null || true
        tail -n +2 "$GAP_R2_REC" 2>/dev/null || true
    } >"$GAP_RECORD"

    samtools faidx "$GAP_AFTER_FA"

    # ----------------------------------------------------------
    # 4. 最终输出 + 再审计
    # ----------------------------------------------------------
    log "[4/4] 写出最终结果并复检..."

    FINAL_FA="${FIN}/${SAMPLE}_telo_gap_filled.fa"
    cp -f "$GAP_AFTER_FA" "$FINAL_FA"
    samtools faidx "$FINAL_FA"
    cp -f "$TELO_RECORD" "${FIN}/telo_fill_record.tsv"
    cp -f "$GAP_RECORD" "${FIN}/gap_fill_record.tsv"

    FINAL_TELO_PREFIX="${FIN}/${SAMPLE}_final_Telo"
    run_telo_explorer "$FINAL_FA" "${SAMPLE}_final_Telo" "$FIN"
    bash "$CHECK_GAP" "$FINAL_FA" >"${FIN}/final_gaps.tsv" 2>"${FIN}/final_gaps.summary.txt"

    $PYTHON3 "$HELPER" audit \
        --fasta "$FINAL_FA" \
        --telo-info "${FINAL_TELO_PREFIX}.telo.info" \
        --gaps-tsv "${FIN}/final_gaps.tsv" \
        --top-n "$TOP_N" \
        --top-list "${FIN}/top${TOP_N}.list" \
        --missing-telo "${FIN}/missing_telo.tsv" \
        --summary "${FIN}/final_audit_summary.txt"

    {
        echo "========== ${SAMPLE} gap/telo pipeline =========="
        echo "input:  $IN_FA"
        echo "ont:    $ONT_FA"
        echo "final:  $FINAL_FA"
        echo "params: BAIT_LEN=$BAIT_LEN GAP_FLANK=$GAP_FLANK"
        echo "        telo alen=$TELO_MIN_ALEN ident=$TELO_MIN_IDENT copies=$TELO_MIN_COPIES"
        echo "        gap  alen=$GAP_MIN_ALEN ident=$GAP_MIN_IDENT"
        echo
        echo "----- initial audit -----"
        cat "${AUD}/audit_summary.txt"
        echo
        echo "----- telo fill -----"
        if [[ -s "$TELO_RECORD" ]]; then
            awk -F'\t' 'NR==1{next} {c[$3]++} END{for (k in c) print k, c[k]}' "$TELO_RECORD"
        fi
        echo
        echo "----- gap fill -----"
        if [[ -s "$GAP_RECORD" ]]; then
            awk -F'\t' 'NR==1{next} {c[$6]++} END{for (k in c) print k, c[k]}' "$GAP_RECORD"
        fi
        echo
        echo "----- final audit -----"
        cat "${FIN}/final_audit_summary.txt"
    } | tee "${FIN}/REPORT.txt"

    log "$SAMPLE 完成。最终序列：$FINAL_FA"
    log "报告：${FIN}/REPORT.txt"
    echo "------------------------------------------------"
done

log "========== 全部样本处理结束 =========="
