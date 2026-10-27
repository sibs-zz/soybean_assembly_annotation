#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# 06_gap_telo_final_fix：智能端粒/gap 补全
#
# 分流：
#   tip_has_telo      → skip
#   buried_near_tip   → ONT 证据分流（默认非纯 motif trim）
#       无支持 → trim（AUTO_TRIM=1）
#       有远端支持 → trial_trim 一次：新 tip motif+quarTeT 达标且 ≤MAX_TRIM_BP 则保留，
#                    否则回滚并 unresolved；裁掉远端写入 clipped_*
#       仅 tip overhang → extend 只试 1 轮
#       证据不足 → unresolved
#   true_missing      → ONT tip 外延
#   deep_internal     → unresolved（不改序列）
#
# FORCE_MOTIF_TRIM=1：跳过 ONT，纯 motif trim（旧行为）
# AUTO_TRIM=0：即使无 ONT 支持也不 trim，改标 unresolved
#
# perfect：仅 quarTeT（-m 50）双端齐全（40 端）+ 0 gap + 无 unresolved
# 达标才替换正式 04_final（替换前备份）；开跑前 -m50 预检通过则保留原输出
# ============================================================

# ======================== 配置区 =============================

SAMPLES=(
#    "BZ_primary"
    "DY11_primary"
    "TH2_primary"
#    "TZX1998_primary"
#    "TZX253_primary"
)

if [[ -n "${SAMPLE:-}" ]]; then
    SAMPLES=("$SAMPLE")
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${ROOT:-/tiandata2/zzh/sjy_assembly/batch2}"
TOP_N=20
THREADS="${THREADS:-64}"
FORCE_RERUN="${FORCE_RERUN:-0}"
FORCE_REPLACE="${FORCE_REPLACE:-0}"
FORCE_IMPERFECT="${FORCE_IMPERFECT:-0}"
APPLY_REVCOMP="${APPLY_REVCOMP:-1}"
FLIP_ONESIDED="${FLIP_ONESIDED:-0}"
FORCE_MOTIF_TRIM="${FORCE_MOTIF_TRIM:-0}"        # 1=跳过 ONT，纯 motif trim
AUTO_TRIM="${AUTO_TRIM:-1}"                      # 1=无 ONT 支持则 trim；0=改标 unresolved
ALLOW_DISTAL_TRIAL_TRIM="${ALLOW_DISTAL_TRIAL_TRIM:-1}"  # distal_supported 试 trim+验收

BAIT_LEN=150000
GAP_FLANK=100000
USE_FULL_ONT=1
TELO_ROUNDS=4
GAP_ROUNDS=4

# 判定标准：仅 quarTeT；motif 参数只用于 buried 结构分流，不参与 perfect
QUARTET_MIN_REPEAT=50   # TeloExplorer -m（默认工具为 100，此处放宽）
MOTIF_TIP_BP=20000
MOTIF_MIN_COPIES=50
INTERNAL_MIN_COPIES=50
NEAR_TIP_BP=500000
MAX_TRIM_BP=500000
BURY_ANCHOR_BP=50000
BURY_MIN_SUPPORT_READS=3
BURY_MIN_EXTEND_READS=1
BURY_MIN_ALEN=500
BURY_MIN_IDENT=0.55

# ONT tip fill（仅 true_missing）
TELO_MIN_ALEN=300
TELO_MIN_IDENT=0.55
TELO_MIN_COPIES=5
TELO_EDGE_BP=15000
TELO_MIN_OVERHANG=30

GAP_MIN_ALEN=300
GAP_MIN_IDENT=0.50
GAP_MAX_FILL=15000

QUARTET="/tiandata2/zzh/tools/quarTeT/quartet.py"
CHECK_GAP="${SCRIPT_DIR}/check_gap.sh"
HELPER="${SCRIPT_DIR}/06_gap_telo_lib.py"
FIXLIB="${SCRIPT_DIR}/06_gap_telo_final_fix_lib.py"
SEQKIT="${SEQKIT:-/home/zzh/anaconda3/envs/py37/bin/seqkit}"
MINIMAP2="${MINIMAP2:-minimap2}"
CONDA_ENV="quarTeTdependencies"

# ============================================================

timestamp() { date "+%F %T"; }
log() { echo "[$(timestamp)] $*"; }
die() { echo "[$(timestamp)] ERROR: $*" >&2; exit 1; }

run_telo_explorer() {
    local fa="$1" prefix="$2" outdir="$3"
    (
        cd "$outdir"
        $PYTHON3 "$QUARTET" TeloExplorer -i "$fa" -c plant -m "$QUARTET_MIN_REPEAT" -p "$prefix"
    )
}

# quarTeT 审计（主路径；-m=$QUARTET_MIN_REPEAT）
# 注意：仅最后一行 echo 路径串到 stdout，供 IFS=| read；日志一律 stderr，避免污染解析
refresh_quartet() {
    local fa="$1" workdir="$2" tag="$3"
    local prefix="${SAMPLE}_${tag}_Telo"
    local gaps="${workdir}/${tag}_gaps.tsv"
    local info="${workdir}/${prefix}.telo.info"
    local missing="${workdir}/missing_telo_${tag}.tsv"
    local toplist="${workdir}/top${TOP_N}_${tag}.list"
    local summary="${workdir}/audit_quartet_${tag}.txt"
    bash "$CHECK_GAP" "$fa" >"$gaps" 2>"${workdir}/${tag}_gaps.summary.txt"
    run_telo_explorer "$fa" "$prefix" "$workdir" >&2
    [[ -s "$info" ]] || die "quarTeT 未产出：$info"
    $PYTHON3 "$HELPER" audit \
        --fasta "$fa" \
        --telo-info "$info" \
        --gaps-tsv "$gaps" \
        --top-n "$TOP_N" \
        --top-list "$toplist" \
        --missing-telo "$missing" \
        --summary "$summary" >&2
    # echo: telo_info|gaps|missing|top_list|summary
    echo "${info}|${gaps}|${missing}|${toplist}|${summary}"
}

run_telo_fill_round() {
    local round="$1" in_fa="$2" telo_info="$3" missing_tsv="$4" top_list="$5" out_fa="$6" record="$7" workdir="$8"
    local bait_fa="${workdir}/telo_baits_r${round}.fa"
    local paf="${workdir}/telo_baits_r${round}.paf"
    local hit_ids="${workdir}/telo_hit_ids_r${round}.txt"
    local hit_fa="${workdir}/telo_hit_reads_r${round}.fa"
    local min_alen="$TELO_MIN_ALEN"
    local min_ident="$TELO_MIN_IDENT"
    local min_copies="$TELO_MIN_COPIES"

    if [[ "$round" -ge 3 ]]; then
        min_alen=200
        min_ident=0.50
        min_copies=3
    fi

    if [[ ! -s "$missing_tsv" ]]; then
        log "端粒 round${round}：无 true_missing，跳过"
        cp -f "$in_fa" "$out_fa"
        printf "scaffold\tside\tstatus\tread\tadded_bp\ttelo_copies\tident\tnote\n" >"$record"
        return 0
    fi

    log "端粒 round${round}：$(wc -l <"$missing_tsv") 个 true_missing (alen>=${min_alen} ident>=${min_ident} copies>=${min_copies})"
    $PYTHON3 "$HELPER" make-telo-baits \
        --fasta "$in_fa" \
        --missing-telo "$missing_tsv" \
        --bait-len "$BAIT_LEN" \
        --out-fa "$bait_fa"

    log "minimap2 map-ont：端粒 bait round${round} <- ONT..."
    "$MINIMAP2" -x map-ont -t "$THREADS" "$bait_fa" "$ONT_FA" >"$paf"

    awk '{print $1}' "$paf" | sort -u >"$hit_ids"
    if [[ ! -s "$hit_ids" ]]; then
        log "端粒 round${round}：无 ONT hit"
        cp -f "$in_fa" "$out_fa"
        printf "scaffold\tside\tstatus\tread\tadded_bp\ttelo_copies\tident\tnote\n" >"$record"
        return 0
    fi
    "$SEQKIT" grep -f "$hit_ids" "$ONT_FA" >"$hit_fa"

    $PYTHON3 "$HELPER" apply-telo \
        --fasta "$in_fa" \
        --telo-info "$telo_info" \
        --top-list "$top_list" \
        --missing-telo "$missing_tsv" \
        --paf "$paf" \
        --reads-fa "$hit_fa" \
        --out-fa "$out_fa" \
        --record "$record" \
        --min-alen "$min_alen" \
        --min-ident "$min_ident" \
        --min-telo-copies "$min_copies" \
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
    local max_fill="$GAP_MAX_FILL"

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
    if [[ "$round" -ge 3 ]]; then
        min_alen=200
        min_ident=0.45
        max_fill=25000
    elif [[ "$round" -ge 2 ]]; then
        min_alen=200
        min_ident=0.48
        max_fill=20000
    fi

    log "gap round${round} (${mode})：${gap_n} 个 gap (alen>=${min_alen} ident>=${min_ident} max_fill=${max_fill})"
    $PYTHON3 "$HELPER" make-gap-baits \
        --fasta "$in_fa" \
        --gaps-tsv "$gaps_tsv" \
        --flank "$GAP_FLANK" \
        --out-fa "$bait_fa" \
        "${unified_flag[@]}"

    log "minimap2 map-ont：gap bait round${round} <- ONT..."
    "$MINIMAP2" -x map-ont -t "$THREADS" "$bait_fa" "$ONT_FA" >"$paf"

    awk '{print $1}' "$paf" | sort -u >"$hit_ids"
    if [[ ! -s "$hit_ids" ]]; then
        log "gap round${round}：无 ONT hit"
        cp -f "$in_fa" "$out_fa"
        printf "scaffold\tgap_index\told_start\told_end\told_length\tstatus\tread\tfill_bp\tident\tnote\n" >"$record"
        return 0
    fi
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
        --max-fill "$max_fill" \
        "${unified_flag[@]}"
}

# -------------------- conda --------------------
if [[ -z "${CONDA_DEFAULT_ENV:-}" || "${CONDA_DEFAULT_ENV}" != "$CONDA_ENV" ]]; then
    command -v conda >/dev/null 2>&1 || die "未找到 conda"
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$CONDA_ENV"
fi
log "conda env: ${CONDA_DEFAULT_ENV} ($(command -v python3))"

PYTHON3="${CONDA_PREFIX}/bin/python3"
[[ -x "$PYTHON3" ]] || die "环境 ${CONDA_ENV} 中找不到 python3：$PYTHON3"
[[ -f "$HELPER" ]] || die "缺少 helper：$HELPER"
[[ -f "$FIXLIB" ]] || die "缺少 fix lib：$FIXLIB"
[[ -f "$CHECK_GAP" ]] || die "缺少 check_gap.sh"
[[ -x "$SEQKIT" || -f "$SEQKIT" ]] || die "缺少 seqkit：$SEQKIT"
command -v "$MINIMAP2" >/dev/null 2>&1 || die "缺少 minimap2"
[[ -f "$QUARTET" ]] || die "缺少 quarTeT：$QUARTET"

log "perfect：quarTeT -m=$QUARTET_MIN_REPEAT only | buried 分流仍用 tip motif 结构"
log "buried：FORCE_MOTIF_TRIM=$FORCE_MOTIF_TRIM AUTO_TRIM=$AUTO_TRIM DISTAL_TRIAL_TRIM=$ALLOW_DISTAL_TRIAL_TRIM support_reads>=$BURY_MIN_SUPPORT_READS"
log "ONT tip：alen>=$TELO_MIN_ALEN ident>=$TELO_MIN_IDENT rounds=$TELO_ROUNDS"
log "gap：alen>=$GAP_MIN_ALEN ident>=$GAP_MIN_IDENT max_fill=$GAP_MAX_FILL rounds=$GAP_ROUNDS"
log "FORCE_RERUN=$FORCE_RERUN FORCE_REPLACE=$FORCE_REPLACE FORCE_IMPERFECT=$FORCE_IMPERFECT"

# 预检通过：追加写入正式 04_final/REPORT.txt（不改序列）
append_keep_report() {
    local report="$1" summary="$2" check_msg="$3" preflight_dir="$4"
    mkdir -p "$(dirname "$report")"
    {
        echo
        echo "========== final_fix preflight $(timestamp) =========="
        echo "quarTeT TeloExplorer -m ${QUARTET_MIN_REPEAT}: PASS"
        echo "action: keep original 04_final sequences (no rewrite)"
        echo "check: ${check_msg}"
        echo "preflight_dir: ${preflight_dir}"
        echo "----- quarTeT -m ${QUARTET_MIN_REPEAT} audit -----"
        if [[ -s "$summary" ]]; then
            cat "$summary"
        else
            echo "(no audit summary)"
        fi
    } >>"$report"
    log "已追加预检通过记录 -> $report"
}

# -------------------- quarTeT -m50 预检（全部样本） --------------------
log "========== quarTeT -m ${QUARTET_MIN_REPEAT} 预检 =========="
NEED_SAMPLES=()
SKIP_KEEP_SAMPLES=()
for SAMPLE in "${SAMPLES[@]}"; do
    sp="${SAMPLE%%_*}"
    FIN06="${ROOT}/06_gap_telo/${SAMPLE}/04_final"
    IN_FINAL_FA="${FIN06}/${SAMPLE}_telo_gap_filled.fa"
    IN_GAPCLOSED="${ROOT}/05_LRgapcloser_out/${SAMPLE}/${SAMPLE}_gapclosed.fa"
    ONT_FA="${ROOT}/ont/${sp}_ONT_clean.fa"
    if [[ "$USE_FULL_ONT" != "1" && -s "${ROOT}/ont/${sp}_ONT_min20000.fa" ]]; then
        ONT_FA="${ROOT}/ont/${sp}_ONT_min20000.fa"
    fi

    if [[ ! -s "$ONT_FA" ]]; then
        log "预检跳过 $SAMPLE：无 ONT"
        continue
    fi
    if [[ ! -s "$IN_FINAL_FA" && ! -s "$IN_GAPCLOSED" ]]; then
        log "预检跳过 $SAMPLE：无输入 fa"
        continue
    fi

    # 尚无 04_final：必须进入下游
    if [[ ! -s "$IN_FINAL_FA" ]]; then
        log "预检 $SAMPLE：无 04_final，列入处理"
        NEED_SAMPLES+=("$SAMPLE")
        continue
    fi

    PRE="${ROOT}/06_gap_telo_final_fix/${SAMPLE}/00_preflight"
    mkdir -p "$PRE"
    log "预检 $SAMPLE：quarTeT TeloExplorer -m ${QUARTET_MIN_REPEAT} ..."
    PRE_GAPS="${PRE}/preflight_gaps.tsv"
    bash "$CHECK_GAP" "$IN_FINAL_FA" >"$PRE_GAPS" 2>"${PRE}/preflight_gaps.summary.txt" || true
    run_telo_explorer "$IN_FINAL_FA" "${SAMPLE}_preflight_Telo" "$PRE" || {
        log "预检 $SAMPLE：quarTeT 失败，列入处理"
        NEED_SAMPLES+=("$SAMPLE")
        continue
    }
    PRE_INFO="${PRE}/${SAMPLE}_preflight_Telo.telo.info"
    PRE_TOP="${PRE}/top${TOP_N}_preflight.list"
    PRE_MISS="${PRE}/missing_telo_preflight.tsv"
    PRE_SUM="${PRE}/audit_quartet_preflight.txt"
    if [[ ! -s "$PRE_INFO" ]]; then
        log "预检 $SAMPLE：无 .telo.info，列入处理"
        NEED_SAMPLES+=("$SAMPLE")
        continue
    fi
    $PYTHON3 "$HELPER" audit \
        --fasta "$IN_FINAL_FA" \
        --telo-info "$PRE_INFO" \
        --gaps-tsv "$PRE_GAPS" \
        --top-n "$TOP_N" \
        --top-list "$PRE_TOP" \
        --missing-telo "$PRE_MISS" \
        --summary "$PRE_SUM" || true
    cat "$PRE_SUM" || true

    PRE_CHECK_MSG=""
    if PRE_CHECK_MSG=$($PYTHON3 "$FIXLIB" check-perfect \
        --mode quartet \
        --fasta "$IN_FINAL_FA" \
        --gaps-tsv "$PRE_GAPS" \
        --top-list "$PRE_TOP" \
        --top-n "$TOP_N" \
        --telo-info "$PRE_INFO" 2>&1)
    then
        echo "$PRE_CHECK_MSG"
        if [[ "$FORCE_IMPERFECT" == "1" ]]; then
            log "预检 $SAMPLE：quarTeT -m${QUARTET_MIN_REPEAT} 已满分，但 FORCE_IMPERFECT=1，仍处理"
            NEED_SAMPLES+=("$SAMPLE")
        else
            log "预检 $SAMPLE：quarTeT -m${QUARTET_MIN_REPEAT} 双端齐全且 0 gap → 可保留原 04_final，不处理"
            append_keep_report \
                "${FIN06}/REPORT.txt" \
                "$PRE_SUM" \
                "$PRE_CHECK_MSG" \
                "$PRE"
            # 同步一份 -m50 预检产物到 04_final，便于复核
            cp -f "$PRE_INFO" "${FIN06}/${SAMPLE}_preflight_m${QUARTET_MIN_REPEAT}.telo.info"
            cp -f "$PRE_SUM" "${FIN06}/preflight_m${QUARTET_MIN_REPEAT}_audit.txt"
            SKIP_KEEP_SAMPLES+=("$SAMPLE")
        fi
    else
        echo "$PRE_CHECK_MSG"
        log "预检 $SAMPLE：quarTeT -m${QUARTET_MIN_REPEAT} 未过 → 进入下游填补"
        NEED_SAMPLES+=("$SAMPLE")
    fi
done

log "========== 预检汇总 =========="
if [[ ${#SKIP_KEEP_SAMPLES[@]} -gt 0 ]]; then
    log "保留原 04_final（不改）：${SKIP_KEEP_SAMPLES[*]}"
else
    log "保留原 04_final：无"
fi
if [[ ${#NEED_SAMPLES[@]} -gt 0 ]]; then
    log "将下游填补：${NEED_SAMPLES[*]}"
else
    log "将下游填补：无（全部预检通过或无可处理样本）"
    log "========== final_fix 全部结束 =========="
    exit 0
fi
SAMPLES=("${NEED_SAMPLES[@]}")

for SAMPLE in "${SAMPLES[@]}"; do
    sp="${SAMPLE%%_*}"
    BASE06="${ROOT}/06_gap_telo/${SAMPLE}"
    FIN06="${BASE06}/04_final"
    IN_FINAL_FA="${FIN06}/${SAMPLE}_telo_gap_filled.fa"
    IN_GAPCLOSED="${ROOT}/05_LRgapcloser_out/${SAMPLE}/${SAMPLE}_gapclosed.fa"
    ONT_FA="${ROOT}/ont/${sp}_ONT_clean.fa"
    if [[ "$USE_FULL_ONT" != "1" && -s "${ROOT}/ont/${sp}_ONT_min20000.fa" ]]; then
        ONT_FA="${ROOT}/ont/${sp}_ONT_min20000.fa"
    fi

    OUT_ROOT="${ROOT}/06_gap_telo_final_fix/${SAMPLE}"
    AUD="${OUT_ROOT}/01_diagnose"
    TELO="${OUT_ROOT}/02_telo_fill"
    GAPD="${OUT_ROOT}/03_gap_fill"
    FIN="${OUT_ROOT}/04_final"

    log "============================================================"
    log "final_fix：$SAMPLE"
    log "============================================================"

    if [[ ! -s "$ONT_FA" ]]; then
        log "跳过 $SAMPLE：ONT 不存在：$ONT_FA"
        continue
    fi

    if [[ -s "$IN_FINAL_FA" ]]; then
        IN_FA="$IN_FINAL_FA"
        log "输入：现有 04_final"
    elif [[ -s "$IN_GAPCLOSED" ]]; then
        IN_FA="$IN_GAPCLOSED"
        log "输入：05 gapclosed（尚无 04_final）"
    else
        log "跳过 $SAMPLE：无可用输入 fa"
        continue
    fi

    mkdir -p "$AUD" "$TELO" "$GAPD" "$FIN" "$OUT_ROOT"

    # 保留进入 final_fix 前的 06 fill 记录，结束时与本次 fill 合并（供 07 polish）
    BASE_TELO="${OUT_ROOT}/00_base_telo_fill_record.tsv"
    BASE_GAP="${OUT_ROOT}/00_base_gap_fill_record.tsv"
    TELO_HDR="scaffold\tside\tstatus\tread\tadded_bp\ttelo_copies\tident\tnote"
    GAP_HDR="scaffold\tgap_index\told_start\told_end\told_length\tstatus\tread\tfill_bp\tident\tnote"
    if [[ -s "${FIN06}/telo_fill_record.tsv" ]]; then
        cp -f "${FIN06}/telo_fill_record.tsv" "$BASE_TELO"
    else
        echo -e "$TELO_HDR" >"$BASE_TELO"
    fi
    if [[ -s "${FIN06}/gap_fill_record.tsv" ]]; then
        cp -f "${FIN06}/gap_fill_record.tsv" "$BASE_GAP"
    else
        echo -e "$GAP_HDR" >"$BASE_GAP"
    fi
    log "base fill records: telo=$(grep -cve '^$' "$BASE_TELO" || true) lines gap=$(grep -cve '^$' "$BASE_GAP" || true) lines"

    if [[ "$FORCE_RERUN" == "1" ]]; then
        log "FORCE_RERUN=1：清除 final_fix 中间结果"
        rm -rf "${AUD:?}/"* "${TELO:?}/"* "${GAPD:?}/"* "${FIN:?}/"*
    fi

    WORK_FA="${OUT_ROOT}/00_input.fa"
    cp -f "$IN_FA" "$WORK_FA"
    samtools faidx "$WORK_FA"

    # ----------------------------------------------------------
    # 1. quarTeT 审计 + 诊断分流（缺失列表以 quarTeT 为准）
    # ----------------------------------------------------------
    log "[1/6] quarTeT 审计（-m ${QUARTET_MIN_REPEAT}）+ 诊断分流..."
    IFS='|' read -r TELO_INFO GAPS_TSV MISSING_TSV TOP_LIST AUDIT_TXT < <(
        refresh_quartet "$WORK_FA" "$AUD" "init"
    )
    if [[ -n "${AUDIT_TXT:-}" && -s "$AUDIT_TXT" ]]; then
        cat "$AUDIT_TXT"
    else
        log "WARN: quarTeT audit summary 未解析到（AUDIT_TXT='${AUDIT_TXT:-}'）"
    fi

    DIAG_TSV="${AUD}/orient_diagnose.tsv"
    REV_LIST="${AUD}/revcomp_targets.txt"
    TRIM_CAND="${AUD}/trim_candidates.tsv"
    TRIM_PLAN="${AUD}/trim_plan.tsv"
    FILL_MISS="${AUD}/fill_missing.tsv"
    UNRESOLVED="${AUD}/unresolved.tsv"
    BURY_DECIDE="${AUD}/buried_decide.tsv"
    DIAG_EXTRA=(--trust-missing)
    [[ "$FLIP_ONESIDED" == "1" ]] && DIAG_EXTRA+=(--flip-onesided)

    $PYTHON3 "$FIXLIB" diagnose-orient \
        --fasta "$WORK_FA" \
        --missing-telo "$MISSING_TSV" \
        --telo-info "$TELO_INFO" \
        --top-list "$TOP_LIST" \
        --top-n "$TOP_N" \
        --out-tsv "$DIAG_TSV" \
        --revcomp-list "$REV_LIST" \
        --trim-plan "$TRIM_CAND" \
        --fill-missing "$FILL_MISS" \
        --tip-bp "$MOTIF_TIP_BP" \
        --tip-min-copies "$MOTIF_MIN_COPIES" \
        --internal-min-copies "$INTERNAL_MIN_COPIES" \
        --near-tip-bp "$NEAR_TIP_BP" \
        --max-trim-bp "$MAX_TRIM_BP" \
        "${DIAG_EXTRA[@]}"

    if [[ "$APPLY_REVCOMP" == "1" && -s "$REV_LIST" ]]; then
        log "revcomp $(wc -l <"$REV_LIST") scaffolds..."
        $PYTHON3 "$FIXLIB" apply-revcomp \
            --fasta "$WORK_FA" \
            --revcomp-list "$REV_LIST" \
            --out-fa "${AUD}/after_revcomp.fa"
        cp -f "${AUD}/after_revcomp.fa" "$WORK_FA"
        samtools faidx "$WORK_FA"
        IFS='|' read -r TELO_INFO GAPS_TSV MISSING_TSV TOP_LIST AUDIT_TXT < <(
            refresh_quartet "$WORK_FA" "$AUD" "post_revcomp"
        )
        $PYTHON3 "$FIXLIB" diagnose-orient \
            --fasta "$WORK_FA" \
            --missing-telo "$MISSING_TSV" \
            --telo-info "$TELO_INFO" \
            --top-list "$TOP_LIST" \
            --out-tsv "$DIAG_TSV" \
            --revcomp-list "$REV_LIST" \
            --trim-plan "$TRIM_CAND" \
            --fill-missing "$FILL_MISS" \
            --tip-bp "$MOTIF_TIP_BP" \
            --tip-min-copies "$MOTIF_MIN_COPIES" \
            --internal-min-copies "$INTERNAL_MIN_COPIES" \
            --near-tip-bp "$NEAR_TIP_BP" \
            --max-trim-bp "$MAX_TRIM_BP" \
            --trust-missing
    fi

    # ----------------------------------------------------------
    # 2. buried：ONT 证据三分流（默认非纯 motif trim）
    # ----------------------------------------------------------
    log "[2/6] buried_near_tip：ONT 证据分流..."
    : >"$UNRESOLVED"
    echo -e "scaffold\tside\tstatus\tspan_reads\tdistal_reads\textend_reads\ttrim_bp\tnote" >"$UNRESOLVED"
    printf "scaffold\tside\ttrim_bp\tkeep_start\tkeep_end\tcluster_start\tcluster_end\tcluster_copies\tnew_tip_copies\tnote\n" >"$TRIM_PLAN"
    BURY_EXTEND="${AUD}/bury_extend.tsv"
    : >"$BURY_EXTEND"
    CLIP_FA="${AUD}/clipped_distal.fa"
    CLIP_TSV="${AUD}/clipped_distal.tsv"
    : >"$CLIP_FA"
    printf "scaffold\tside\ttrim_bp\tclipped_start\tclipped_end\tclipped_len\tnote\n" >"$CLIP_TSV"
    TRIAL_ACCEPT="${AUD}/trial_trim_accept.tsv"
    printf "scaffold\tside\tstatus\ttip_copies\tquartet_ok\ttrim_bp\tnote\n" >"$TRIAL_ACCEPT"

    cand_n=0
    if [[ -s "$TRIM_CAND" ]]; then
        cand_n=$(($(wc -l <"$TRIM_CAND") - 1))
    fi

    if [[ "$cand_n" -le 0 ]]; then
        log "无 buried 候选，跳过 ONT bury 检查"
        cp -f "$WORK_FA" "${AUD}/after_trim.fa"
        printf "scaffold\tside\tstatus\ttrim_bp\told_len\tnew_len\tnote\n" >"${AUD}/trim_record.tsv"
    elif [[ "$FORCE_MOTIF_TRIM" == "1" ]]; then
        log "FORCE_MOTIF_TRIM=1：跳过 ONT，按 motif 候选直接 trim"
        cp -f "$TRIM_CAND" "$TRIM_PLAN"
        $PYTHON3 "$FIXLIB" trim-buried \
            --fasta "$WORK_FA" \
            --trim-plan "$TRIM_PLAN" \
            --out-fa "${AUD}/after_trim.fa" \
            --record "${AUD}/trim_record.tsv" \
            --max-trim-bp "$MAX_TRIM_BP"
        cp -f "${AUD}/after_trim.fa" "$WORK_FA"
        samtools faidx "$WORK_FA"
    else
        BAIT_FA="${AUD}/buried_baits.fa"
        BAIT_COORDS="${AUD}/buried_baits.coords.tsv"
        BAIT_PAF="${AUD}/buried_baits.paf"
        BAIT_IDS="${AUD}/buried_hit_ids.txt"
        BAIT_READS="${AUD}/buried_hit_reads.fa"

        $PYTHON3 "$FIXLIB" make-buried-baits \
            --fasta "$WORK_FA" \
            --candidates "$TRIM_CAND" \
            --out-fa "$BAIT_FA" \
            --out-coords "$BAIT_COORDS" \
            --anchor-bp "$BURY_ANCHOR_BP"

        if [[ ! -s "$BAIT_FA" ]]; then
            log "buried bait 为空，全部标 unresolved"
            awk -F'\t' 'NR==1{next} {print $1"\t"$2"\tunresolved\t0\t0\t0\t"$3"\tno_bait"}' \
                "$TRIM_CAND" >>"$UNRESOLVED"
            cp -f "$WORK_FA" "${AUD}/after_trim.fa"
            printf "scaffold\tside\tstatus\ttrim_bp\told_len\tnew_len\tnote\n" >"${AUD}/trim_record.tsv"
        else
            log "minimap2 map-ont：buried bait <- ONT..."
            "$MINIMAP2" -x map-ont -t "$THREADS" "$BAIT_FA" "$ONT_FA" >"$BAIT_PAF"
            awk '{print $1}' "$BAIT_PAF" | sort -u >"$BAIT_IDS"
            if [[ -s "$BAIT_IDS" ]]; then
                "$SEQKIT" grep -f "$BAIT_IDS" "$ONT_FA" >"$BAIT_READS"
            else
                : >"$BAIT_READS"
            fi

            $PYTHON3 "$FIXLIB" decide-buried-ont \
                --candidates "$TRIM_CAND" \
                --coords "$BAIT_COORDS" \
                --paf "$BAIT_PAF" \
                --reads-fa "$BAIT_READS" \
                --out-trim-plan "$TRIM_PLAN" \
                --out-extend "$BURY_EXTEND" \
                --out-unresolved "$UNRESOLVED" \
                --out-decide "$BURY_DECIDE" \
                --min-alen "$BURY_MIN_ALEN" \
                --min-ident "$BURY_MIN_IDENT" \
                --min-support-reads "$BURY_MIN_SUPPORT_READS" \
                --min-extend-reads "$BURY_MIN_EXTEND_READS" \
                --min-overhang "$TELO_MIN_OVERHANG" \
                --min-telo-copies "$TELO_MIN_COPIES" \
                --auto-trim "$AUTO_TRIM" \
                --allow-distal-trial-trim "$ALLOW_DISTAL_TRIAL_TRIM" \
                --min-new-tip-copies "$MOTIF_MIN_COPIES" \
                --max-trim-bp "$MAX_TRIM_BP"

            [[ -s "$BURY_DECIDE" ]] && cat "$BURY_DECIDE" || true

            # 试 trim 前备份，验收失败可回滚
            cp -f "$WORK_FA" "${AUD}/before_trim.fa"
            CLIP_FA="${AUD}/clipped_distal.fa"
            CLIP_TSV="${AUD}/clipped_distal.tsv"
            $PYTHON3 "$FIXLIB" trim-buried \
                --fasta "$WORK_FA" \
                --trim-plan "$TRIM_PLAN" \
                --out-fa "${AUD}/after_trim.fa" \
                --record "${AUD}/trim_record.tsv" \
                --max-trim-bp "$MAX_TRIM_BP" \
                --out-clipped-fa "$CLIP_FA" \
                --out-clipped-tsv "$CLIP_TSV"

            # distal trial_trim：quarTeT + motif 验收；不达标回滚该 scaffold
            TRIAL_N=0
            if [[ -s "$BURY_DECIDE" ]]; then
                TRIAL_N=$(awk -F'\t' 'NR>1 && $3=="trial_trim"{c++} END{print c+0}' "$BURY_DECIDE")
            fi
            if [[ "$TRIAL_N" -gt 0 ]]; then
                log "trial_trim ${TRIAL_N} 端：quarTeT -m ${QUARTET_MIN_REPEAT} 验收..."
                IFS='|' read -r _TELO_INFO _GAPS _MISS _TOP _AUDIT < <(
                    refresh_quartet "${AUD}/after_trim.fa" "$AUD" "trial_trim"
                )
                $PYTHON3 "$FIXLIB" accept-trial-trim \
                    --pre-fa "${AUD}/before_trim.fa" \
                    --post-fa "${AUD}/after_trim.fa" \
                    --decide "$BURY_DECIDE" \
                    --telo-info "$_TELO_INFO" \
                    --out-fa "${AUD}/after_trim.fa" \
                    --out-accept "$TRIAL_ACCEPT" \
                    --out-unresolved "$UNRESOLVED" \
                    --tip-bp "$MOTIF_TIP_BP" \
                    --min-copies "$MOTIF_MIN_COPIES" \
                    --clipped-tsv "$CLIP_TSV" \
                    --clipped-fa "$CLIP_FA" \
                    --out-clipped-tsv "$CLIP_TSV" \
                    --out-clipped-fa "$CLIP_FA"
                [[ -s "$TRIAL_ACCEPT" ]] && cat "$TRIAL_ACCEPT" || true
            fi

            cp -f "${AUD}/after_trim.fa" "$WORK_FA"
            samtools faidx "$WORK_FA"

            # tip-overhang-only extend：只并入一次，供第 1 轮 fill
            if [[ -s "$BURY_EXTEND" ]]; then
                cat "$BURY_EXTEND" >>"$FILL_MISS"
                sort -u "$FILL_MISS" -o "$FILL_MISS"
                log "buried tip-overhang extend（仅 1 轮）：$(wc -l <"$BURY_EXTEND") 端"
            fi
        fi
    fi

    TRIM_REC="${AUD}/trim_record.tsv"
    if [[ -s "$TRIM_REC" ]]; then
        awk -F'\t' 'NR==1 || $3=="OK"' "$TRIM_REC" || true
    fi
    if [[ -s "$UNRESOLVED" ]] && [[ $(($(wc -l <"$UNRESOLVED") - 1)) -gt 0 ]]; then
        log "unresolved 端：$(($(wc -l <"$UNRESOLVED") - 1))"
        cat "$UNRESOLVED"
    fi
    # deep_internal / diagnose 直接标 unresolved 的端并入清单
    if [[ -s "$DIAG_TSV" ]]; then
        awk -F'\t' 'NR>1 && ($7=="deep_internal" || $8=="unresolved") {
            print $1"\t"$2"\tunresolved\t0\t0\t0\t"$9"\t"$7":"$10
        }' "$DIAG_TSV" >>"$UNRESOLVED" || true
    fi

    IFS='|' read -r TELO_INFO GAPS_TSV MISSING_TSV TOP_LIST AUDIT_TXT < <(
        refresh_quartet "$WORK_FA" "$AUD" "post_trim"
    )
    # 刷新 true_missing；bury_extend 只在首轮 fill 用一次
    FILL_DIAG="${AUD}/fill_missing_post_trim.tsv"
    $PYTHON3 "$FIXLIB" diagnose-orient \
        --fasta "$WORK_FA" \
        --missing-telo "$MISSING_TSV" \
        --telo-info "$TELO_INFO" \
        --top-list "$TOP_LIST" \
        --out-tsv "${AUD}/orient_diagnose_post_trim.tsv" \
        --revcomp-list "${AUD}/revcomp_post_trim.txt" \
        --trim-plan "${AUD}/trim_candidates_post_trim.tsv" \
        --fill-missing "$FILL_DIAG" \
        --tip-bp "$MOTIF_TIP_BP" \
        --tip-min-copies "$MOTIF_MIN_COPIES" \
        --internal-min-copies "$INTERNAL_MIN_COPIES" \
        --near-tip-bp "$NEAR_TIP_BP" \
        --max-trim-bp "$MAX_TRIM_BP" \
        --trust-missing
    {
        cat "$FILL_DIAG" 2>/dev/null || true
        if [[ -s "$BURY_EXTEND" ]]; then cat "$BURY_EXTEND"; fi
    } | awk 'NF' | sort -u >"$FILL_MISS" || true
    cat "${AUD}/audit_quartet_post_trim.txt" || true

    # ----------------------------------------------------------
    # 3. ONT 补端：true_missing 多轮；bury_extend 仅第 1 轮
    #    quarTeT 20/20 或 missing 无进展则停
    # ----------------------------------------------------------
    log "[3/6] ONT 补端粒（最多 ${TELO_ROUNDS} 轮；extend 仅 r1；停轮看 quarTeT）..."
    CUR_FA="$WORK_FA"
    echo -e "scaffold\tside\tstatus\tread\tadded_bp\ttelo_copies\tident\tnote" >"${TELO}/telo_fill_record.tsv"
    PREV_MISS=-1

    for ((r=1; r<=TELO_ROUNDS; r++)); do
        tag="telo_r${r}"
        IFS='|' read -r TELO_INFO GAPS_TSV MISSING_TSV TOP_LIST AUDIT_TXT < <(
            refresh_quartet "$CUR_FA" "$TELO" "$tag"
        )
        # quarTeT 已满分：提前结束
        if [[ -s "$AUDIT_TXT" ]] && grep -qE $'top[0-9]*_both_telomere\t'"${TOP_N}/${TOP_N}"$'\tYES' "$AUDIT_TXT"; then
            log "quarTeT 已 ${TOP_N}/${TOP_N}（-m ${QUARTET_MIN_REPEAT}），跳过后续端粒轮次"
            cp -f "$CUR_FA" "${TELO}/after_telo_r${r}.fa"
            break
        fi

        $PYTHON3 "$FIXLIB" diagnose-orient \
            --fasta "$CUR_FA" \
            --missing-telo "$MISSING_TSV" \
            --telo-info "$TELO_INFO" \
            --top-list "$TOP_LIST" \
            --out-tsv "${TELO}/orient_${tag}.tsv" \
            --revcomp-list "${TELO}/revcomp_${tag}.txt" \
            --trim-plan "${TELO}/trim_${tag}.tsv" \
            --fill-missing "${TELO}/fill_${tag}.diag.tsv" \
            --tip-bp "$MOTIF_TIP_BP" \
            --tip-min-copies "$MOTIF_MIN_COPIES" \
            --internal-min-copies "$INTERNAL_MIN_COPIES" \
            --near-tip-bp "$NEAR_TIP_BP" \
            --max-trim-bp "$MAX_TRIM_BP" \
            --trust-missing

        # r1：true_missing ∪ bury_extend；之后只用 true_missing
        if [[ "$r" -eq 1 ]]; then
            {
                cat "${TELO}/fill_${tag}.diag.tsv" 2>/dev/null || true
                if [[ -s "$BURY_EXTEND" ]]; then cat "$BURY_EXTEND"; fi
            } | awk 'NF' | sort -u >"${TELO}/fill_${tag}.tsv" || true
        else
            cp -f "${TELO}/fill_${tag}.diag.tsv" "${TELO}/fill_${tag}.tsv"
            : >"$BURY_EXTEND"  # 确保不再反复外延
        fi

        CUR_MISS=0
        if [[ -s "${TELO}/fill_${tag}.tsv" ]]; then
            CUR_MISS=$(wc -l <"${TELO}/fill_${tag}.tsv")
        fi

        if [[ "$CUR_MISS" -le 0 ]]; then
            log "无待补端，端粒轮次结束"
            cp -f "$CUR_FA" "${TELO}/after_telo_r${r}.fa"
            break
        fi
        if [[ "$PREV_MISS" -ge 0 && "$CUR_MISS" -ge "$PREV_MISS" && "$r" -gt 1 ]]; then
            log "missing 无进展（${PREV_MISS}→${CUR_MISS}），停止端粒轮次"
            cp -f "$CUR_FA" "${TELO}/after_telo_r${r}.fa"
            break
        fi

        OUT_FA="${TELO}/after_telo_r${r}.fa"
        REC="${TELO}/telo_fill_r${r}.tsv"
        run_telo_fill_round "$r" "$CUR_FA" "$TELO_INFO" "${TELO}/fill_${tag}.tsv" \
            "$TOP_LIST" "$OUT_FA" "$REC" "$TELO"
        samtools faidx "$OUT_FA"
        tail -n +2 "$REC" >>"${TELO}/telo_fill_record.tsv" || true
        CUR_FA="$OUT_FA"
        PREV_MISS=$CUR_MISS
    done
    TELO_AFTER_FA="$CUR_FA"

    # ----------------------------------------------------------
    # 4. gap
    # ----------------------------------------------------------
    log "[4/6] ONT 填 gap（最多 ${GAP_ROUNDS} 轮）..."
    CUR_FA="$TELO_AFTER_FA"
    echo -e "scaffold\tgap_index\told_start\told_end\told_length\tstatus\tread\tfill_bp\tident\tnote" >"${GAPD}/gap_fill_record.tsv"

    for ((r=1; r<=GAP_ROUNDS; r++)); do
        GAPS_NOW="${GAPD}/gaps_before_r${r}.tsv"
        bash "$CHECK_GAP" "$CUR_FA" >"$GAPS_NOW" 2>"${GAPD}/gaps_before_r${r}.summary.txt"
        gap_n=$(($(wc -l <"$GAPS_NOW") - 1))
        if [[ "$gap_n" -le 0 ]]; then
            log "gap 已清零"
            cp -f "$CUR_FA" "${GAPD}/after_gap_r${r}.fa"
            break
        fi
        if (( r % 2 == 1 )); then mode=unified; else mode=split; fi
        OUT_FA="${GAPD}/after_gap_r${r}.fa"
        REC="${GAPD}/gap_fill_r${r}.tsv"
        run_gap_fill_round "$r" "$mode" "$CUR_FA" "$GAPS_NOW" "$OUT_FA" "$REC" "$GAPD"
        samtools faidx "$OUT_FA"
        tail -n +2 "$REC" >>"${GAPD}/gap_fill_record.tsv" || true
        CUR_FA="$OUT_FA"
    done
    GAP_AFTER_FA="$CUR_FA"

    # ----------------------------------------------------------
    # 5. final + quarTeT perfect
    # ----------------------------------------------------------
    log "[5/6] 写出 final_fix/04_final（quarTeT -m ${QUARTET_MIN_REPEAT} 复检）..."
    FINAL_FA="${FIN}/${SAMPLE}_telo_gap_filled.fa"
    cp -f "$GAP_AFTER_FA" "$FINAL_FA"
    samtools faidx "$FINAL_FA"
    # 本次 final_fix 单独留底；正式记录 = 原06 ∪ 本次（同 key 以本次为准）
    cp -f "${TELO}/telo_fill_record.tsv" "${FIN}/telo_fill_record.final_fix_only.tsv"
    cp -f "${GAPD}/gap_fill_record.tsv" "${FIN}/gap_fill_record.final_fix_only.tsv"
    $PYTHON3 "$FIXLIB" merge-fill-records \
        --kind telo \
        --base "$BASE_TELO" \
        --new "${TELO}/telo_fill_record.tsv" \
        --out "${FIN}/telo_fill_record.tsv"
    $PYTHON3 "$FIXLIB" merge-fill-records \
        --kind gap \
        --base "$BASE_GAP" \
        --new "${GAPD}/gap_fill_record.tsv" \
        --out "${FIN}/gap_fill_record.tsv"
    cp -f "$DIAG_TSV" "${FIN}/orient_diagnose.tsv" 2>/dev/null || true
    cp -f "$TRIM_REC" "${FIN}/trim_record.tsv" 2>/dev/null || true
    cp -f "$TRIM_PLAN" "${FIN}/trim_plan.tsv" 2>/dev/null || true
    cp -f "$TRIM_CAND" "${FIN}/trim_candidates.tsv" 2>/dev/null || true
    cp -f "$UNRESOLVED" "${FIN}/unresolved.tsv" 2>/dev/null || true
    cp -f "$BURY_DECIDE" "${FIN}/buried_decide.tsv" 2>/dev/null || true
    [[ -s "${CLIP_TSV:-}" ]] && cp -f "$CLIP_TSV" "${FIN}/clipped_distal.tsv" || true
    [[ -s "${CLIP_FA:-}" ]] && cp -f "$CLIP_FA" "${FIN}/clipped_distal.fa" || true
    [[ -s "${TRIAL_ACCEPT:-}" ]] && cp -f "$TRIAL_ACCEPT" "${FIN}/trial_trim_accept.tsv" || true

    IFS='|' read -r TELO_INFO GAPS_TSV MISSING_TSV TOP_LIST AUDIT_TXT < <(
        refresh_quartet "$FINAL_FA" "$FIN" "final"
    )
    # refresh_quartet tag=final → ${SAMPLE}_final_Telo.* / final_gaps.tsv
    cp -f "$AUDIT_TXT" "${FIN}/final_audit_summary.txt"
    cp -f "$MISSING_TSV" "${FIN}/missing_telo_quartet.tsv"
    cp -f "$TOP_LIST" "${FIN}/top${TOP_N}.list"
    cp -f "$TOP_LIST" "${FIN}/top${TOP_N}_quartet.list"
    [[ "$GAPS_TSV" == "${FIN}/final_gaps.tsv" ]] || cp -f "$GAPS_TSV" "${FIN}/final_gaps.tsv"

    PERFECT=0
    if $PYTHON3 "$FIXLIB" check-perfect \
        --mode quartet \
        --fasta "$FINAL_FA" \
        --gaps-tsv "${FIN}/final_gaps.tsv" \
        --top-list "${FIN}/top${TOP_N}.list" \
        --top-n "$TOP_N" \
        --unresolved "${FIN}/unresolved.tsv" \
        --telo-info "${FIN}/${SAMPLE}_final_Telo.telo.info"
    then
        PERFECT=1
    fi

    {
        echo "========== ${SAMPLE} gap/telo FINAL FIX =========="
        echo "input:  $IN_FA"
        echo "ont:    $ONT_FA"
        echo "work:   $OUT_ROOT"
        echo "final:  $FINAL_FA"
        echo "perfect_criterion: quarTeT(-m $QUARTET_MIN_REPEAT) ${TOP_N} scaffolds both ends (40 tips) + 0 gaps + no unresolved"
        echo "buried: FORCE_MOTIF_TRIM=$FORCE_MOTIF_TRIM AUTO_TRIM=$AUTO_TRIM DISTAL_TRIAL_TRIM=$ALLOW_DISTAL_TRIAL_TRIM"
        echo "params: MAX_TRIM_BP=$MAX_TRIM_BP TELO_ROUNDS=$TELO_ROUNDS GAP_ROUNDS=$GAP_ROUNDS"
        echo "perfect=$PERFECT force_replace=$FORCE_REPLACE"
        echo
        echo "----- orient diagnose (init) -----"
        [[ -s "$DIAG_TSV" ]] && cat "$DIAG_TSV" || echo "(none)"
        echo
        echo "----- buried decide -----"
        [[ -s "$BURY_DECIDE" ]] && cat "$BURY_DECIDE" || echo "(none)"
        echo
        echo "----- trial_trim accept -----"
        [[ -s "${TRIAL_ACCEPT:-}" ]] && cat "$TRIAL_ACCEPT" || echo "(none)"
        echo
        echo "----- clipped distal (accepted trims) -----"
        if [[ -s "${CLIP_TSV:-}" ]]; then
            cat "$CLIP_TSV"
            echo "(fasta backup: ${FIN}/clipped_distal.fa)"
        else
            echo "(none)"
        fi
        echo
        echo "----- trim record -----"
        [[ -s "$TRIM_REC" ]] && cat "$TRIM_REC" || echo "(none)"
        echo
        echo "----- unresolved -----"
        [[ -s "$UNRESOLVED" ]] && cat "$UNRESOLVED" || echo "(none)"
        echo
        echo "----- quarTeT final (-m $QUARTET_MIN_REPEAT) -----"
        cat "${FIN}/final_audit_summary.txt"
    } | tee "${FIN}/REPORT.txt"

    # ----------------------------------------------------------
    # 6. 替换正式 04_final（仅 perfect=40端+0gap 或 FORCE_REPLACE）
    # ----------------------------------------------------------
    log "[6/6] 是否替换 ${FIN06} ..."
    if [[ "$PERFECT" -eq 1 || "$FORCE_REPLACE" == "1" ]]; then
        mkdir -p "$FIN06"
        if [[ -n "$(ls -A "$FIN06" 2>/dev/null || true)" ]]; then
            BAK="${BASE06}/04_final.bak_$(date +%Y%m%d_%H%M%S)"
            log "备份旧 04_final -> $BAK"
            mkdir -p "$BAK"
            find "$FIN06" -maxdepth 1 -type f -exec cp -a {} "$BAK"/ \;
        fi
        log "写入正式 04_final（40端+0gap perfect=$PERFECT）"
        # 只拷文件，跳过 quarTeT 产生的 tmp/ 等目录（避免 set -e 被 cp 非零退出打断）
        find "$FIN" -maxdepth 1 -type f -exec cp -f {} "$FIN06"/ \;
    else
        log "$SAMPLE：未达 quarTeT -m${QUARTET_MIN_REPEAT} 40端+0gap，保留旧 04_final；结果在 $FIN"
    fi

    log "$SAMPLE 完成。报告：${FIN}/REPORT.txt"
    echo "------------------------------------------------"
done

log "========== final_fix 全部结束 =========="
