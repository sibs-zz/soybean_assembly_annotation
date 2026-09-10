#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# LAI：读取 01_EDTA.sh 产物，全基因组 LAI + 自动汇总
#
# 依赖（均在 ${sp}_TE/，来自 01_EDTA.sh）：
#   -genome  ${sp}_Chr.fa          # 必须用原始 ID（Chr01…）；.mod 被 EDTA 改成 _J*，
#                                  # 与 pass.list / RepeatMasker .out 的 Chr* 对不上会空跑失败
#   -intact  ${sp}_Chr.fa.mod.EDTA.raw/LTR/${sp}_Chr.fa.mod.pass.list
#   -all     ${sp}_Chr.fa.out
#
# 用法：
#   bash 02_LAI.sh                     # 前台跑 LAI + 写 summary
#   STAGE=summary bash 02_LAI.sh       # 仅从已有 .LAI 结果汇总
#   BG=1 bash 02_LAI.sh                # 后台提交（不等待、不写 summary）
#   BG=1 WAIT=1 bash 02_LAI.sh         # 后台提交，等待完成后写 summary
#   FORCE_RERUN=1 bash 02_LAI.sh       # 强制重跑 LAI
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
BG="${BG:-0}"
WAIT="${WAIT:-0}"
# lai | summary | all
STAGE="${STAGE:-all}"

CONDA_EDTA="${CONDA_EDTA:-/home/zzh/anaconda3/envs/EDTA}"
LAI_BIN="${LAI_BIN:-${CONDA_EDTA}/share/LTR_retriever/LAI}"
SUMMARY_TSV="${SCRIPT_DIR}/LAI_summary.tsv"

# ============================================================

timestamp() { date "+%F %T"; }
log() { echo "[$(timestamp)] $*"; }
die() { echo "[$(timestamp)] ERROR: $*" >&2; exit 1; }

declare -a LAI_PIDS=()
declare -a LAI_SAMPLES=()

resolve_edta_inputs() {
    local sp="$1"
    local te_dir="${SCRIPT_DIR}/${sp}_TE"
    local genome_base="${sp}_Chr.fa"
    # LAI -genome 用原始 Chr.fa（Chr* ID），与 RM .out / pass.list 一致
    local genome="$genome_base"
    local intact="${genome_base}.mod.EDTA.raw/LTR/${genome_base}.mod.pass.list"
    local rm_out="${genome_base}.out"

    printf '%s\n%s\n%s\n%s\n' "$te_dir" "$genome" "$intact" "$rm_out"
}

edta_inputs_exist() {
    local sp="$1"
    local te_dir genome intact rm_out

    mapfile -t _paths < <(resolve_edta_inputs "$sp")
    te_dir="${_paths[0]}"
    genome="${_paths[1]}"
    intact="${_paths[2]}"
    rm_out="${_paths[3]}"

    [[ -d "$te_dir" ]] || return 1
    [[ -s "${te_dir}/${genome}" && -s "${te_dir}/${intact}" && -s "${te_dir}/${rm_out}" ]]
}

check_edta_ready() {
    local sample="$1"
    local sp="${sample%%_*}"
    local te_dir genome intact rm_out

    mapfile -t _paths < <(resolve_edta_inputs "$sp")
    te_dir="${_paths[0]}"
    genome="${_paths[1]}"
    intact="${_paths[2]}"
    rm_out="${_paths[3]}"

    if [[ ! -d "$te_dir" ]]; then
        log "跳过 $sample：缺少目录 $te_dir（请先跑 01_EDTA.sh）"
        return 1
    fi

    local missing=0
    for f in "$genome" "$intact" "$rm_out"; do
        if [[ ! -s "${te_dir}/${f}" ]]; then
            log "跳过 $sample：缺少 01_EDTA 产物 ${te_dir}/${f}"
            missing=1
        fi
    done
    [[ "$missing" -eq 0 ]]
}

clean_lai_intermediates() {
    local te_dir="$1"
    local rm_out="$2"
    rm -f "${te_dir}/${rm_out}.LAI" \
        "${te_dir}/${rm_out}.LAI".* \
        "${te_dir}/${rm_out}.LAI.LTR.fa" \
        "${te_dir}/${rm_out}.LAI.LTRlist" \
        "${te_dir}/${sp}_all.LAI.log" 2>/dev/null || true
}

clean_lai_intermediates() {
    local te_dir="$1"
    local rm_out="$2"
    local logf="$3"
    rm -f "${te_dir}/${rm_out}.LAI" \
        "${te_dir}/${rm_out}.LAI".* \
        "${logf}" 2>/dev/null || true
}

run_lai() {
    local sample="$1"
    local sp="${sample%%_*}"
    local te_dir genome intact rm_out logf lai_result

    mapfile -t _paths < <(resolve_edta_inputs "$sp")
    te_dir="${_paths[0]}"
    genome="${_paths[1]}"
    intact="${_paths[2]}"
    rm_out="${_paths[3]}"
    logf="${te_dir}/${sp}_all.LAI.log"
    lai_result="${te_dir}/${rm_out}.LAI"

    if [[ -s "$lai_result" && "$FORCE_RERUN" != "1" ]]; then
        log "已有 LAI 结果，跳过：$lai_result"
        return 0
    fi

    # 上次失败留下的空 LTR.fa / LTRlist 也清掉
    clean_lai_intermediates "$te_dir" "$rm_out" "$logf"

    local nseq
    nseq="$(grep -c '^>' "${te_dir}/${genome}" || true)"
    log "LAI 启动：$sample  scaffolds=$nseq"
    log "  genome : ${te_dir}/${genome}"
    log "  intact : ${te_dir}/${intact}"
    log "  -all   : ${te_dir}/${rm_out}"
    log "  log    : $logf"

    (
        cd "$te_dir"
        env PATH="${CONDA_EDTA}/bin:${PATH}" CONDA_PREFIX="$CONDA_EDTA" \
            "$LAI_BIN" \
            -genome "$genome" \
            -intact "$intact" \
            -all "$rm_out" \
            -t "$THREADS" \
            >"$logf" 2>&1
    )

    [[ -s "$lai_result" ]] || die "LAI 未生成结果文件：$lai_result（见 $logf）"
    log "LAI 完成：$sample -> $lai_result"
}

run_lai_bg() {
    local sample="$1"
    local sp="${sample%%_*}"
    local te_dir genome intact rm_out logf lai_result

    mapfile -t _paths < <(resolve_edta_inputs "$sp")
    te_dir="${_paths[0]}"
    genome="${_paths[1]}"
    intact="${_paths[2]}"
    rm_out="${_paths[3]}"
    logf="${te_dir}/${sp}_all.LAI.log"
    lai_result="${te_dir}/${rm_out}.LAI"

    if [[ -s "$lai_result" && "$FORCE_RERUN" != "1" ]]; then
        log "已有 LAI 结果，跳过：$lai_result"
        return 0
    fi

    clean_lai_intermediates "$te_dir" "$rm_out" "$logf"

    log "LAI 后台提交：$sample  日志：$logf"
    (
        cd "$te_dir"
        env PATH="${CONDA_EDTA}/bin:${PATH}" CONDA_PREFIX="$CONDA_EDTA" \
            "$LAI_BIN" \
            -genome "$genome" \
            -intact "$intact" \
            -all "$rm_out" \
            -t "$THREADS" \
            >"$logf" 2>&1
    ) &
    LAI_PIDS+=("$!")
    LAI_SAMPLES+=("$sample")
    log "LAI 已推入后台 PID=${LAI_PIDS[-1]}  sample=$sample"
}

parse_lai_row() {
    local lai_file="$1"
    # 输出: From To Intact Total raw_LAI LAI
    awk -F'\t' '
        NR == 1 { next }
        $1 == "whole_genome" {
            printf "%s\t%s\t%s\t%s\t%s\t%s\n", $2, $3, $4, $5, $6, $7
            exit
        }
    ' "$lai_file"
}

parse_ltr_identity() {
    local logf="$1"
    if [[ -s "$logf" ]]; then
        grep -m1 'The identity of LTR sequences:' "$logf" \
            | grep -oE '[0-9]+\.[0-9]+' \
            | head -1 || true
    fi
}

collect_sample_summary() {
    local sample="$1"
    local sp="${sample%%_*}"
    local te_dir genome rm_out logf lai_result
    local from to intact total raw_lai lai ltr_id genome_mb nseq status

    mapfile -t _paths < <(resolve_edta_inputs "$sp")
    te_dir="${_paths[0]}"
    genome="${_paths[1]}"
    rm_out="${_paths[3]}"
    logf="${te_dir}/${sp}_all.LAI.log"
    lai_result="${te_dir}/${rm_out}.LAI"

    status="ok"
    from=""; to=""; intact=""; total=""; raw_lai=""; lai=""; ltr_id=""
    genome_mb=""; nseq=""

    if [[ ! -s "$lai_result" ]]; then
        status="missing_lai"
    else
        IFS=$'\t' read -r from to intact total raw_lai lai < <(parse_lai_row "$lai_result" || true)
        ltr_id="$(parse_ltr_identity "$logf")"
        if [[ -s "${te_dir}/${genome}" ]]; then
            nseq="$(grep -c '^>' "${te_dir}/${genome}" || echo 0)"
            genome_mb="$(awk -v t="$to" 'BEGIN { if (t == "") print ""; else printf "%.2f", t/1e6 }')"
        fi
        [[ -n "$raw_lai" && -n "$lai" ]] || status="parse_error"
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$sample" "$sp" "$status" "$nseq" "$genome_mb" \
        "$intact" "$total" "$raw_lai" "$lai" "$ltr_id" "$lai_result" "$logf"
}

write_summary() {
    local tmp
    tmp="$(mktemp)"

    {
        printf 'sample\tsp\tstatus\tn_scaffolds\tgenome_mb\tintact_frac\ttotal_frac\traw_LAI\tLAI\tltr_identity_pct\tlai_result\tlai_log\n'
        for SAMPLE in "${SAMPLES[@]}"; do
            sp="${SAMPLE%%_*}"
            if edta_inputs_exist "$sp"; then
                collect_sample_summary "$SAMPLE"
            else
                printf '%s\t%s\t%s\t\t\t\t\t\t\t\t\t\n' "$SAMPLE" "$sp" "skipped_no_edta"
            fi
        done
    } >"$tmp"

    mv "$tmp" "$SUMMARY_TSV"
    log "LAI summary 已写入：$SUMMARY_TSV"
}

print_summary() {
    [[ -s "$SUMMARY_TSV" ]] || die "缺少 summary 文件：$SUMMARY_TSV"

    echo "------------------------------------------------"
    log "LAI summary"
    echo "------------------------------------------------"
    column -t -s $'\t' "$SUMMARY_TSV" 2>/dev/null || cat "$SUMMARY_TSV"
    echo "------------------------------------------------"

    awk -F'\t' '
        NR == 1 { next }
        $3 == "ok" {
            printf "  %s: raw_LAI=%s  LAI=%s  LTR_id=%s%%  genome=%s Mb\n", $1, $8, $9, $10, $5
        }
        NR > 1 && $3 != "ok" {
            printf "  %s: %s\n", $1, $3
        }
    ' "$SUMMARY_TSV"
}

[[ -x "$LAI_BIN" || -f "$LAI_BIN" ]] || die "找不到 LAI：$LAI_BIN"
[[ -d "$CONDA_EDTA" ]] || die "缺少 EDTA conda：$CONDA_EDTA"

log "SCRIPT_DIR=$SCRIPT_DIR"
log "STAGE=$STAGE THREADS=$THREADS BG=$BG WAIT=$WAIT FORCE_RERUN=$FORCE_RERUN"
log "LAI_BIN=$LAI_BIN"
log "SUMMARY_TSV=$SUMMARY_TSV"

case "$STAGE" in
    lai|all)
        for SAMPLE in "${SAMPLES[@]}"; do
            log "============================================================"
            log "样本：$SAMPLE"
            log "============================================================"

            check_edta_ready "$SAMPLE" || continue

            if [[ "$BG" == "1" ]]; then
                run_lai_bg "$SAMPLE"
            else
                run_lai "$SAMPLE"
            fi
        done

        if [[ "$BG" == "1" && ${#LAI_PIDS[@]} -gt 0 ]]; then
            if [[ "$WAIT" == "1" ]]; then
                log "等待后台 LAI 完成..."
                for i in "${!LAI_PIDS[@]}"; do
                    pid="${LAI_PIDS[$i]}"
                    samp="${LAI_SAMPLES[$i]}"
                    if wait "$pid"; then
                        log "后台 LAI 完成：$samp (PID=$pid)"
                    else
                        die "后台 LAI 失败：$samp (PID=$pid)"
                    fi
                done
            else
                log "LAI 已按样本推入后台。完成后汇总："
                echo "  WAIT=1 STAGE=summary bash ${SCRIPT_DIR}/02_LAI.sh"
                for SAMPLE in "${SAMPLES[@]}"; do
                    sp="${SAMPLE%%_*}"
                    echo "  tail -f ${SCRIPT_DIR}/${sp}_TE/${sp}_all.LAI.log"
                done
                log "========== 结束 =========="
                exit 0
            fi
        fi
        ;;
    summary)
        ;;
    *)
        die "未知 STAGE=$STAGE（lai|summary|all）"
        ;;
esac

write_summary
print_summary
log "========== 结束 =========="
