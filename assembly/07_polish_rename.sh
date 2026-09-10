#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# 07_polish_rename：填补区 HiFi racon 抛光 + Wm82 改名定向
#
# 输入：06_gap_telo/${SAMPLE}/04_final/${SAMPLE}_telo_gap_filled.fa
# 输出：07_polish_rename/${SAMPLE}/03_final/${SAMPLE}_Chr.fa
#
# 用法：
#   bash 07_polish_rename.sh
# ============================================================

# ======================== 配置区 =============================

SAMPLES=(
    "KD35_primary"
#    "N2L_primary"
)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${ROOT:-/tiandata2/zzh/sjy_assembly}"
THREADS="${THREADS:-32}"
FLANK="${FLANK:-20000}"
FORCE_RERUN="${FORCE_RERUN:-0}"

REF_WM82="${REF_WM82:-/home/zzh/ref/Wm82.a2.v1.fasta}"
REF_ZH13="${REF_ZH13:-/home/zzh/ref/ZH13.v2.fasta}"

HELPER="${SCRIPT_DIR}/07_polish_rename_lib.py"
MINIMAP2="${MINIMAP2:-minimap2}"
RACON="${RACON:-/home/zzh/anaconda3/envs/qc_synteny/bin/racon}"
SEQKIT="${SEQKIT:-/home/zzh/anaconda3/envs/py37/bin/seqkit}"
PYTHON3="${PYTHON3:-python3}"

# ============================================================

timestamp() { date "+%F %T"; }
log() { echo "[$(timestamp)] $*"; }
die() { echo "[$(timestamp)] ERROR: $*" >&2; exit 1; }

hifi_for_sample() {
    local sample="$1"
    local sp="${sample%%_*}"
    case "$sp" in
        KD35) echo "${ROOT}/hifi/bc2288_hifi.filt.fasta" ;;
        N2L)  echo "${ROOT}/hifi/bc2287_hifi.filt.fasta" ;;
        *) die "未知样本前缀，无法匹配 HiFi：$sample" ;;
    esac
}

[[ -f "$HELPER" ]] || die "缺少 helper：$HELPER"
[[ -x "$RACON" || -f "$RACON" ]] || die "缺少 racon：$RACON"
command -v "$MINIMAP2" >/dev/null 2>&1 || die "缺少 minimap2"
[[ -s "$REF_WM82" ]] || die "缺少 Wm82：$REF_WM82"
[[ -s "$REF_ZH13" ]] || die "缺少 ZH13：$REF_ZH13"
[[ -x "$SEQKIT" || -f "$SEQKIT" ]] || die "缺少 seqkit：$SEQKIT"

log "THREADS=$THREADS FLANK=$FLANK FORCE_RERUN=$FORCE_RERUN"
log "RACON=$RACON"

for SAMPLE in "${SAMPLES[@]}"; do
    IN_FA="${ROOT}/06_gap_telo/${SAMPLE}/04_final/${SAMPLE}_telo_gap_filled.fa"
    TELO_REC="${ROOT}/06_gap_telo/${SAMPLE}/04_final/telo_fill_record.tsv"
    GAP_REC="${ROOT}/06_gap_telo/${SAMPLE}/04_final/gap_fill_record.tsv"
    HIFI_FA="$(hifi_for_sample "$SAMPLE")"
    OUT_ROOT="${ROOT}/07_polish_rename/${SAMPLE}"
    POL="${OUT_ROOT}/01_polish"
    MAP="${OUT_ROOT}/02_map"
    FIN="${OUT_ROOT}/03_final"

    log "============================================================"
    log "开始处理：$SAMPLE"
    log "============================================================"

    if [[ ! -s "$IN_FA" ]]; then
        log "跳过 $SAMPLE：缺少 06 最终 fa：$IN_FA"
        continue
    fi
    [[ -s "$HIFI_FA" ]] || die "$SAMPLE HiFi 不存在：$HIFI_FA"
    [[ -s "$TELO_REC" ]] || die "缺少 telo record：$TELO_REC"
    [[ -s "$GAP_REC" ]] || die "缺少 gap record：$GAP_REC"

    mkdir -p "$POL" "$MAP" "$FIN"
    if [[ "$FORCE_RERUN" == "1" ]]; then
        log "FORCE_RERUN=1：清除中间结果"
        rm -rf "${POL:?}/"* "${MAP:?}/"* "${FIN:?}/"*
    fi

    log "IN_FA=$IN_FA"
    log "HIFI=$HIFI_FA"

    # ----------------------------------------------------------
    # 1. 抽填补窗口 + HiFi racon
    # ----------------------------------------------------------
    log "[1/3] 局部 HiFi racon 抛光填补区..."

    WIN_FA="${POL}/fill_windows.fa"
    WIN_BED="${POL}/fill_windows.bed"
    WIN_PAF="${POL}/fill_windows.hifi.paf"
    WIN_IDS="${POL}/fill_windows.hit_ids.txt"
    WIN_READS="${POL}/fill_windows.hit_reads.fa"
    WIN_POL="${POL}/fill_windows.racon.fa"
    POLISHED_FA="${POL}/${SAMPLE}_polished.fa"
    POLISH_REC="${POL}/polish_record.tsv"

    if [[ ! -s "$POLISHED_FA" ]]; then
        $PYTHON3 "$HELPER" make-windows \
            --fasta "$IN_FA" \
            --telo-record "$TELO_REC" \
            --gap-record "$GAP_REC" \
            --flank "$FLANK" \
            --out-fa "$WIN_FA" \
            --out-bed "$WIN_BED"

        WIN_N=$(($(grep -c '^>' "$WIN_FA" || true)))
        if [[ "$WIN_N" -eq 0 ]]; then
            log "无 OK 填补窗口，跳过 polish，直接拷贝"
            cp -f "$IN_FA" "$POLISHED_FA"
            printf "window_id\tscaffold\tstart\tend\told_len\tnew_len\tstatus\tnote\n" >"$POLISH_REC"
        else
            MM2_LOG="${POL}/fill_windows.minimap2.log"
            if [[ ! -s "$WIN_PAF" ]]; then
                log "minimap2 map-hifi：填补窗口 <- HiFi（进度见 $MM2_LOG）..."
                "$MINIMAP2" -x map-hifi -t "$THREADS" "$WIN_FA" "$HIFI_FA" \
                    >"$WIN_PAF" 2>"$MM2_LOG"
            else
                log "复用已有窗口 PAF：$WIN_PAF"
            fi
            [[ -s "$WIN_PAF" ]] || die "窗口 PAF 为空：$WIN_PAF"

            awk '{print $1}' "$WIN_PAF" | sort -u >"$WIN_IDS"
            # hit_reads 必须与当前 PAF 一致；PAF 更新后强制重建
            if [[ ! -s "$WIN_READS" || "$WIN_PAF" -nt "$WIN_READS" ]]; then
                log "seqkit 提取命中 HiFi reads..."
                "$SEQKIT" grep -f "$WIN_IDS" "$HIFI_FA" >"${WIN_READS}.tmp"
                mv -f "${WIN_READS}.tmp" "$WIN_READS"
            fi
            [[ -s "$WIN_READS" ]] || die "未提取到命中 HiFi：$WIN_READS"

            if [[ ! -s "$WIN_POL" ]]; then
                log "racon 抛光窗口..."
                "$RACON" -t "$THREADS" "$WIN_READS" "$WIN_PAF" "$WIN_FA" \
                    >"${WIN_POL}.tmp" 2>"${POL}/fill_windows.racon.log"
                [[ -s "${WIN_POL}.tmp" ]] || die "racon 输出为空，见 ${POL}/fill_windows.racon.log"
                mv -f "${WIN_POL}.tmp" "$WIN_POL"
            else
                log "复用已有 racon 窗口：$WIN_POL"
            fi
            $PYTHON3 "$HELPER" apply-windows \
                --fasta "$IN_FA" \
                --polished-windows "$WIN_POL" \
                --windows-bed "$WIN_BED" \
                --out-fa "$POLISHED_FA" \
                --record "$POLISH_REC"
        fi
    else
        log "复用已有抛光结果：$POLISHED_FA"
    fi
    samtools faidx "$POLISHED_FA"

    # ----------------------------------------------------------
    # 2. 比对 Wm82（主）+ ZH13（校验）
    # ----------------------------------------------------------
    log "[2/3] minimap2 asm5 -> Wm82 / ZH13..."

    PAF_WM82="${MAP}/${SAMPLE}_to_Wm82.paf"
    PAF_ZH13="${MAP}/${SAMPLE}_to_ZH13.paf"

    if [[ ! -s "$PAF_WM82" ]]; then
        "$MINIMAP2" -x asm5 -t "$THREADS" --secondary=no "$REF_WM82" "$POLISHED_FA" >"$PAF_WM82"
    fi
    if [[ ! -s "$PAF_ZH13" ]]; then
        "$MINIMAP2" -x asm5 -t "$THREADS" --secondary=no "$REF_ZH13" "$POLISHED_FA" >"$PAF_ZH13"
    fi

    # ----------------------------------------------------------
    # 3. 改名定向
    # ----------------------------------------------------------
    log "[3/3] 按 Wm82 改名定向（ZH13 校验）..."

    FINAL_FA="${FIN}/${SAMPLE}_Chr.fa"
    MAP_TSV="${FIN}/chromosome_rename.tsv"
    ZH13_TSV="${FIN}/zh13_mapping.tsv"
    ZH13_CHK="${FIN}/zh13_check.tsv"

    $PYTHON3 "$HELPER" rename \
        --fasta "$POLISHED_FA" \
        --primary-paf "$PAF_WM82" \
        --check-paf "$PAF_ZH13" \
        --out-fa "$FINAL_FA" \
        --mapping-out "$MAP_TSV" \
        --check-out "$ZH13_TSV" \
        --check-report "$ZH13_CHK" \
        --restrict-to-first-20 \
        --require-complete-20

    samtools faidx "$FINAL_FA"

    {
        echo "========== ${SAMPLE} polish + rename =========="
        echo "input:     $IN_FA"
        echo "hifi:      $HIFI_FA"
        echo "polished:  $POLISHED_FA"
        echo "final:     $FINAL_FA"
        echo "wm82:      $REF_WM82"
        echo "zh13:      $REF_ZH13"
        echo
        echo "----- polish windows -----"
        if [[ -s "$POLISH_REC" ]]; then
            awk -F'\t' 'NR==1{next}{c[$7]++} END{for(k in c) print k,c[k]}' "$POLISH_REC"
        fi
        echo
        echo "----- rename (Wm82) -----"
        awk -F'\t' 'NR==1{next}{print $1" -> "$2" "$4" "$12}' "$MAP_TSV"
        echo
        echo "----- ZH13 check conflicts -----"
        if [[ -s "$ZH13_CHK" ]]; then
            awk -F'\t' 'NR==1 || $6!="YES"' "$ZH13_CHK"
        fi
    } | tee "${FIN}/REPORT.txt"

    log "$SAMPLE 完成：$FINAL_FA"
    echo "------------------------------------------------"
done

log "========== 全部样本处理结束 =========="
