#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# 09_run_hic_heatmap_cooler：发表级 Hi-C 热图
#
# 输入（承接 07）：
#   07_polish_rename/${SAMPLE}/03_final/${SAMPLE}_Chr.fa
# Hi-C reads：
#   hic/${sp}_hic_1.fq.gz / ${sp}_hic_2.fq.gz
#
# 流程：
#   bwa mem -> pairtools parse/sort/dedup/select
#   -> cooler matrix -> balance -> heatmap (PDF/PNG)
#
# 用法：
#   bash 09_run_hic_heatmap_cooler.sh
#   FORCE_RERUN=1 bash 09_run_hic_heatmap_cooler.sh
#   CLEAN_INTERMEDIATE_PAIRS=0 bash 09_run_hic_heatmap_cooler.sh  # 保留全部 pairs
#   （自动 conda activate hic_tools）
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
OUT_ROOT="${ROOT}/09_run_hic_heatmap"

THREADS="${THREADS:-32}"
PAIRTOOLS_THREADS="${PAIRTOOLS_THREADS:-32}"
BIN_SIZE="${BIN_SIZE:-500000}"
MIN_MAPQ="${MIN_MAPQ:-10}"
FORCE_RERUN="${FORCE_RERUN:-0}"
# 出图成功后删除 parsed/sorted/dedup/valid.allchr，保留 valid.ChrOnly
CLEAN_INTERMEDIATE_PAIRS="${CLEAN_INTERMEDIATE_PAIRS:-1}"

CONDA_ENV="${CONDA_ENV:-hic_tools}"
CONDA_HIC="${CONDA_HIC:-/home/zzh/anaconda3/envs/${CONDA_ENV}}"
PLOT_SCRIPT="${PLOT_SCRIPT:-${SCRIPT_DIR}/../plot_cooler_hic_heatmap_arial.py}"

# ============================================================

timestamp() { date "+%F %T"; }
log() { echo "[$(timestamp)] $*"; }
die() { echo "[$(timestamp)] ERROR: $*" >&2; exit 1; }

is_nonempty_file() { [[ -s "$1" ]]; }
is_valid_gzip() { [[ -s "$1" ]] && gzip -t "$1" >/dev/null 2>&1; }
is_valid_cool() { [[ -s "$1" ]] && cooler info "$1" >/dev/null 2>&1; }

count_pairs_body() {
    zcat "$1" | grep -v '^#' | wc -l
}

hic_for_sample() {
    local sample="$1"
    local sp="${sample%%_*}"
    echo "${ROOT}/hic/${sp}_hic_1.fq.gz"
    echo "${ROOT}/hic/${sp}_hic_2.fq.gz"
}

resolve_genome() {
    local sample="$1"
    local fa07="${ROOT}/07_polish_rename/${sample}/03_final/${sample}_Chr.fa"
    if [[ -s "$fa07" ]]; then
        echo "$fa07"
    else
        echo ""
    fi
}

bin_label() {
    if [[ "$BIN_SIZE" -ge 1000000 ]]; then
        echo "$((BIN_SIZE / 1000000))Mb"
    else
        echo "$((BIN_SIZE / 1000))kb"
    fi
}

ensure_plot_script() {
    if [[ -s "$PLOT_SCRIPT" ]]; then
        return 0
    fi
    log "写入绘图脚本：$PLOT_SCRIPT"
    cat > "$PLOT_SCRIPT" <<'PY'
#!/usr/bin/env python3
import argparse
import numpy as np
import cooler
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

def setup_font():
    available = {f.name for f in fm.fontManager.ttflist}
    if "Arial" in available:
        use_font = "Arial"
    elif "Liberation Sans" in available:
        use_font = "Liberation Sans"
    else:
        use_font = "DejaVu Sans"
    mpl.rcParams["font.family"] = use_font
    mpl.rcParams["font.sans-serif"] = [use_font]
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42
    mpl.rcParams["svg.fonttype"] = "none"
    mpl.rcParams["axes.unicode_minus"] = False
    print(f"[INFO] Using matplotlib font: {use_font}")

def fetch_matrix(c, chrom1, chrom2=None, use_balance=True):
    if chrom2 is None:
        chrom2 = chrom1
    try:
        return c.matrix(balance=use_balance, sparse=False).fetch(chrom1, chrom2)
    except Exception as e:
        if use_balance:
            print(f"[WARNING] Balanced fetch failed for {chrom1},{chrom2}: {e}. Fallback raw.")
            return c.matrix(balance=False, sparse=False).fetch(chrom1, chrom2)
        raise

def main():
    setup_font()
    ap = argparse.ArgumentParser()
    ap.add_argument("--cool", required=True)
    ap.add_argument("--chromsizes", required=True)
    ap.add_argument("--out_pdf", required=True)
    ap.add_argument("--out_png", required=True)
    ap.add_argument("--title", default="")
    ap.add_argument("--vmax_percentile", type=float, default=99.5)
    ap.add_argument("--figsize", default="8,8")
    ap.add_argument("--dpi", type=int, default=600)
    ap.add_argument("--cmap", default="Reds")
    ap.add_argument("--raw", action="store_true")
    args = ap.parse_args()

    use_balance = not args.raw
    chroms = []
    with open(args.chromsizes) as f:
        for line in f:
            if line.strip():
                chroms.append(line.split()[0])

    c = cooler.Cooler(args.cool)
    bin_counts = []
    for chrom in chroms:
        bin_counts.append(fetch_matrix(c, chrom, use_balance=use_balance).shape[0])

    total = sum(bin_counts)
    full = np.zeros((total, total), dtype=float)
    starts = np.cumsum([0] + bin_counts[:-1])
    ends = np.cumsum(bin_counts)

    for i, ci in enumerate(chroms):
        for j, cj in enumerate(chroms):
            block = fetch_matrix(c, ci, cj, use_balance=use_balance)
            full[starts[i]:ends[i], starts[j]:ends[j]] = block

    full = np.nan_to_num(full, nan=0.0, posinf=0.0, neginf=0.0)
    if use_balance:
        plotmat = full
        cbar_label = "Normalized Hi-C contact intensity"
    else:
        plotmat = np.log1p(full)
        cbar_label = "log1p(contact count)"

    positive = plotmat[plotmat > 0]
    vmax = float(np.percentile(positive, args.vmax_percentile)) if positive.size else 1.0
    if vmax <= 0:
        vmax = 1.0

    w, h = [float(x) for x in args.figsize.split(",")]
    fig, ax = plt.subplots(figsize=(w, h))
    im = ax.imshow(plotmat, cmap=args.cmap, vmin=0, vmax=vmax,
                   interpolation="nearest", origin="upper")
    centers = [(s + e) / 2 for s, e in zip(starts, ends)]
    ax.set_xticks(centers)
    ax.set_xticklabels(chroms, rotation=90, fontsize=9)
    ax.set_yticks(centers)
    ax.set_yticklabels(chroms, fontsize=9)
    for b in ends[:-1]:
        ax.axhline(b - 0.5, color="black", linewidth=0.15, alpha=0.65)
        ax.axvline(b - 0.5, color="black", linewidth=0.15, alpha=0.65)
    if args.title and args.title != "NONE":
        ax.set_title(args.title, fontsize=12)
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cbar.set_label(cbar_label, fontsize=9)
    plt.tight_layout()
    plt.savefig(args.out_pdf, dpi=args.dpi)
    plt.savefig(args.out_png, dpi=args.dpi)
    plt.close()

if __name__ == "__main__":
    main()
PY
    chmod +x "$PLOT_SCRIPT"
}

# -------------------- 激活 hic_tools conda 环境 --------------------
if [[ -z "${CONDA_DEFAULT_ENV:-}" || "${CONDA_DEFAULT_ENV}" != "$CONDA_ENV" ]]; then
    command -v conda >/dev/null 2>&1 || die "未找到 conda"
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$CONDA_ENV"
fi
log "conda env: ${CONDA_DEFAULT_ENV} ($(command -v python))"

# 强制使用该环境下的工具，避免 PATH 里其他版本抢先
export PATH="${CONDA_PREFIX}/bin:${PATH}"
[[ -x "${CONDA_PREFIX}/bin/bwa" ]] || die "环境 ${CONDA_ENV} 中找不到 bwa"
[[ -x "${CONDA_PREFIX}/bin/pairtools" ]] || die "环境 ${CONDA_ENV} 中找不到 pairtools"
[[ -x "${CONDA_PREFIX}/bin/cooler" ]] || die "环境 ${CONDA_ENV} 中找不到 cooler"

for exe in bwa samtools pairtools cooler python gzip awk; do
    command -v "$exe" >/dev/null 2>&1 || die "缺少工具：$exe（请确认已激活 ${CONDA_ENV}）"
done
ensure_plot_script

BIN_LABEL="$(bin_label)"
log "THREADS=$THREADS PAIRTOOLS_THREADS=$PAIRTOOLS_THREADS BIN_SIZE=$BIN_SIZE ($BIN_LABEL)"
log "MIN_MAPQ=$MIN_MAPQ FORCE_RERUN=$FORCE_RERUN"
log "OUT_ROOT=$OUT_ROOT"
log "CONDA_ENV=$CONDA_ENV PREFIX=$CONDA_PREFIX"

for SAMPLE in "${SAMPLES[@]}"; do
    sp="${SAMPLE%%_*}"
    GENOME="$(resolve_genome "$SAMPLE")"
    mapfile -t HIC_PAIR < <(hic_for_sample "$SAMPLE")
    R1="${HIC_PAIR[0]}"
    R2="${HIC_PAIR[1]}"

    SAMPLE_OUT="${OUT_ROOT}/${SAMPLE}"
    IDX="${SAMPLE_OUT}/00_index"
    PAIRS="${SAMPLE_OUT}/02_pairs"
    COOLDIR="${SAMPLE_OUT}/03_cool"
    FIG="${SAMPLE_OUT}/04_figures"
    LOGD="${SAMPLE_OUT}/logs"
    TMPD="${SAMPLE_OUT}/tmp_pairtools"
    MPL="${SAMPLE_OUT}/mpl_cache"
    mkdir -p "$IDX" "$PAIRS" "$COOLDIR" "$FIG" "$LOGD" "$TMPD" "$MPL"
    export MPLCONFIGDIR="$MPL"

    PREFIX="${sp}_hic"
    ALL_CHROMSIZES="${IDX}/${PREFIX}.all.chrom.sizes"
    CHR_CHROMSIZES="${IDX}/${PREFIX}.ChrOnly.chrom.sizes"
    PARSED_PAIRS="${PAIRS}/${PREFIX}.parsed.pairs.gz"
    SORTED_PAIRS="${PAIRS}/${PREFIX}.sorted.pairs.gz"
    DEDUP_PAIRS="${PAIRS}/${PREFIX}.dedup.pairs.gz"
    DEDUP_STATS="${PAIRS}/${PREFIX}.dedup.stats.txt"
    VALID_ALL_PAIRS="${PAIRS}/${PREFIX}.valid.allchr.pairs.gz"
    VALID_CHR_PAIRS="${PAIRS}/${PREFIX}.valid.ChrOnly.pairs.gz"
    VALID_COUNT="${PAIRS}/${PREFIX}.valid.ChrOnly.count.txt"
    DIST_DIAG="${PAIRS}/${PREFIX}.distance_diagnostics.txt"
    COOL_FILE="${COOLDIR}/${PREFIX}.${BIN_SIZE}.cool"
    BALANCE_LOG="${LOGD}/${PREFIX}.${BIN_SIZE}.cooler_balance.log"
    BALANCE_DONE="${LOGD}/${PREFIX}.${BIN_SIZE}.balance.done"

    FIG_P995_PDF="${FIG}/${PREFIX}.HiC_heatmap.${BIN_LABEL}.balanced.p99.5.pdf"
    FIG_P995_PNG="${FIG}/${PREFIX}.HiC_heatmap.${BIN_LABEL}.balanced.p99.5.png"
    FIG_P95_PDF="${FIG}/${PREFIX}.HiC_heatmap.${BIN_LABEL}.balanced.p95.pdf"
    FIG_P95_PNG="${FIG}/${PREFIX}.HiC_heatmap.${BIN_LABEL}.balanced.p95.png"
    FIG_P90_PDF="${FIG}/${PREFIX}.HiC_heatmap.${BIN_LABEL}.balanced.p90.pdf"
    FIG_P90_PNG="${FIG}/${PREFIX}.HiC_heatmap.${BIN_LABEL}.balanced.p90.png"

    log "============================================================"
    log "开始 Hi-C 热图：$SAMPLE"
    log "============================================================"

    if [[ -z "$GENOME" ]]; then
        log "跳过 $SAMPLE：缺少 07 最终 Chr.fa"
        continue
    fi
    [[ -s "$R1" && -s "$R2" ]] || die "$SAMPLE Hi-C 缺失：$R1 / $R2"
    log "GENOME=$GENOME"
    log "R1=$R1"
    log "R2=$R2"
    log "OUT=$SAMPLE_OUT"

    if [[ "$FORCE_RERUN" == "1" ]]; then
        log "FORCE_RERUN=1：清除中间 pairs/cool/figures"
        rm -rf "${PAIRS:?}/"* "${COOLDIR:?}/"* "${FIG:?}/"* "${LOGD:?}/"* "${TMPD:?}/"*
    fi

    # ---------- 1. index ----------
    log "[1] Index genome"
    if [[ ! -s "${GENOME}.fai" ]]; then
        samtools faidx "$GENOME"
    else
        log "[SKIP] faidx 已存在"
    fi
    if [[ -s "${GENOME}.bwt" || -s "${GENOME}.0123" ]]; then
        log "[SKIP] bwa index 已存在"
    else
        bwa index "$GENOME"
    fi

    # ---------- 2. chrom.sizes ----------
    log "[2] chrom.sizes"
    cut -f1,2 "${GENOME}.fai" > "$ALL_CHROMSIZES"
    awk '$1 ~ /^Chr[0-9][0-9]$/' "${GENOME}.fai" | sort -k1,1V | cut -f1,2 > "$CHR_CHROMSIZES"
    N_CHR=$(wc -l < "$CHR_CHROMSIZES")
    log "ChrOnly = ${N_CHR} 条"
    cat "$CHR_CHROMSIZES"
    [[ "$N_CHR" -eq 20 ]] || log "[WARNING] ChrOnly 不是 20 条，请检查染色体命名"

    # ---------- 3. align + parse ----------
    log "[3] bwa mem | pairtools parse"
    if is_valid_gzip "$PARSED_PAIRS"; then
        log "[SKIP] parsed pairs 已存在"
    else
        rm -f "$PARSED_PAIRS"
        bwa mem -5SP -t "$THREADS" "$GENOME" "$R1" "$R2" 2>"${LOGD}/${PREFIX}.bwa.log" \
        | pairtools parse \
            --chroms-path "$ALL_CHROMSIZES" \
            --min-mapq "$MIN_MAPQ" \
            --walks-policy 5unique \
            --drop-sam \
            --nproc-in "$PAIRTOOLS_THREADS" \
            --nproc-out "$PAIRTOOLS_THREADS" \
            -o "$PARSED_PAIRS"
    fi

    # ---------- 4. sort ----------
    log "[4] sort pairs"
    if is_valid_gzip "$SORTED_PAIRS"; then
        log "[SKIP] sorted pairs 已存在"
    else
        rm -f "$SORTED_PAIRS"
        pairtools sort --nproc "$PAIRTOOLS_THREADS" --tmpdir "$TMPD" \
            "$PARSED_PAIRS" -o "$SORTED_PAIRS"
    fi

    # ---------- 5. dedup ----------
    log "[5] dedup"
    if is_valid_gzip "$DEDUP_PAIRS" && [[ -s "$DEDUP_STATS" ]]; then
        log "[SKIP] dedup 已存在"
    else
        rm -f "$DEDUP_PAIRS"
        pairtools dedup \
            --nproc-in "$PAIRTOOLS_THREADS" \
            --nproc-out "$PAIRTOOLS_THREADS" \
            --mark-dups \
            --output-stats "$DEDUP_STATS" \
            "$SORTED_PAIRS" \
            -o "$DEDUP_PAIRS"
    fi

    # ---------- 6. valid UU ----------
    log "[6] select UU"
    if is_valid_gzip "$VALID_ALL_PAIRS"; then
        log "[SKIP] valid allchr 已存在"
    else
        rm -f "$VALID_ALL_PAIRS"
        pairtools select '(pair_type == "UU")' "$DEDUP_PAIRS" -o "$VALID_ALL_PAIRS"
    fi

    # ---------- 7. ChrOnly ----------
    log "[7] select Chr01-Chr20"
    if is_valid_gzip "$VALID_CHR_PAIRS"; then
        log "[SKIP] valid ChrOnly 已存在"
    else
        rm -f "$VALID_CHR_PAIRS"
        pairtools select \
            '(chrom1.startswith("Chr") and chrom2.startswith("Chr"))' \
            "$VALID_ALL_PAIRS" \
            -o "$VALID_CHR_PAIRS"
    fi

    # ---------- 8. format check ----------
    log "[8] pairs 格式检查"
    set +o pipefail
    FIRST_DATA_LINE=$(zcat "$VALID_CHR_PAIRS" | awk '!/^#/ {print; exit}')
    set -o pipefail
    [[ -n "$FIRST_DATA_LINE" ]] || die "ChrOnly pairs 为空：$VALID_CHR_PAIRS"
    SECOND_COL=$(echo "$FIRST_DATA_LINE" | awk '{print $2}')
    [[ "$SECOND_COL" =~ ^Chr ]] || die "pairs 第2列不是染色体名：$FIRST_DATA_LINE"

    # ---------- 9. count ----------
    log "[9] 统计 valid ChrOnly pairs"
    if [[ -s "$VALID_COUNT" && "$FORCE_RERUN" != "1" ]]; then
        log "[SKIP] count 已存在：$(cat "$VALID_COUNT")"
    else
        count_pairs_body "$VALID_CHR_PAIRS" > "$VALID_COUNT"
        log "valid ChrOnly pairs = $(cat "$VALID_COUNT")"
    fi

    # ---------- 10. distance diagnostics ----------
    log "[10] 接触距离诊断"
    if [[ -s "$DIST_DIAG" && "$FORCE_RERUN" != "1" ]]; then
        log "[SKIP] distance diagnostics 已存在"
        cat "$DIST_DIAG"
    else
        zcat "$VALID_CHR_PAIRS" | awk '
        !/^#/{
            n++; chr1=$2; pos1=$3; chr2=$4; pos2=$5;
            if(chr1==chr2){
                cis++; d=pos1-pos2; if(d<0)d=-d;
                if(d<1000) d1k++;
                else if(d<10000) d10k++;
                else if(d<100000) d100k++;
                else if(d<1000000) d1m++;
                else if(d<10000000) d10m++;
                else dlong++;
            } else { trans++; }
        }
        END{
            if(n==0){ print "total_pairs", 0; exit }
            print "total_pairs", n;
            print "cis", cis, "cis_ratio", cis/n;
            print "trans", trans, "trans_ratio", trans/n;
            print "cis_<1kb", d1k, d1k/n;
            print "cis_1kb_10kb", d10k, d10k/n;
            print "cis_10kb_100kb", d100k, d100k/n;
            print "cis_100kb_1Mb", d1m, d1m/n;
            print "cis_1Mb_10Mb", d10m, d10m/n;
            print "cis_>10Mb", dlong, dlong/n;
        }' > "$DIST_DIAG"
        cat "$DIST_DIAG"
    fi

    # ---------- 12. cooler ----------
    log "[12] cooler cload"
    if is_valid_cool "$COOL_FILE"; then
        log "[SKIP] cool 已存在：$COOL_FILE"
    else
        rm -f "$COOL_FILE"
        cooler cload pairs \
            -c1 2 -p1 3 -c2 4 -p2 5 \
            "${CHR_CHROMSIZES}:${BIN_SIZE}" \
            "$VALID_CHR_PAIRS" \
            "$COOL_FILE"
        is_valid_cool "$COOL_FILE" || die "cooler 生成失败：$COOL_FILE"
    fi

    # ---------- 13. balance ----------
    log "[13] cooler balance"
    if [[ -s "$BALANCE_DONE" && "$FORCE_RERUN" != "1" ]]; then
        log "[SKIP] balance 已完成"
    else
        if cooler balance "$COOL_FILE" >"$BALANCE_LOG" 2>&1; then
            echo "balance_done $(date '+%F %T')" > "$BALANCE_DONE"
            log "balance 完成"
        else
            log "[WARNING] balance 失败，将用 raw 矩阵绘图。见 $BALANCE_LOG"
        fi
    fi

    # ---------- 15. plot ----------
    log "[15] 绘制热图"
    run_plot_if_needed() {
        local pct="$1" pdf="$2" png="$3"
        if [[ -s "$pdf" && -s "$png" && "$FORCE_RERUN" != "1" ]]; then
            log "[SKIP] 图已存在：$pdf"
            return 0
        fi
        python "$PLOT_SCRIPT" \
            --cool "$COOL_FILE" \
            --chromsizes "$CHR_CHROMSIZES" \
            --out_pdf "$pdf" \
            --out_png "$png" \
            --title "" \
            --vmax_percentile "$pct" \
            --figsize 8,8 \
            --dpi 600
    }
    run_plot_if_needed 99.5 "$FIG_P995_PDF" "$FIG_P995_PNG"
    run_plot_if_needed 95   "$FIG_P95_PDF"  "$FIG_P95_PNG"
    run_plot_if_needed 90   "$FIG_P90_PDF"  "$FIG_P90_PNG"

    # ---------- 16. 清理中间大文件 ----------
    if [[ "$CLEAN_INTERMEDIATE_PAIRS" == "1" ]]; then
        if [[ -s "$FIG_P995_PDF" && -s "$FIG_P95_PDF" && -s "$FIG_P90_PDF" && -s "$COOL_FILE" ]]; then
            log "[16] 出图成功，删除中间 pairs（保留 valid.ChrOnly）"
            rm -f "$PARSED_PAIRS" "$SORTED_PAIRS" "$DEDUP_PAIRS" "$VALID_ALL_PAIRS"
            rm -rf "${TMPD:?}/"*
            log "已删：parsed / sorted / dedup / valid.allchr / tmp"
            log "保留：$VALID_CHR_PAIRS"
        else
            log "[16] 图或 cool 不完整，跳过清理"
        fi
    else
        log "[16] CLEAN_INTERMEDIATE_PAIRS=0，保留全部 pairs"
    fi

    {
        echo "========== ${SAMPLE} Hi-C heatmap =========="
        echo "genome: $GENOME"
        echo "R1:     $R1"
        echo "R2:     $R2"
        echo "cool:   $COOL_FILE"
        echo "bin:    $BIN_SIZE ($BIN_LABEL)"
        echo "kept_pairs: $VALID_CHR_PAIRS"
        echo
        echo "----- valid ChrOnly count -----"
        cat "$VALID_COUNT"
        echo
        echo "----- distance diagnostics -----"
        cat "$DIST_DIAG"
        echo
        echo "----- figures -----"
        echo "$FIG_P995_PDF"
        echo "$FIG_P95_PDF"
        echo "$FIG_P90_PDF"
    } | tee "${SAMPLE_OUT}/REPORT.txt"

    log "$SAMPLE 完成。报告：${SAMPLE_OUT}/REPORT.txt"
    echo "------------------------------------------------"
done

log "========== 全部样本 Hi-C 热图结束 =========="
log "输出根目录：$OUT_ROOT"
