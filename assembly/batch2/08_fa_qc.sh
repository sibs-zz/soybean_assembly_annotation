#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# 08_fa_qc：发表级基因组 QC
#   1) Merqury（HiFi k-mer）：QV + completeness + spectra-cn
#   2) BUSCO offline（embryophyta_odb10）
#   3) QUAST-LG（vs Wm82）
#   4) vs Wm82 / ZH13 点图（minimap2 asm5 + matplotlib；可选 nucmer）
#
# 输入优先：
#   07_polish_rename/${SAMPLE}/03_final/${SAMPLE}_Chr.fa
# 回退：
#   06_gap_telo/${SAMPLE}/04_final/${SAMPLE}_telo_gap_filled.fa
#
# 用法：
#   bash 08_fa_qc.sh
#   FORCE_RERUN=1 bash 08_fa_qc.sh
#   USE_NUCMER=1 bash 08_fa_qc.sh          # 额外跑 nucmer 点图（更慢）
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
OUT_ROOT="${ROOT}/08_fa_qc"
THREADS="${THREADS:-64}"
FORCE_RERUN="${FORCE_RERUN:-0}"
USE_NUCMER="${USE_NUCMER:-0}"

# 开关：1=跑，0=跳过
RUN_MERQURY="${RUN_MERQURY:-1}"
RUN_BUSCO="${RUN_BUSCO:-1}"
RUN_QUAST="${RUN_QUAST:-1}"
RUN_DOTPLOT="${RUN_DOTPLOT:-1}"

REF_WM82="/home/zzh/ref/Wm82.a2.v1.fasta"
REF_ZH13="/home/zzh/ref/ZH13.v2.fasta"

# Merqury / meryl
MERYL_K="${MERYL_K:-21}"
CONDA_MERQURY="/home/zzh/anaconda3/envs/qc_merqury"
MERQURY_HOME="${CONDA_MERQURY}/share/merqury"

# BUSCO（offline，参考 forother/fc/run.sh）
CONDA_BUSCO="/home/zzh/anaconda3/envs/qc_busco"
BUSCO_LINEAGE="/tiandata2/zzh/forother/fc/busco_downloads/lineages/embryophyta_odb10"
BUSCO_DOWNLOAD_PATH="/tiandata2/zzh/forother/fc/busco_downloads"

# QUAST
CONDA_QUAST="/home/zzh/anaconda3/envs/qc_quast"

# synteny / plot
CONDA_SYNTENY="/home/zzh/anaconda3/envs/qc_synteny"
MINIMAP2="${MINIMAP2:-${CONDA_SYNTENY}/bin/minimap2}"
SEQKIT="${SEQKIT:-/home/zzh/anaconda3/envs/py37/bin/seqkit}"
PLOT_PYTHON="${PLOT_PYTHON:-/home/zzh/anaconda3/envs/py37/bin/python}"
HELPER="${SCRIPT_DIR}/../08_fa_qc_lib.py"

# ============================================================

timestamp() { date "+%F %T"; }
log() { echo "[$(timestamp)] $*"; }
die() { echo "[$(timestamp)] ERROR: $*" >&2; exit 1; }

hifi_for_sample() {
    local sample="$1"
    local sp="${sample%%_*}"
    local fq="${ROOT}/hifi/${sp}_hifi.fq.gz"
    local fa="${ROOT}/hifi/${sp}_hifi.fasta"
    if [[ -s "$fq" ]]; then
        echo "$fq"
    elif [[ -s "$fa" ]]; then
        echo "$fa"
    else
        die "未知样本或缺少 HiFi：$sample（期望 $fq 或 $fa）"
    fi
}

resolve_input_fa() {
    local sample="$1"
    local fa07="${ROOT}/07_polish_rename/${sample}/03_final/${sample}_Chr.fa"
    local fa06="${ROOT}/06_gap_telo/${sample}/04_final/${sample}_telo_gap_filled.fa"
    if [[ -s "$fa07" ]]; then
        echo "$fa07"
    elif [[ -s "$fa06" ]]; then
        echo "$fa06"
    else
        echo ""
    fi
}

with_env_path() {
    # 把 env/bin 放到 PATH 最前，避免 /opt/anaconda3/bin 抢 python
    local env_prefix="$1"
    shift
    env PATH="${env_prefix}/bin:${PATH}" \
        CONDA_PREFIX="${env_prefix}" \
        MERQURY="${MERQURY_HOME}" \
        "$@"
}

extract_chr_fa() {
    local in_fa="$1" out_fa="$2" pattern="$3"
    if [[ -s "$out_fa" && -s "${out_fa}.fai" && "$FORCE_RERUN" != "1" ]]; then
        return 0
    fi
    "$SEQKIT" grep -r -p "$pattern" "$in_fa" >"${out_fa}.tmp"
    mv -f "${out_fa}.tmp" "$out_fa"
    samtools faidx "$out_fa"
}

run_merqury() {
    local sample="$1" fa="$2" hifi="$3" outdir="$4"
    local sp="${sample%%_*}"
    local meryl_root="${OUT_ROOT}/meryl_dbs"
    local meryl_db="${meryl_root}/${sp}_hifi_k${MERYL_K}.meryl"
    local prefix="${sample}_merqury"
    mkdir -p "$meryl_root" "$outdir"

    if [[ -s "${outdir}/${prefix}.qv" && "$FORCE_RERUN" != "1" ]]; then
        log "已有 Merqury QV，跳过：${outdir}/${prefix}.qv"
        return 0
    fi

    if [[ ! -d "$meryl_db" ]]; then
        log "meryl count k=${MERYL_K} <- $hifi"
        with_env_path "$CONDA_MERQURY" \
            meryl k="$MERYL_K" count threads="$THREADS" "$hifi" output "$meryl_db"
    else
        log "复用已有 meryl：$meryl_db"
    fi

    (
        cd "$outdir"
        # 清理旧中间，保留目录
        find . -maxdepth 1 -type f -name "${prefix}*" -delete 2>/dev/null || true
        find . -maxdepth 1 -type d -name '*.meryl' ! -name '.' -exec rm -rf {} + 2>/dev/null || true
        local fa_link="${sample}.fa"
        ln -sfn "$fa" "$fa_link"
        log "merqury.sh ..."
        with_env_path "$CONDA_MERQURY" \
            merqury.sh "$meryl_db" "$fa_link" "$prefix"
    )
    log "Merqury 完成：$outdir"
    ls -lh "$outdir"/*.qv "$outdir"/*.completeness.stats 2>/dev/null || true
}

run_busco() {
    local sample="$1" fa="$2" outdir="$3"
    local run_name="${sample}_busco"
    local summary="" rc=0
    mkdir -p "$outdir"

    if [[ "$FORCE_RERUN" == "1" ]]; then
        rm -rf "${outdir}/${run_name}" "${outdir}/busco.log"
    fi
    # 以 short_summary*.txt 为准；JSON 写挂也不挡下游
    summary="$(find "${outdir}/${run_name}" -name 'short_summary*.txt' 2>/dev/null | head -1 || true)"
    if [[ -n "$summary" && -s "$summary" ]]; then
        log "已有 BUSCO short_summary，跳过：$summary"
        return 0
    fi

    [[ -d "$BUSCO_LINEAGE" ]] || die "缺少 BUSCO lineage：$BUSCO_LINEAGE"

    log "BUSCO offline <- $BUSCO_LINEAGE"
    set +e
    (
        cd "$outdir"
        with_env_path "$CONDA_BUSCO" \
            python "${CONDA_BUSCO}/bin/busco" \
                -m genome \
                -i "$fa" \
                -l "$BUSCO_LINEAGE" \
                --offline \
                --download_path "$BUSCO_DOWNLOAD_PATH" \
                -o "$run_name" \
                -c "$THREADS"
    ) >"${outdir}/busco.log" 2>&1
    rc=$?
    set -e

    summary="$(find "${outdir}/${run_name}" -name 'short_summary*.txt' 2>/dev/null | head -1 || true)"
    if [[ -n "$summary" && -s "$summary" ]]; then
        if [[ "$rc" -ne 0 ]]; then
            log "BUSCO 退出码=$rc（可忽略：已有 short_summary），继续下游：$summary"
        else
            log "BUSCO 完成：见 ${outdir}/busco.log"
        fi
        find "$outdir" -name 'short_summary*.txt' -print -exec ls -lh {} \; 2>/dev/null | head -10 || true
        return 0
    fi
    die "BUSCO 失败且无 short_summary.txt（rc=$rc）：见 ${outdir}/busco.log"
}

run_quast() {
    local sample="$1" fa="$2" ref="$3" outdir="$4"
    mkdir -p "$outdir"
    if [[ -s "${outdir}/report.tsv" && "$FORCE_RERUN" != "1" ]]; then
        log "已有 QUAST 结果，跳过：$outdir"
        return 0
    fi
    [[ -s "$ref" ]] || die "缺少 QUAST 参考：$ref"

    log "QUAST-LG vs $(basename "$ref")"
    with_env_path "$CONDA_QUAST" \
        python "${CONDA_QUAST}/bin/quast-lg.py" \
            -o "$outdir" \
            -r "$ref" \
            -t "$THREADS" \
            -l "$sample" \
            --min-contig 1000 \
            "$fa" \
            >"${outdir}/quast.log" 2>&1

    log "QUAST 完成：$outdir"
    ls -lh "${outdir}/report.tsv" "${outdir}/report.txt" 2>/dev/null || true
}

run_dotplot_minimap() {
    local sample="$1" query_chr="$2" ref_chr="$3" ref_label="$4" outdir="$5"
    local paf="${outdir}/${sample}_vs_${ref_label}.paf"
    local png="${outdir}/${sample}_vs_${ref_label}.dotplot.png"

    mkdir -p "$outdir"
    if [[ -s "$png" && "$FORCE_RERUN" != "1" ]]; then
        log "已有点图，跳过：$png"
        return 0
    fi

    log "minimap2 asm5：${sample} <- ${ref_label}"
    "$MINIMAP2" -x asm5 -t "$THREADS" --secondary=no \
        "$ref_chr" "$query_chr" >"$paf"

    "$PLOT_PYTHON" "$HELPER" dotplot \
        --paf "$paf" \
        --query-fai "${query_chr}.fai" \
        --ref-fai "${ref_chr}.fai" \
        --out "$png" \
        --title "${sample} vs ${ref_label}" \
        --query-label "$sample" \
        --ref-label "$ref_label" \
        --min-alen 20000 \
        --min-mapq 20

    log "点图完成：$png"
}

run_dotplot_nucmer() {
    local sample="$1" query_chr="$2" ref_chr="$3" ref_label="$4" outdir="$5"
    local prefix="${outdir}/${sample}_vs_${ref_label}_nucmer"
    local png="${prefix}.dotplot.png"

    mkdir -p "$outdir"
    if [[ -s "$png" && "$FORCE_RERUN" != "1" ]]; then
        log "已有 nucmer 点图，跳过：$png"
        return 0
    fi

    log "nucmer：${sample} <- ${ref_label}（较慢）"
    with_env_path "$CONDA_SYNTENY" \
        nucmer --maxmatch -l 100 -c 500 -t "$THREADS" \
            -p "$prefix" "$ref_chr" "$query_chr"

    with_env_path "$CONDA_SYNTENY" \
        delta-filter -1 -q -r "${prefix}.delta" >"${prefix}.1delta"

    with_env_path "$CONDA_SYNTENY" \
        show-coords -T -H -c -l "${prefix}.1delta" >"${prefix}.coords"

    "$PLOT_PYTHON" "$HELPER" dotplot \
        --coords "${prefix}.coords" \
        --query-fai "${query_chr}.fai" \
        --ref-fai "${ref_chr}.fai" \
        --out "$png" \
        --title "${sample} vs ${ref_label} (nucmer)" \
        --query-label "$sample" \
        --ref-label "$ref_label" \
        --min-alen 20000

    log "nucmer 点图完成：$png"
}

# -------------------- 预检 --------------------
[[ -f "$HELPER" ]] || die "缺少 helper：$HELPER"
[[ -x "$SEQKIT" || -f "$SEQKIT" ]] || die "缺少 seqkit：$SEQKIT"
[[ -x "$PLOT_PYTHON" ]] || die "缺少绘图 python：$PLOT_PYTHON"
[[ -x "$MINIMAP2" ]] || command -v "$MINIMAP2" >/dev/null 2>&1 || die "缺少 minimap2：$MINIMAP2"
[[ -d "$CONDA_MERQURY" ]] || die "缺少 conda env：$CONDA_MERQURY"
[[ -d "$CONDA_BUSCO" ]] || die "缺少 conda env：$CONDA_BUSCO"
[[ -d "$CONDA_QUAST" ]] || die "缺少 conda env：$CONDA_QUAST"
[[ -d "$CONDA_SYNTENY" ]] || die "缺少 conda env：$CONDA_SYNTENY"
[[ -s "$REF_WM82" ]] || die "缺少 Wm82：$REF_WM82"
[[ -s "$REF_ZH13" ]] || die "缺少 ZH13：$REF_ZH13"
[[ -x "${CONDA_MERQURY}/bin/meryl" ]] || die "缺少 meryl"
[[ -f "${CONDA_BUSCO}/bin/busco" ]] || die "缺少 busco"
[[ -f "${CONDA_QUAST}/bin/quast-lg.py" ]] || die "缺少 quast-lg.py"
command -v samtools >/dev/null 2>&1 || die "缺少 samtools"

mkdir -p "$OUT_ROOT"

# 参考染色体子集（缓存）
REF_CACHE="${OUT_ROOT}/ref_chr"
mkdir -p "$REF_CACHE"
WM82_CHR="${REF_CACHE}/Wm82.Chr01-20.fa"
ZH13_CHR="${REF_CACHE}/ZH13.Chr01-20.fa"
extract_chr_fa "$REF_WM82" "$WM82_CHR" '^Wm82\.Chr([0-9]+)$'
extract_chr_fa "$REF_ZH13" "$ZH13_CHR" '^ZH13\.Chr([0-9]+)$'

log "THREADS=$THREADS FORCE_RERUN=$FORCE_RERUN USE_NUCMER=$USE_NUCMER"
log "RUN: merqury=$RUN_MERQURY busco=$RUN_BUSCO quast=$RUN_QUAST dotplot=$RUN_DOTPLOT"
log "SAMPLES: ${SAMPLES[*]}"

for SAMPLE in "${SAMPLES[@]}"; do
    log "============================================================"
    log "开始 QC：$SAMPLE"
    log "============================================================"

    IN_FA="$(resolve_input_fa "$SAMPLE")"
    if [[ -z "$IN_FA" ]]; then
        log "跳过 $SAMPLE：缺少 07 Chr.fa 与 06 final fa"
        continue
    fi
    HIFI_FA="$(hifi_for_sample "$SAMPLE")"
    [[ -s "$HIFI_FA" ]] || die "$SAMPLE HiFi 不存在：$HIFI_FA"

    SAMPLE_OUT="${OUT_ROOT}/${SAMPLE}"
    PREP="${SAMPLE_OUT}/00_prep"
    MERQ="${SAMPLE_OUT}/01_merqury"
    BUSC="${SAMPLE_OUT}/02_busco"
    QUAS="${SAMPLE_OUT}/03_quast_Wm82"
    DOT="${SAMPLE_OUT}/04_dotplot"
    FIN="${SAMPLE_OUT}/05_summary"
    mkdir -p "$PREP" "$MERQ" "$BUSC" "$QUAS" "$DOT" "$FIN"

    if [[ "$FORCE_RERUN" == "1" ]]; then
        log "FORCE_RERUN=1：清除 $SAMPLE 中间结果（保留共享 meryl_dbs）"
        rm -rf "${MERQ:?}/"* "${BUSC:?}/"* "${QUAS:?}/"* "${DOT:?}/"* "${FIN:?}/"* "${PREP:?}/"*
    fi

    log "IN_FA=$IN_FA"
    log "HIFI=$HIFI_FA"
    log "OUT=$SAMPLE_OUT"

    WORK_FA="${PREP}/${SAMPLE}.fa"
    CHR_FA="${PREP}/${SAMPLE}.Chr01-20.fa"
    cp -f "$IN_FA" "$WORK_FA"
    samtools faidx "$WORK_FA"

    # 用 fai 判定，避免 seqkit|grep 管道在部分环境下假阴性
    if [[ -s "${WORK_FA}.fai" ]] && awk '$1 ~ /^Chr[0-9]+$/ {found=1} END{exit !found}' "${WORK_FA}.fai"; then
        extract_chr_fa "$WORK_FA" "$CHR_FA" '^Chr([0-9]+)$'
        log "点图用染色体子集：$(awk '$1 ~ /^Chr[0-9]+$/' "${CHR_FA}.fai" | wc -l) 条 Chr*"
    else
        log "未检测到 Chr* 命名，点图改用全序列"
        cp -f "$WORK_FA" "$CHR_FA"
        samtools faidx "$CHR_FA"
    fi

    STATS_TSV="${PREP}/seqkit_stats.tsv"
    "$SEQKIT" stats -a -T "$WORK_FA" >"$STATS_TSV"
    log "seqkit stats："
    cat "$STATS_TSV"

    if [[ "$RUN_MERQURY" == "1" ]]; then
        run_merqury "$SAMPLE" "$WORK_FA" "$HIFI_FA" "$MERQ"
    else
        log "跳过 Merqury"
    fi

    if [[ "$RUN_BUSCO" == "1" ]]; then
        run_busco "$SAMPLE" "$WORK_FA" "$BUSC"
    else
        log "跳过 BUSCO"
    fi

    if [[ "$RUN_QUAST" == "1" ]]; then
        run_quast "$SAMPLE" "$WORK_FA" "$WM82_CHR" "$QUAS"
        run_quast "$SAMPLE" "$WORK_FA" "$ZH13_CHR" "${SAMPLE_OUT}/03_quast_ZH13"
    else
        log "跳过 QUAST"
    fi

    if [[ "$RUN_DOTPLOT" == "1" ]]; then
        run_dotplot_minimap "$SAMPLE" "$CHR_FA" "$WM82_CHR" "Wm82" "${DOT}/Wm82"
        run_dotplot_minimap "$SAMPLE" "$CHR_FA" "$ZH13_CHR" "ZH13" "${DOT}/ZH13"
        if [[ "$USE_NUCMER" == "1" ]]; then
            run_dotplot_nucmer "$SAMPLE" "$CHR_FA" "$WM82_CHR" "Wm82" "${DOT}/Wm82_nucmer"
            run_dotplot_nucmer "$SAMPLE" "$CHR_FA" "$ZH13_CHR" "ZH13" "${DOT}/ZH13_nucmer"
        fi
    else
        log "跳过点图"
    fi

    log "汇总 REPORT..."
    "$PLOT_PYTHON" "$HELPER" summarize \
        --sample "$SAMPLE" \
        --stats "$STATS_TSV" \
        --merqury-dir "$MERQ" \
        --busco-dir "$BUSC" \
        --quast-dir "$QUAS" \
        --out-tsv "${FIN}/qc_metrics.tsv" \
        --out-txt "${FIN}/REPORT.txt"

    {
        echo -e "key\tvalue"
        echo -e "input_fa\t$IN_FA"
        echo -e "hifi\t$HIFI_FA"
        echo -e "wm82_chr\t$WM82_CHR"
        echo -e "zh13_chr\t$ZH13_CHR"
        echo -e "merqury_dir\t$MERQ"
        echo -e "busco_dir\t$BUSC"
        echo -e "quast_dir\t$QUAS"
        echo -e "dotplot_dir\t$DOT"
    } >"${FIN}/paths.tsv"

    log "$SAMPLE QC 完成"
    log "  REPORT : ${FIN}/REPORT.txt"
    log "  metrics: ${FIN}/qc_metrics.tsv"
    echo "------------------------------------------------"
done

log "========== 全部样本 QC 结束 =========="
log "输出根目录：$OUT_ROOT"
