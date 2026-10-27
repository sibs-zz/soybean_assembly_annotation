#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# 10_run_centro：quarTeT CentroMiner 着丝粒定位
#
# 输入（承接 07）：
#   07_polish_rename/${SAMPLE}/03_final/${SAMPLE}_Chr.fa
#
# TE（quarTeT --TE，必须 gff3；CentroMiner 实际只用其中的 LTR）：
#   优先  annot/${sp}_TE/${sp}_Chr.fa.mod.EDTA.TEanno.gff3
#         （EDTA --anno 1 全基因组注释；当前 01_EDTA 默认未开 --anno，通常没有）
#   回退  由 annot/${sp}_TE/${sp}_Chr.fa.out（RepeatMasker）转成 *.out.gff3
#         ← 你现在有的、也是该用的同源 TE 全基因组注释
#   不要用 intact.gff3：只有完整 TE，LTR 很少，着丝粒侧翼密度会严重偏低
#
# 输出：
#   10_run_centro/${SAMPLE}/
#
# 用法：
#   bash 10_run_centro.sh                 # 默认要求 TE gff 已就绪
#   FORCE_RERUN=1 bash 10_run_centro.sh   # 覆盖已有 TRF dat 重跑
#   REQUIRE_TE=0 bash 10_run_centro.sh    # 允许无 TE（不推荐）
# ============================================================

# ======================== 配置区 =============================

SAMPLES=(
    "BZ_primary"
    "DY11_primary"
    "TH2_primary"
    "TZX1998_primary"
    "TZX253_primary"
)

ROOT="${ROOT:-/tiandata2/zzh/sjy_assembly/batch2}"
OUT_ROOT="${ROOT}/10_run_centro"
ANNOT_ROOT="${ROOT}/annot"

THREADS="${THREADS:-32}"
FORCE_RERUN="${FORCE_RERUN:-0}"
# 1=没有可用 TE gff 就跳过；0=允许裸跑基因组
REQUIRE_TE="${REQUIRE_TE:-1}"

QUARTET="${QUARTET:-/tiandata2/zzh/tools/quarTeT/quartet.py}"
CONDA_ENV="${CONDA_ENV:-quarTeTdependencies}"
# RM.out → gff3（EDTA 自带）
RM2GFF="${RM2GFF:-/home/zzh/anaconda3/envs/EDTA/share/EDTA/bin/make_gff3_with_RMout.pl}"

# ============================================================

timestamp() { date "+%F %T"; }
log() { echo "[$(timestamp)] $*"; }
die() { echo "[$(timestamp)] ERROR: $*" >&2; exit 1; }

resolve_genome() {
    local sample="$1"
    local fa07="${ROOT}/07_polish_rename/${sample}/03_final/${sample}_Chr.fa"
    if [[ -s "$fa07" ]]; then
        echo "$fa07"
    else
        echo ""
    fi
}

# 从 RepeatMasker .out 生成 quarTeT 可用的 gff3
ensure_rmout_gff() {
    local rm_out="$1"
    local gff="${rm_out}.gff3"

    if [[ -s "$gff" && "$FORCE_RERUN" != "1" ]]; then
        echo "$gff"
        return 0
    fi
    [[ -s "$rm_out" ]] || { echo ""; return 1; }
    [[ -f "$RM2GFF" ]] || die "缺少 RM→gff3 脚本：$RM2GFF"

    log "将 RepeatMasker .out 转为 gff3：$rm_out -> $gff" >&2
    (
        cd "$(dirname "$rm_out")"
        perl "$RM2GFF" "$(basename "$rm_out")"
    )
    # make_gff3_with_RMout.pl 写到 ${rm_out}.gff3
    [[ -s "$gff" ]] || die "RM.out → gff3 失败：$gff"
    echo "$gff"
}

# TE gff 解析顺序：
#   1) EDTA TEanno.gff3（官方推荐，需 --anno 1）
#   2) RepeatMasker .out 转 gff3（01_EDTA 已有产物）
resolve_te_gff() {
    local sample="$1"
    local sp="${sample%%_*}"
    local te_dir="${ANNOT_ROOT}/${sp}_TE"
    local teanno="${te_dir}/${sp}_Chr.fa.mod.EDTA.TEanno.gff3"
    local rm_out="${te_dir}/${sp}_Chr.fa.out"

    if [[ -s "$teanno" ]]; then
        echo "$teanno"
        return 0
    fi

    if [[ -s "$rm_out" ]]; then
        ensure_rmout_gff "$rm_out"
        return $?
    fi

    echo ""
}

# -------------------- 激活 conda 环境 --------------------
if [[ -z "${CONDA_DEFAULT_ENV:-}" || "${CONDA_DEFAULT_ENV}" != "$CONDA_ENV" ]]; then
    command -v conda >/dev/null 2>&1 || die "未找到 conda"
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$CONDA_ENV"
fi
log "conda env: ${CONDA_DEFAULT_ENV} ($(command -v python3))"

PYTHON3="${CONDA_PREFIX}/bin/python3"
[[ -x "$PYTHON3" ]] || die "环境 ${CONDA_ENV} 中找不到 python3：$PYTHON3"
[[ -f "$QUARTET" ]] || die "缺少 quarTeT：$QUARTET"

log "THREADS=$THREADS FORCE_RERUN=$FORCE_RERUN REQUIRE_TE=$REQUIRE_TE"
log "OUT_ROOT=$OUT_ROOT CONDA_ENV=$CONDA_ENV"

for SAMPLE in "${SAMPLES[@]}"; do
    sp="${SAMPLE%%_*}"
    GENOME="$(resolve_genome "$SAMPLE")"
    TE_GFF="$(resolve_te_gff "$SAMPLE")"
    OUT="${OUT_ROOT}/${SAMPLE}"
    PREFIX="${SAMPLE}_Centro"

    log "============================================================"
    log "开始着丝粒定位：$SAMPLE"
    log "============================================================"
    log "GENOME=$GENOME"
    log "TE_GFF=$TE_GFF"
    log "OUT=$OUT"

    if [[ -z "$GENOME" ]]; then
        log "跳过 $SAMPLE：找不到 Chr.fa（需先完成 07_polish_rename）"
        continue
    fi

    if [[ -z "$TE_GFF" ]]; then
        msg="找不到可用 TE gff（无 TEanno.gff3，也无 ${ANNOT_ROOT}/${sp}_TE/${sp}_Chr.fa.out）"
        if [[ "$REQUIRE_TE" == "1" ]]; then
            log "跳过 $SAMPLE：$msg"
            continue
        fi
        log "警告：$msg；REQUIRE_TE=0，将不带 --TE 继续"
    fi

    mkdir -p "$OUT"

    overwrite_flag=()
    if [[ "$FORCE_RERUN" == "1" ]]; then
        overwrite_flag=(--overwrite)
        log "FORCE_RERUN=1：使用 --overwrite"
    fi

    te_flag=()
    if [[ -n "$TE_GFF" ]]; then
        te_flag=(--TE "$TE_GFF")
        log "使用 EDTA TE 注释：$TE_GFF"
    fi

    log "[1] quarTeT CentroMiner -> ${OUT}/${PREFIX}.*"
    (
        cd "$OUT"
        "$PYTHON3" "$QUARTET" CentroMiner \
            -i "$GENOME" \
            "${te_flag[@]}" \
            -t "$THREADS" \
            -p "$PREFIX" \
            "${overwrite_flag[@]}"
    )

    log "完成 $SAMPLE → $OUT"
done

log "全部完成。"
