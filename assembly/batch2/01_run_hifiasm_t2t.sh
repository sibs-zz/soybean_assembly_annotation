#!/usr/bin/env bash
set -euo pipefail

# batch2: HiFi + ONT-UL + Hi-C 联合组装（hifiasm）
# 样本输出：01_assembly_out/${SAMPLE}.fasta（如 BZ_primary.fasta）
#
# 依赖：hifiasm（默认 conda env: nanopore）
#
# 用法：
#   bash 01_run_hifiasm_t2t.sh
#   nohup bash 01_run_hifiasm_t2t.sh > 01_run_hifiasm_t2t.log 2>&1 &

ROOT="${ROOT:-/tiandata2/zzh/sjy_assembly/batch2}"
cd "$ROOT"

THREADS="${THREADS:-64}"
CONDA_ENV="${CONDA_ENV:-nanopore}"
CONDA_BIN="${CONDA_BIN:-/home/zzh/anaconda3/envs/${CONDA_ENV}/bin}"

SAMPLES=(
    "BZ_primary"
    "DY11_primary"
    "TH2_primary"
    "TZX1998_primary"
    "TZX253_primary"
)

HIC_DIR="${ROOT}/hic"
ONT_DIR="${ROOT}/ont"
HIFI_DIR="${ROOT}/hifi"
OUT_DIR="${ROOT}/01_assembly_out"

export PATH="${CONDA_BIN}:$PATH"
command -v hifiasm >/dev/null 2>&1 || {
    echo "ERROR: 找不到 hifiasm（期望在 ${CONDA_BIN}）" >&2
    exit 1
}
echo "[INFO] PATH tools: hifiasm=$(command -v hifiasm) ($(hifiasm -v 2>&1 | head -1))"

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

for SAMPLE in "${SAMPLES[@]}"; do
    sp="${SAMPLE%%_*}"
    HIFI="${HIFI_DIR}/${sp}_hifi.fq.gz"
    # 部分样本若只有 fasta，回退
    if [[ ! -s "$HIFI" ]]; then
        HIFI="${HIFI_DIR}/${sp}_hifi.fasta"
    fi
    ONT="${ONT_DIR}/${sp}_ONT_clean.fq.gz"
    HIC1="${HIC_DIR}/${sp}_hic_1.fq.gz"
    HIC2="${HIC_DIR}/${sp}_hic_2.fq.gz"
    PREFIX="${sp}_T2T"
    PRIMARY_FA="${SAMPLE}.fasta"

    echo "========== 开始 ${SAMPLE} 的 T2T 联合组装 =========="
    [[ -s "$HIFI" ]] || { echo "ERROR: 缺少 HiFi：$HIFI" >&2; exit 1; }
    [[ -s "$ONT"  ]] || { echo "ERROR: 缺少 ONT：$ONT" >&2; exit 1; }
    [[ -s "$HIC1" ]] || { echo "ERROR: 缺少 Hi-C R1：$HIC1" >&2; exit 1; }
    [[ -s "$HIC2" ]] || { echo "ERROR: 缺少 Hi-C R2：$HIC2" >&2; exit 1; }

    if [[ -s "$PRIMARY_FA" && "${FORCE_RERUN:-0}" != "1" ]]; then
        echo "[skip] 已有 $PRIMARY_FA"
        continue
    fi

    hifiasm \
        -o "$PREFIX" \
        -t "$THREADS" \
        --ul "$ONT" \
        --h1 "$HIC1" \
        --h2 "$HIC2" \
        "$HIFI"

    awk '/^S/{print ">"$2"\n"$3}' "${PREFIX}.hic.p_ctg.gfa" > "$PRIMARY_FA"
    ls -lh "$PRIMARY_FA"
    echo ">>> $SAMPLE 完成"
    echo "------------------------------------------------"
done

echo "========== 组装全部完成！ =========="
