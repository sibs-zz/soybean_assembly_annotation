#!/usr/bin/env bash
set -euo pipefail

# 从盐胁迫 RNA-seq 里取出可用于基因组注释/组装的文库：
#   - 全部 CK（对照，非 salt）
#   - TZX1998-Flo、TZX-253-Pod（补组织）
# 质控参数与 rna.sh 一致，输出到 rna_clean/

RAW_ROOT="/tiandata4/sjy_salt_rnaseq/raw"
OUT_DIR="/tiandata2/zzh/sjy_assembly/rna_clean"
THREADS="${THREADS:-16}"
FORCE_RERUN="${FORCE_RERUN:-0}"

mkdir -p "$OUT_DIR"

FASTP="$(command -v fastp || true)"
[[ -n "$FASTP" ]] || FASTP="/opt/anaconda3/bin/fastp"
[[ -x "$FASTP" ]] || { echo "ERROR: 找不到 fastp" >&2; exit 1; }

# 输出样本名：TZX-253-Pod -> TZX253-Pod，与已有 TZX253-FL 一致
out_name_for() {
    local name="$1"
    if [[ "$name" == "TZX-253-Pod" ]]; then
        echo "TZX253-Pod"
    else
        echo "$name"
    fi
}

is_usable() {
    local name="$1"
    [[ "$name" == CK-* || "$name" == "TZX1998-Flo" || "$name" == "TZX-253-Pod" ]]
}

echo "========== 盐胁迫批次：筛选可用 RNA-seq 并质控 =========="
echo "RAW=$RAW_ROOT"
echo "OUT=$OUT_DIR"
echo "THREADS=$THREADS"
echo "fastp=$FASTP"

n_ok=0
n_skip=0
n_fail=0

mapfile -t SAMPLE_DIRS < <(find "$RAW_ROOT" -path '*/RawData/*' -mindepth 3 -maxdepth 3 -type d | sort)

for d in "${SAMPLE_DIRS[@]}"; do
    name="$(basename "$d")"
    is_usable "$name" || continue

    r1="$(find "$d" -maxdepth 1 -name '*_1.fq.gz' | head -1)"
    r2="$(find "$d" -maxdepth 1 -name '*_2.fq.gz' | head -1)"
    if [[ -z "$r1" || -z "$r2" || ! -s "$r1" || ! -s "$r2" ]]; then
        echo "ERROR: $name 缺少成对 fq.gz" >&2
        n_fail=$((n_fail + 1))
        continue
    fi

    sample="$(out_name_for "$name")"
    o1="${OUT_DIR}/${sample}_clean_1.fq.gz"
    o2="${OUT_DIR}/${sample}_clean_2.fq.gz"

    if [[ -s "$o1" && -s "$o2" && "$FORCE_RERUN" != "1" ]]; then
        echo ">>> 跳过已有：$sample"
        n_skip=$((n_skip + 1))
        continue
    fi

    echo ">>> 正在处理转录组: $sample ..."
    echo "    R1=$r1"
    echo "    R2=$r2"

    "$FASTP" \
        -i "$r1" \
        -I "$r2" \
        -o "$o1" \
        -O "$o2" \
        --thread "$THREADS" \
        --detect_adapter_for_pe \
        --length_required 50 \
        --json "${OUT_DIR}/${sample}.fastp.json" \
        --html "${OUT_DIR}/${sample}.fastp.html"

    echo "  - $sample 质控完成！"
    echo "------------------------------------------------"
    n_ok=$((n_ok + 1))
done

echo "========== 质控结束：新跑=$n_ok  已存在跳过=$n_skip  失败=$n_fail =========="
echo "输出目录：$OUT_DIR"
