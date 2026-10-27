#!/usr/bin/env bash
set -euo pipefail

# batch2: YaHS 染色体挂载
# 输入：01_assembly_out/${SAMPLE}.fasta + hic/${sp}_hic_{1,2}.fq.gz
# 输出：02_yahs_results/
#
# 依赖：bwa / samtools / samblaster（默认 conda env: beagle）

ROOT="${ROOT:-/tiandata2/zzh/sjy_assembly/batch2}"
cd "$ROOT"

YAHS="${YAHS:-/tiandata2/zzh/tools/yahs/yahs}"
HIC_DIR="${ROOT}/hic"
THREADS="${THREADS:-64}"
CONDA_ENV="${CONDA_ENV:-beagle}"
CONDA_BIN="${CONDA_BIN:-/home/zzh/anaconda3/envs/${CONDA_ENV}/bin}"

SAMPLES=(
    "BZ_primary"
    "DY11_primary"
    "TH2_primary"
    "TZX1998_primary"
    "TZX253_primary"
)

# 保证 samblaster/bwa/samtools 来自同一环境（避免当前 shell 缺工具）
export PATH="${CONDA_BIN}:$PATH"
for cmd in bwa samtools samblaster; do
    command -v "$cmd" >/dev/null 2>&1 || {
        echo "ERROR: 找不到 $cmd（期望在 ${CONDA_BIN}）" >&2
        exit 1
    }
done
[[ -x "$YAHS" || -f "$YAHS" ]] || { echo "ERROR: 找不到 yahs：$YAHS" >&2; exit 1; }

echo "[INFO] PATH tools: bwa=$(command -v bwa) samtools=$(command -v samtools) samblaster=$(command -v samblaster)"

mkdir -p 02_yahs_results
cd 02_yahs_results

for prefix in "${SAMPLES[@]}"; do
    echo "========== 正在处理样本: $prefix =========="

    sp="${prefix%%_*}"
    HIC1="${HIC_DIR}/${sp}_hic_1.fq.gz"
    HIC2="${HIC_DIR}/${sp}_hic_2.fq.gz"
    REF="../01_assembly_out/${prefix}.fasta"

    [[ -s "$REF" ]]  || { echo "ERROR: 缺少组装：$REF" >&2; exit 1; }
    [[ -s "$HIC1" ]] || { echo "ERROR: 缺少 Hi-C R1：$HIC1" >&2; exit 1; }
    [[ -s "$HIC2" ]] || { echo "ERROR: 缺少 Hi-C R2：$HIC2" >&2; exit 1; }

    if [[ -f "${REF}.bwt" && -f "${REF}.sa" ]]; then
        echo ">>> [1/3] 已有 BWA 索引，跳过：${REF}.bwt"
        [[ -s "${REF}.fai" ]] || samtools faidx "$REF"
    else
        echo ">>> [1/3] 正在构建 BWA 索引..."
        bwa index "$REF"
        samtools faidx "$REF"
    fi

    BAM="${prefix}_hic.nameSort.markdup.bam"

    if [[ -s "$BAM" ]]; then
        echo ">>> [2/3] 已有 BAM，跳过：$BAM"
    else
        echo ">>> [2/3] 正在进行 Hi-C 比对并标记 duplicates..."
        bwa mem -5SP -t "$THREADS" "$REF" "$HIC1" "$HIC2" \
        | samblaster \
        | samtools view -@ "$THREADS" -b - \
        | samtools sort -@ "$THREADS" -n -o "$BAM" -
    fi

    YAHS_BIN="${prefix}_yahs.bin"
    if [[ -s "$YAHS_BIN" ]]; then
        echo ">>> [3/3] 已有 YaHS 结果，跳过：$YAHS_BIN"
    else
        echo ">>> [3/3] 正在运行 YaHS 进行染色体挂载..."
        "$YAHS" "$REF" "$BAM" -o "${prefix}_yahs"
    fi

    echo ">>> $prefix 处理完成！"
    echo "------------------------------------------------"
done

echo "========== 所有任务运行结束！ =========="
