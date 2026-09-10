#!/usr/bin/env bash
# 从 PacBio HiFi BAM 导出 fastq.gz 和 fasta，供 hifiasm 等下游使用。
#
# 样本对应（与 01_run_hifiasm_t2t.sh 一致）：
#   bc2287 -> N2L
#   bc2288 -> KD35
#
# 用法:
#   bash bam2fastq_fasta.sh
#   FORCE=1 bash bam2fastq_fasta.sh
#   nohup bash bam2fastq_fasta.sh > bam2fastq_fasta.log 2>&1 &
set -euo pipefail

SRC_DIR="/tiandata1/Original_sequence_data/re-seq/SJY_assembly/m84074_260211_170138_s2/hifi_reads"
RUN_PREFIX="m84074_260211_170138_s2.hifi_reads"
OUT_DIR="/tiandata2/zzh/sjy_assembly/hifi"
THREADS="${THREADS:-32}"
FORCE="${FORCE:-0}"

mkdir -p "${OUT_DIR}"

# barcode -> 样本名（仅用于日志）
declare -A SAMPLE_NAME=(
    [bc2287]="N2L"
    [bc2288]="KD35"
)

need_cmd() {
    command -v "$1" >/dev/null 2>&1 || { echo "ERROR: $1 not found in PATH" >&2; exit 1; }
}

need_cmd samtools
need_cmd pigz

convert_one() {
    local bc="$1"
    local sample="${SAMPLE_NAME[$bc]}"
    local bam="${SRC_DIR}/${RUN_PREFIX}.${bc}.bam"
    local fq="${OUT_DIR}/${bc}_hifi.filt.fastq.gz"
    local fa="${OUT_DIR}/${bc}_hifi.filt.fasta"
    local stats="${OUT_DIR}/${bc}_hifi.stats"

    echo "========== ${sample} (${bc}) =========="

    if [ ! -s "${bam}" ]; then
        echo "ERROR: missing BAM: ${bam}" >&2
        return 1
    fi

    # BAM -> FASTQ.GZ（pigz -9 确保真正 gzip 压缩）
    if [ -s "${fq}" ] && [ "${FORCE}" != "1" ]; then
        echo "[skip] FASTQ exists: ${fq}"
    else
        echo "[1/3] samtools fastq -> ${fq}"
        rm -f "${fq}"
        samtools fastq -@ "${THREADS}" "${bam}" | pigz -9 -p "${THREADS}" -c > "${fq}"
        pigz -t "${fq}"
    fi

    # BAM -> FASTA
    if [ -s "${fa}" ] && [ "${FORCE}" != "1" ]; then
        echo "[skip] FASTA exists: ${fa}"
    else
        echo "[2/3] samtools fasta -> ${fa}"
        rm -f "${fa}"
        samtools fasta -@ "${THREADS}" "${bam}" > "${fa}"
    fi

    # 统计
    echo "[3/3] seqkit stats -> ${stats}"
    if command -v seqkit >/dev/null 2>&1; then
        seqkit stats -T "${fq}" > "${stats}"
        cat "${stats}"
    else
        echo "WARN: seqkit not found, skip stats"
        samtools idxstats "${bam}" | awk '{s+=$3} END{print "reads:", s}'
    fi

    echo ">>> done: ${fq}"
    echo ">>> done: ${fa}"
}

main() {
    echo "SRC: ${SRC_DIR}"
    echo "OUT: ${OUT_DIR}"
    echo "THREADS=${THREADS} FORCE=${FORCE}"
    echo ""

    local fail=0
    for bc in bc2287 bc2288; do
        convert_one "${bc}" || fail=1
        echo ""
    done

    if [ "${fail}" -ne 0 ]; then
        echo "========== 存在失败样本 ==========" >&2
        exit 1
    fi

    echo "========== 全部完成 =========="
    ls -lah "${OUT_DIR}"/*.{fq.gz,fasta,stats} 2>/dev/null || ls -lah "${OUT_DIR}/"
}

main "$@"
