#!/usr/bin/env bash
# batch2: HiFi BAM -> fq.gz + fasta；HiC 基础 fastp。
#
# 样本: BZ / DY11 / TH2 / TZX1998 / TZX253
# TH2 有两次 HiFi（bc2219 + bc2318），分别转换后 cat 合并为 TH2_hifi.{fq.gz,fasta}
#
# 用法:
#   bash 01_prep_hifi_hic.sh
#   FORCE=1 bash 01_prep_hifi_hic.sh
#   nohup bash 01_prep_hifi_hic.sh > 01_prep_hifi_hic.log 2>&1 &
set -euo pipefail

BATCH2="${BATCH2:-/tiandata2/zzh/sjy_assembly/batch2}"
HIFI_OUT="${BATCH2}/hifi"
HIC_OUT="${BATCH2}/hic"
THREADS="${THREADS:-16}"
FORCE="${FORCE:-0}"

# 系统 /usr/local/bin/samtools 1.3.1 的 fastq 不支持 -@，优先用 beagle 1.18
SAMTOOLS="${SAMTOOLS:-/home/zzh/anaconda3/envs/beagle/bin/samtools}"
FASTP="${FASTP:-/opt/anaconda3/bin/fastp}"

SAMPLES=(BZ DY11 TH2 TZX1998 TZX253)

need_cmd() {
    command -v "$1" >/dev/null 2>&1 || { echo "ERROR: $1 not found in PATH" >&2; exit 1; }
}

need_file() {
    [ -x "$1" ] || { echo "ERROR: not executable: $1" >&2; exit 1; }
}

need_file "${SAMTOOLS}"
need_file "${FASTP}"
need_cmd pigz

mkdir -p "${HIFI_OUT}" "${HIC_OUT}"

# ---------- HiFi: BAM -> fq.gz / fasta ----------
bam_to_fq() {
    local bam="$1"
    local fq="$2"

    if [ ! -s "${bam}" ]; then
        echo "ERROR: missing BAM: ${bam}" >&2
        return 1
    fi

    if [ -s "${fq}" ] && [ "$(stat -c%s "${fq}")" -gt 1000000 ] && [ "${FORCE}" != "1" ]; then
        echo "[skip] ${fq}"
        return 0
    fi

    echo "[hifi] ${SAMTOOLS} fastq | pigz -9 -> ${fq}"
    rm -f "${fq}"
    "${SAMTOOLS}" fastq -@ "${THREADS}" "${bam}" | pigz -9 -p "${THREADS}" -c > "${fq}"
    pigz -t "${fq}"
    if [ "$(stat -c%s "${fq}")" -lt 1000000 ]; then
        echo "ERROR: output too small (likely empty): ${fq}" >&2
        return 1
    fi
}

bam_to_fa() {
    local bam="$1"
    local fa="$2"

    if [ ! -s "${bam}" ]; then
        echo "ERROR: missing BAM: ${bam}" >&2
        return 1
    fi

    if [ -s "${fa}" ] && [ "$(stat -c%s "${fa}")" -gt 1000000 ] && [ "${FORCE}" != "1" ]; then
        echo "[skip] ${fa}"
        return 0
    fi

    echo "[hifi] ${SAMTOOLS} fasta -> ${fa}"
    rm -f "${fa}"
    "${SAMTOOLS}" fasta -@ "${THREADS}" "${bam}" > "${fa}"
    if [ "$(stat -c%s "${fa}")" -lt 1000000 ]; then
        echo "ERROR: output too small (likely empty): ${fa}" >&2
        return 1
    fi
}

process_hifi_sample() {
    local sample="$1"
    local pb_dir="${BATCH2}/${sample}/pb"
    local out_fq="${HIFI_OUT}/${sample}_hifi.fq.gz"
    local out_fa="${HIFI_OUT}/${sample}_hifi.fasta"

    echo "========== HiFi: ${sample} =========="

    if [ "${sample}" = "TH2" ]; then
        local bam1="${pb_dir}/hifi_reads.bc2219.bam"
        local bam2="${pb_dir}/hifi_reads.bc2318.bam"
        local part1_fq="${HIFI_OUT}/TH2_bc2219_hifi.fq.gz"
        local part2_fq="${HIFI_OUT}/TH2_bc2318_hifi.fq.gz"
        local part1_fa="${HIFI_OUT}/TH2_bc2219_hifi.fasta"
        local part2_fa="${HIFI_OUT}/TH2_bc2318_hifi.fasta"

        bam_to_fq "${bam1}" "${part1_fq}"
        bam_to_fq "${bam2}" "${part2_fq}"
        bam_to_fa "${bam1}" "${part1_fa}"
        bam_to_fa "${bam2}" "${part2_fa}"

        if [ -s "${out_fq}" ] && [ "$(stat -c%s "${out_fq}")" -gt 1000000 ] && [ "${FORCE}" != "1" ]; then
            echo "[skip] ${out_fq}"
        else
            echo "[hifi] cat ${part1_fq} + ${part2_fq} -> ${out_fq}"
            rm -f "${out_fq}"
            cat "${part1_fq}" "${part2_fq}" > "${out_fq}"
            pigz -t "${out_fq}"
        fi

        if [ -s "${out_fa}" ] && [ "$(stat -c%s "${out_fa}")" -gt 1000000 ] && [ "${FORCE}" != "1" ]; then
            echo "[skip] ${out_fa}"
        else
            echo "[hifi] cat ${part1_fa} + ${part2_fa} -> ${out_fa}"
            rm -f "${out_fa}"
            cat "${part1_fa}" "${part2_fa}" > "${out_fa}"
        fi

        if [ "${KEEP_PARTS:-0}" != "1" ]; then
            rm -f "${part1_fq}" "${part2_fq}" "${part1_fa}" "${part2_fa}"
        fi
    else
        local bam
        bam=$(ls "${pb_dir}"/hifi_reads.bc*.bam 2>/dev/null | head -1)
        if [ -z "${bam}" ]; then
            echo "ERROR: no BAM under ${pb_dir}" >&2
            return 1
        fi
        bam_to_fq "${bam}" "${out_fq}"
        bam_to_fa "${bam}" "${out_fa}"
    fi

    echo ">>> done: ${out_fq}"
    echo ">>> done: ${out_fa}"
}

# ---------- HiC: fastp ----------
process_hic_sample() {
    local sample="$1"
    local r1="${BATCH2}/${sample}/hic/${sample}_R1.fq.gz"
    local r2="${BATCH2}/${sample}/hic/${sample}_R2.fq.gz"
    local o1="${HIC_OUT}/${sample}_hic_1.fq.gz"
    local o2="${HIC_OUT}/${sample}_hic_2.fq.gz"
    local json="${HIC_OUT}/${sample}_hic.fastp.json"
    local html="${HIC_OUT}/${sample}_hic.fastp.html"

    echo "========== HiC: ${sample} =========="

    if [ ! -s "${r1}" ] || [ ! -s "${r2}" ]; then
        echo "ERROR: missing HiC: ${r1} / ${r2}" >&2
        return 1
    fi

    if [ -s "${o1}" ] && [ -s "${o2}" ] && [ "$(stat -c%s "${o1}")" -gt 1000000 ] && [ "${FORCE}" != "1" ]; then
        echo "[skip] ${o1} / ${o2}"
        return 0
    fi

    echo "[hic] fastp -> ${o1} / ${o2}"
    rm -f "${o1}" "${o2}"
    "${FASTP}" \
        -i "${r1}" \
        -I "${r2}" \
        -o "${o1}" \
        -O "${o2}" \
        --thread "${THREADS}" \
        --n_base_limit 5 \
        --length_required 50 \
        --json "${json}" \
        --html "${html}"

    echo ">>> done: ${o1} / ${o2}"
}

# ---------- main ----------
echo "BATCH2=${BATCH2}"
echo "HIFI_OUT=${HIFI_OUT}"
echo "HIC_OUT=${HIC_OUT}"
echo "SAMTOOLS=${SAMTOOLS} ($("${SAMTOOLS}" --version | head -1))"
echo "THREADS=${THREADS} FORCE=${FORCE}"
echo ""

fail=0
for sample in "${SAMPLES[@]}"; do
    process_hifi_sample "${sample}" || fail=1
    echo ""
done

for sample in "${SAMPLES[@]}"; do
    process_hic_sample "${sample}" || fail=1
    echo ""
done

if [ "${fail}" -ne 0 ]; then
    echo "========== 存在失败样本 ==========" >&2
    exit 1
fi

echo "========== 全部完成 =========="
ls -lh "${HIFI_OUT}"/*.{fq.gz,fasta} 2>/dev/null || true
ls -lh "${HIC_OUT}"/*_hic_*.fq.gz 2>/dev/null || true
