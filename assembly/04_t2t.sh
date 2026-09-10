#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Juicebox 人工纠错后，用 juicer post 生成最终 scaffolds
#
# 输入（Juicebox 导出 review.assembly 后）：
#   03_JBAT_output/${sample}_JBAT.review.assembly
#   03_JBAT_output/${sample}_JBAT.liftover.agp
#   01_assembly_out/${sample}.fasta
#
# 输出：
#   04_t2t_output/${sample}_T2T.FINAL.fa
#   04_t2t_output/${sample}_T2T.FINAL.agp
# ============================================================

# ======================== 配置区 =============================

YAHS_JUICER="/tiandata2/zzh/tools/yahs/juicer"

SAMPLES=(
    "KD35_primary"
#    "N2L_primary"
)

ASM_DIR="/tiandata2/zzh/sjy_assembly/01_assembly_out"
JBAT_DIR="/tiandata2/zzh/sjy_assembly/03_JBAT_output"
OUT_DIR="/tiandata2/zzh/sjy_assembly/04_t2t_output"

# ============================================================

mkdir -p "$OUT_DIR"

for sample in "${SAMPLES[@]}"; do
    echo "========== 正在处理样本: $sample =========="

    review="${JBAT_DIR}/${sample}_JBAT.review.assembly"
    liftover="${JBAT_DIR}/${sample}_JBAT.liftover.agp"
    contigs="${ASM_DIR}/${sample}.fasta"
    out_prefix="${OUT_DIR}/${sample}_T2T"

    [[ -s "$review" ]]   || { echo "缺少或为空: $review"; exit 1; }
    [[ -s "$liftover" ]] || { echo "缺少或为空: $liftover"; exit 1; }
    [[ -s "$contigs" ]]  || { echo "缺少或为空: $contigs"; exit 1; }

    echo ">>> juicer post -> ${out_prefix}.FINAL.{fa,agp}"
    "$YAHS_JUICER" post \
        -o "$out_prefix" \
        "$review" \
        "$liftover" \
        "$contigs"

    ls -lh "${out_prefix}.FINAL.fa" "${out_prefix}.FINAL.agp"
    echo ">>> $sample 完成"
    echo "------------------------------------------------"
done

echo "========== 全部转换完成 =========="
