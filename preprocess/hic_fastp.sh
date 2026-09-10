#!/usr/bin/env bash
set -e

echo "========== 开始进行 Hi-C 数据的极简安全质控 =========="

# 遍历你目录下的两个样本前缀
for sample in KD35 N2L; do
    echo ">>> 正在处理: $sample ..."
    
    # 运行 fastp，只进行基础去接头和质量控制
    fastp \
        -i "${sample}_R1.fq.gz" \
        -I "${sample}_R2.fq.gz" \
        -o "${sample}_clean_1.fq.gz" \
        -O "${sample}_clean_2.fq.gz" \
        --thread 16 \
        --n_base_limit 5 \
        --length_required 50 \
        --json "${sample}_hic.fastp.json" \
        --html "${sample}_hic.fastp.html"
    
    echo "  - $sample 质控完成！"
    echo "------------------------------------------------"
done

echo "========== Hi-C 预处理全部完成！ =========="
