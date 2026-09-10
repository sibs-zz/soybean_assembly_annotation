#!/usr/bin/env bash
set -e

# ================= 配置区 =================
YAHS="/tiandata2/zzh/tools/yahs/yahs"
HIC_DIR="/tiandata2/zzh/sjy_assembly/hic"
THREADS=64  # 根据服务器核心数调整

# 这里的样本列表可以根据你 BUSCO 的评估结果进行增删
# 建议先跑每个样本最强的那个单倍型
SAMPLES=("KD35_primary" "N2L_primary") 
# ==========================================

mkdir -p 02_yahs_results
cd 02_yahs_results

for prefix in "${SAMPLES[@]}"; do
    echo "========== 正在处理样本: $prefix =========="
    
    # 确定对应的 Hi-C 原始数据前缀（KD35_hap1 对应 KD35 的 Hi-C 数据）
    sp=$(echo $prefix | cut -d'_' -f1)
    HIC1="${HIC_DIR}/${sp}_clean_1.fq.gz"
    HIC2="${HIC_DIR}/${sp}_clean_2.fq.gz"
    REF="../01_assembly_out/${prefix}.fasta"

    # 1. 构建索引
    if [ ! -f "${REF}.bwt" ]; then
        echo ">>> [1/3] 正在构建 BWA 索引..."
        bwa index "$REF"
        samtools faidx "$REF"
    fi

    # 2. Hi-C 比对 (关键：-5SP 参数)
    # YaHS 开发者建议使用 -5SP 来保留那些跨越连接点的嵌合 Read 信息
    BAM="${prefix}_hic.nameSort.markdup.bam"

    if [ ! -f "$BAM" ]; then
        echo ">>> [2/3] 正在进行 Hi-C 比对并标记 duplicates..."

        bwa mem -5SP -t "$THREADS" "$REF" "$HIC1" "$HIC2"  \
        | samblaster \
        | samtools view -@ "$THREADS" -b - \
        | samtools sort -@ "$THREADS" -n -o "$BAM" -
    fi

    # 3. 运行 YaHS 挂载
    echo ">>> [3/3] 正在运行 YaHS 进行染色体挂载..."
    $YAHS "$REF" "$BAM" -o "${prefix}_yahs"

    echo ">>> $prefix 处理完成！"
    echo "------------------------------------------------"
done

echo "========== 所有任务运行结束！ =========="
