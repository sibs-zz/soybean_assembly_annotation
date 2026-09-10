#!/usr/bin/env bash
set -e

# ================= 配置区 =================
# 设置线程数（请根据你服务器的实际 CPU 核心数拉满，比如 64 或 128）
THREADS=64

# 定义数据路径
HIC_DIR="/tiandata2/zzh/sjy_assembly/hic"
ONT_DIR="ont"
HIFI_DIR="hifi"

# ⚠️ 注意：请核对这里的 HiFi barcode 对应关系是否正确！
HIFI_N2L="${HIFI_DIR}/bc2287_hifi.filt.fastq.gz"
HIFI_KD35="${HIFI_DIR}/bc2288_hifi.filt.fastq.gz"
# ==========================================

# 创建专门的组装输出目录
mkdir -p 01_assembly_out
cd 01_assembly_out

echo "========== 开始 KD35 的 T2T 联合组装 =========="
# hifiasm 参数解析：
# -o: 输出文件前缀
# -t: 线程数
# --ul: 传入 ONT 超长数据，用于解决组装图里的死结（打通着丝粒等）
# --h1 / --h2: 传入 Hi-C 双端数据，用于分型（Phasing）和长距离挂载连接
hifiasm \
    -o KD35_T2T \
    -t $THREADS \
    --ul "../${ONT_DIR}/KD35_ONT_clean.fq.gz" \
    --h1 "${HIC_DIR}/KD35_clean_1.fq.gz" \
    --h2 "${HIC_DIR}/KD35_clean_2.fq.gz" \
    "../${HIFI_KD35}"

echo "========== 开始 N2L 的 T2T 联合组装 =========="
hifiasm \
    -o N2L_T2T \
    -t $THREADS \
    --ul "../${ONT_DIR}/N2L_ONT_clean.fq.gz" \
    --h1 "${HIC_DIR}/N2L_clean_1.fq.gz" \
    --h2 "${HIC_DIR}/N2L_clean_2.fq.gz" \
    "../${HIFI_N2L}"

# 提取 KD35 
awk '/^S/{print ">"$2"\n"$3}' KD35_T2T.hic.p_ctg.gfa > KD35_primary.fasta

# 提取 N2L 
awk '/^S/{print ">"$2"\n"$3}' N2L_T2T.hic.p_ctg.gfa > N2L_primary.fasta

echo "========== 组装全部完成！ =========="
