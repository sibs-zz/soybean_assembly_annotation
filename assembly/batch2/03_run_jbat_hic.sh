#!/usr/bin/env bash
set -Eeuo pipefail

# ============================================================
# YaHS scaffold结果转换为Juicebox JBAT格式
#
# 输入文件：
#   ${sample}_yahs.bin
#   ${sample}_yahs_scaffolds_final.agp
#   ../${sample}.fasta.fai
#
# 输出文件：
#   ${sample}_JBAT.assembly
#   ${sample}_JBAT.txt
#   ${sample}_JBAT.hic
#   ${sample}_JBAT.chrom.sizes
#   ${sample}_JBAT.log
#
# 运行位置：
#   建议在YaHS输出文件所在目录运行
# ============================================================


# ======================== 配置区 =============================

JUICER_TOOLS="/tiandata2/zzh/tools/juicer_tools.1.9.9_jcuda.0.8.jar"
YAHS_JUICER="/tiandata2/zzh/tools/yahs/juicer"

ROOT="${ROOT:-/tiandata2/zzh/sjy_assembly/batch2}"

SAMPLES=(
    "BZ_primary"
    "DY11_primary"
    "TH2_primary"
    "TZX1998_primary"
    "TZX253_primary"
)

# fasta.fai所在目录
FAI_DIR="${ROOT}/01_assembly_out"
YAHS_DIR="${ROOT}/02_yahs_results"

# 输出目录
OUT_DIR="${ROOT}/03_JBAT_output"

# Java内存
JAVA_XMS="8G"
JAVA_XMX="64G"

# Java临时目录
JAVA_TMPDIR="${OUT_DIR}/java_tmp"

# 已存在结果是否覆盖：
# 0 = 不覆盖，发现已有hic时跳过
# 1 = 删除旧结果并重新生成
FORCE_OVERWRITE=0

# juicer_tools是否限制分辨率
#
# 0：使用juicer_tools默认分辨率，兼容性最高，推荐
# 1：仅生成指定分辨率，可减少运行时间和文件大小
#
# 对于JBAT人工纠错，默认分辨率通常更稳妥。
LIMIT_RESOLUTIONS=0

# 仅当LIMIT_RESOLUTIONS=1时生效
RESOLUTIONS="2500000,1000000,500000,250000,100000,50000,25000"

# ============================================================


timestamp() {
    date "+%F %T"
}


log() {
    echo "[$(timestamp)] $*"
}


die() {
    echo "[$(timestamp)] ERROR: $*" >&2
    exit 1
}


# 捕获脚本错误，显示出错行号
trap 'echo "[$(timestamp)] ERROR: 脚本在第 ${LINENO} 行运行失败。" >&2' ERR


# ====================== 环境检查 =============================

[[ -x "$YAHS_JUICER" ]] || die "找不到或不能执行YaHS juicer程序：$YAHS_JUICER"
[[ -f "$JUICER_TOOLS" ]] || die "找不到juicer_tools jar：$JUICER_TOOLS"

command -v java >/dev/null 2>&1 || die "系统中未找到java"
command -v awk  >/dev/null 2>&1 || die "系统中未找到awk"
command -v grep >/dev/null 2>&1 || die "系统中未找到grep"

mkdir -p "$OUT_DIR"
mkdir -p "$JAVA_TMPDIR"

log "Java版本："
java -version 2>&1 | head -n 3

log "输出目录：$OUT_DIR"
log "Java最大内存：$JAVA_XMX"


# ====================== 主循环 ===============================

for sample in "${SAMPLES[@]}"; do

    log "============================================================"
    log "开始处理样本：$sample"
    log "============================================================"

    # YaHS输入文件
    bin_file="${YAHS_DIR}/${sample}_yahs.bin"
    agp_file="${YAHS_DIR}/${sample}_yahs_scaffolds_final.agp"
    fai_file="${FAI_DIR}/${sample}.fasta.fai"

    # 输出文件前缀
    out_prefix="${OUT_DIR}/${sample}_JBAT"

    log_file="${out_prefix}.log"
    chrom_sizes="${out_prefix}.chrom.sizes"
    hic_file="${out_prefix}.hic"
    txt_file="${out_prefix}.txt"
    assembly_file="${out_prefix}.assembly"

    # --------------------------------------------------------
    # 1. 检查输入
    # --------------------------------------------------------

    [[ -s "$bin_file" ]] || die "$sample 缺少或文件为空：$bin_file"
    [[ -s "$agp_file" ]] || die "$sample 缺少或文件为空：$agp_file"
    [[ -s "$fai_file" ]] || die "$sample 缺少或文件为空：$fai_file"

    log "输入文件检查通过："
    ls -lh "$bin_file" "$agp_file" "$fai_file"

    # --------------------------------------------------------
    # 2. 处理已有结果
    # --------------------------------------------------------

    if [[ -s "$hic_file" && "$FORCE_OVERWRITE" -eq 0 ]]; then
        log "发现已有hic文件，跳过该样本："
        ls -lh "$hic_file"
        log "如需重新生成，将FORCE_OVERWRITE设为1。"
        continue
    fi

    if [[ "$FORCE_OVERWRITE" -eq 1 ]]; then
        log "删除该样本已有JBAT中间结果。"

        rm -f \
            "${out_prefix}.assembly" \
            "${out_prefix}.txt" \
            "${out_prefix}.hic" \
            "${out_prefix}.chrom.sizes" \
            "${out_prefix}.log"
    fi

    # --------------------------------------------------------
    # 3. YaHS juicer pre
    # --------------------------------------------------------

    log "[1/3] 运行YaHS juicer pre，生成JBAT文件。"

    "$YAHS_JUICER" pre \
        -a \
        -o "$out_prefix" \
        "$bin_file" \
        "$agp_file" \
        "$fai_file" \
        > "$log_file" 2>&1

    # 检查YaHS输出
    [[ -s "$log_file" ]] || die "$sample 未生成YaHS日志：$log_file"
    [[ -s "$txt_file" ]] || die "$sample 未生成JBAT contact文件：$txt_file"
    [[ -s "$assembly_file" ]] || die "$sample 未生成JBAT assembly文件：$assembly_file"

    log "YaHS JBAT文件生成完成："
    ls -lh "$txt_file" "$assembly_file" "$log_file"

    # --------------------------------------------------------
    # 4. 提取虚拟assembly长度
    # --------------------------------------------------------

    log "[2/3] 从PRE_C_SIZE记录中提取虚拟assembly长度。"

    # -a模式下，Juicebox人工纠错使用一个连续的虚拟染色体assembly。
    # PRE_C_SIZE行最后一列应当是assembly总长度。
    #
    # 先取最后一条PRE_C_SIZE记录，避免日志中意外出现重复记录。
    pre_c_size="$(
        awk '
            /PRE_C_SIZE/ {
                value = $NF
            }
            END {
                if (value != "") {
                    print value
                }
            }
        ' "$log_file"
    )"

    [[ -n "$pre_c_size" ]] || {
        echo "---------- YaHS日志末尾 ----------" >&2
        tail -n 50 "$log_file" >&2
        echo "----------------------------------" >&2
        die "$sample 日志中未找到PRE_C_SIZE。"
    }

    # 检查是否为纯整数
    if [[ ! "$pre_c_size" =~ ^[0-9]+$ ]]; then
        grep "PRE_C_SIZE" "$log_file" >&2 || true
        die "$sample 提取到的PRE_C_SIZE不是整数：$pre_c_size"
    fi

    # 长度不应为0
    if (( pre_c_size <= 0 )); then
        die "$sample PRE_C_SIZE异常：$pre_c_size"
    fi

    printf "assembly\t%s\n" "$pre_c_size" > "$chrom_sizes"

    log "虚拟assembly长度：$pre_c_size bp"
    cat "$chrom_sizes"

    # --------------------------------------------------------
    # 5. juicer_tools生成hic
    # --------------------------------------------------------

    log "[3/3] 运行juicer_tools pre，生成hic热图。"

    juicer_cmd=(
        java
        "-Xms${JAVA_XMS}"
        "-Xmx${JAVA_XMX}"
        "-Djava.awt.headless=true"
        "-Djava.io.tmpdir=${JAVA_TMPDIR}"
        -jar "$JUICER_TOOLS"
        pre
    )

    if [[ "$LIMIT_RESOLUTIONS" -eq 1 ]]; then
        juicer_cmd+=(
            -r "$RESOLUTIONS"
        )

        log "限制hic分辨率为：$RESOLUTIONS"
    else
        log "使用juicer_tools默认分辨率集合。"
    fi

    juicer_cmd+=(
        "$txt_file"
        "$hic_file"
        "$chrom_sizes"
    )

    # 打印实际命令，便于复现和排错
    printf '[%s] 执行命令：' "$(timestamp)"
    printf '%q ' "${juicer_cmd[@]}"
    printf '\n'

    "${juicer_cmd[@]}" \
        >> "${out_prefix}.juicer_tools.log" 2>&1

    # --------------------------------------------------------
    # 6. 结果检查
    # --------------------------------------------------------

    [[ -s "$hic_file" ]] || {
        echo "---------- juicer_tools日志末尾 ----------" >&2
        tail -n 100 "${out_prefix}.juicer_tools.log" >&2 || true
        echo "-----------------------------------------" >&2
        die "$sample hic文件未成功生成。"
    }

    hic_size_bytes=$(stat -c "%s" "$hic_file")

    # 一个正常的大豆hic一般明显大于几MB。
    # 这里只做非常宽松的完整性检查，避免产生空壳hic。
    if (( hic_size_bytes < 1048576 )); then
        die "$sample hic文件小于1 MB，可能生成不完整：$hic_file"
    fi

    log "$sample 转换完成。"
    ls -lh \
        "$assembly_file" \
        "$chrom_sizes" \
        "$hic_file" \
        "$txt_file"

    log "Juicebox中应同时加载："
    log "  热图文件：$hic_file"
    log "  组装文件：$assembly_file"

done


log "============================================================"
log "所有样本处理完成。"
log "下一步：在Juicebox中打开.hic并导入对应的.assembly文件进行人工纠错。"
log "============================================================"