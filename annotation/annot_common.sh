#!/usr/bin/env bash
# 由 04–09 脚本 source。不要直接跑。
# 不 activate conda；按步骤把对应 env 的 bin 放到 PATH 最前。

ANACONDA="${ANACONDA:-/home/zzh/anaconda3/envs}"

ENV_ALIGNER="${ENV_ALIGNER:-${ANACONDA}/aligner}"
ENV_ANNP="${ENV_ANNP:-${ANACONDA}/annp}"
ENV_BRAKER="${ENV_BRAKER:-${ANACONDA}/braker3}"
ENV_GMAP="${ENV_GMAP:-${ANACONDA}/annot_gmap}"
ENV_TRINITY="${ENV_TRINITY:-${ANACONDA}/annot_trinity}"
ENV_TRANSDECODER="${ENV_TRANSDECODER:-${ANACONDA}/annot_transdecoder}"
ENV_PASA="${ENV_PASA:-${ANACONDA}/annot_pasa}"
ENV_MIKADO="${ENV_MIKADO:-${ANACONDA}/annot_mikado}"
ENV_EVM="${ENV_EVM:-${ANACONDA}/annot_evm}"
ENV_EDTA="${ENV_EDTA:-${ANACONDA}/EDTA}"
ENV_NANOPORE="${ENV_NANOPORE:-${ANACONDA}/nanopore}"

export PASAHOME="${PASAHOME:-${ENV_PASA}/opt/pasa-2.5.3}"
export EVM_HOME="${EVM_HOME:-${ENV_EVM}/opt/evidencemodeler-2.1.0}"

PROT_DB="${PROT_DB:-/tiandata2/zzh/sjy_assembly/anno/Viridiplantae.fa}"
PROT_DMND="${PROT_DMND:-/tiandata2/zzh/sjy_assembly/anno/Viridiplantae.dmnd}"
MIKADO_SCORING="${MIKADO_SCORING:-${ENV_MIKADO}/lib/python3.11/site-packages/Mikado/configuration/scoring_files/plant.yaml}"

STRINGTIE="${STRINGTIE:-${ENV_ALIGNER}/bin/stringtie}"
GFFREAD="${GFFREAD:-${ENV_ANNP}/bin/gffread}"
DIAMOND="${DIAMOND:-${ENV_ANNP}/bin/diamond}"
MINIPROT="${MINIPROT:-${ENV_ALIGNER}/bin/miniprot}"
SAMTOOLS="${SAMTOOLS:-${ENV_ALIGNER}/bin/samtools}"
SEQKIT="${SEQKIT:-${ENV_NANOPORE}/bin/seqkit}"
CDHIT="${CDHIT:-${ENV_EDTA}/bin/cd-hit-est}"
GTF2GFF="${GTF2GFF:-${ENV_BRAKER}/bin/gtf2gff.pl}"

timestamp() { date "+%F %T"; }
log()       { echo "[$(timestamp)] $*"; }
die()       { echo "[$(timestamp)] ERROR: $*" >&2; exit 1; }

# 把指定 env 的 bin 插到 PATH 最前（可叠加多个）
use_env() {
    local name="$1"
    local bin="${ANACONDA}/${name}/bin"
    [[ -d "$bin" ]] || die "找不到 conda env bin：$bin（环境名 $name）"
    case ":$PATH:" in
        *":${bin}:"*) ;;
        *) export PATH="${bin}:${PATH}" ;;
    esac
    hash -r
}

require_bin() {
    local exe="$1"
    local hint="${2:-}"
    command -v "$exe" >/dev/null 2>&1 || die "找不到命令 $exe${hint:+  ($hint)}"
}

# Trinity 2.15.2 的 salmon_runner.pl 仍用 1.x 参数 --minAssignedFrags（salmon 2 已删除）。
# 只检查 PATH 上的 salmon；09 polish 走 annp，不受 annot_trinity wrapper 影响。
require_salmon_1x() {
    require_bin salmon "annot_trinity env"
    local ver
    ver=$(salmon --version 2>&1 | head -1)
    log "salmon: $ver  ($(command -v salmon))"
    if ! [[ "$ver" =~ ^salmon[[:space:]]+1\. ]]; then
        die "Trinity 2.15.2 需要 salmon 1.x（--minAssignedFrags），当前是：$ver"
    fi
}

# Trinity --full_cleanup 会把 ${tr_dir}/Trinity-GG.fasta
# 挪到 ${tr_dir}.Trinity-GG.fasta，工作目录里不再留 fasta。
link_trinity_gg_fasta() {
    local tr_dir="$1"
    local tr_fa="${tr_dir}/Trinity-GG.fasta"
    local cleaned="${tr_dir}.Trinity-GG.fasta"
    [[ -s "$tr_fa" ]] && return 0
    if [[ -s "$cleaned" ]]; then
        mkdir -p "$tr_dir"
        ln -sfn "$cleaned" "$tr_fa"
        log "链接 Trinity --full_cleanup 产物：$tr_fa -> $cleaned"
        return 0
    fi
    return 1
}

# portcullis full 会做 sklearn 自训练；annot_mikado 没有 sklearn，C++ 报
# "basic_string: construction from null is not valid"。拆成 prep/junc/filter --no_ml。
# 已有 2-junc/*.junctions.tab 则跳过 ident，只做规则过滤。
run_portcullis_no_ml() {
    local genome="$1" bam="$2" pc_dir="$3" threads="$4" max_intron="$5" logf="${6:-/dev/stderr}"
    local prep="${pc_dir}/1-prep"
    local tab=""
    mkdir -p "$pc_dir"
    tab=$(find "$pc_dir" -name '*.junctions.tab' 2>/dev/null | head -1 || true)
    if [[ -z "${tab:-}" || ! -s "$tab" ]]; then
        log "Portcullis prep -> $prep"
        portcullis prep -t "$threads" -o "$prep" "$genome" "$bam" >>"$logf" 2>&1
        log "Portcullis junc"
        portcullis junc -t "$threads" --orientation UNKNOWN --strandedness unstranded \
            -o "${pc_dir}/2-junc/portcullis" "$prep" >>"$logf" 2>&1
        tab=$(find "$pc_dir" -name '*.junctions.tab' 2>/dev/null | head -1 || true)
    else
        log "已有 Portcullis junctions tab，跳过 ident：$tab"
    fi
    [[ -s "${tab:-}" ]] || return 1
    [[ -d "$prep" ]] || return 1
    log "Portcullis filter --no_ml (canonical C,S min_cov=2 max_intron=$max_intron)"
    mkdir -p "${pc_dir}/3-filt"
    portcullis filter --no_ml --canonical C,S --min_cov 2 --max_length "$max_intron" \
        -t "$threads" -o "${pc_dir}/3-filt/portcullis" "$prep" "$tab" >>"$logf" 2>&1
}

resolve_genome() {
    local sample="$1"
    local sp="${sample%%_*}"
    local g="${ROOT}/07_polish_rename/${sample}/03_final/${sample}_Chr.fa"
    if [[ ! -s "$g" ]]; then
        g="${SCRIPT_DIR}/${sp}_TE/${sp}_Chr.fa.masked"
    fi
    [[ -s "$g" ]] || die "缺少基因组：$g"
    echo "$g"
}

resolve_bam() {
    local sp="$1"
    local bam="${SCRIPT_DIR}/${sp}_RNA/${sp}_RNA_sorted.bam"
    echo "$bam"
}
