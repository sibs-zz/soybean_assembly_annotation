#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# 07  EvidenceModeler：ab initio + 蛋白 + 转录本 → 共识基因
#
# 环境：annot_evm   export EVM_HOME=.../opt/evidencemodeler-2.1.0
#       aligner     miniprot
#       braker3     gtf2gff.pl
#
# 证据：
#   ABINITIO   AUGUSTUS / GeneMark（来自 03 BRAKER 工作目录）
#   OTHER      BRAKER TSEBRA、Mikado loci
#   TRANSCRIPT PASA assemblies（06）
#   PROTEIN    miniprot(Viridiplantae vs genome)
#
# 输入依赖：03 braker 预测 + 05 mikado（可选）+ 06 pasa
# 输出：
#   ${sp}_models/07_evm/${sp}.evm.gff3
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SAMPLES=(
    "KD35_primary"
    "N2L_primary"
)
ROOT="/tiandata2/zzh/sjy_assembly"
# shellcheck source=/dev/null
source "${SCRIPT_DIR}/annot_common.sh"

THREADS="${THREADS:-16}"
FORCE_RERUN="${FORCE_RERUN:-0}"
FILTER_SAMPLE="${SAMPLE:-}"
SEGMENT_SIZE="${SEGMENT_SIZE:-500000}"
OVERLAP_SIZE="${OVERLAP_SIZE:-50000}"

[[ -d "$EVM_HOME" ]] || die "EVM_HOME 不存在：$EVM_HOME"
[[ -x "$MINIPROT" ]] || die "找不到 miniprot：$MINIPROT"
[[ -s "$PROT_DB" ]] || die "缺少蛋白库：$PROT_DB"
[[ -f "${SCRIPT_DIR}/prep_evm_evidence.py" ]] || die "缺少 prep_evm_evidence.py"

log "EVM_HOME=$EVM_HOME THREADS=$THREADS"

find_first() {
    local f
    for f in "$@"; do
        [[ -s "$f" ]] && { echo "$f"; return 0; }
    done
    return 1
}

for SAMPLE in "${SAMPLES[@]}"; do
    [[ -z "$FILTER_SAMPLE" || "$SAMPLE" == "$FILTER_SAMPLE" ]] || continue
    sp="${SAMPLE%%_*}"
    genome="$(resolve_genome "$SAMPLE")"
    braker="${SCRIPT_DIR}/BRAKER_${sp}"
    outdir="${SCRIPT_DIR}/${sp}_models/07_evm"
    evm_gff="${outdir}/${sp}.evm.gff3"

    log "============================================================"
    log "样本：$SAMPLE"
    log "============================================================"

    aug=$(find_first \
        "${braker}/augustus.hints.gtf" \
        "${braker}/Augustus/augustus.hints.gtf" \
        "${braker}/augustus.gtf" || true)
    gm=$(find_first \
        "${braker}/GeneMark-ETP/genemark.gtf" \
        "${braker}/genemark.gtf" \
        "${braker}/GeneMark-ETP/data/genemark.gtf" || true)
    brk=$(find_first "${braker}/braker.gff3" "${braker}/braker.gtf" || true)
    pasa_gff=$(find_first \
        "${SCRIPT_DIR}/${sp}_models/06_pasa/${sp}.pasa_assemblies.gff3" || true)
    mikado=$(find_first \
        "${SCRIPT_DIR}/${sp}_models/05_mikado/${sp}.mikado.loci.gff3" || true)

    [[ -s "${brk:-}" ]] || { log "跳过 $SAMPLE：缺少 BRAKER GFF（等 03）"; continue; }
    [[ -s "${pasa_gff:-}" ]] || { log "跳过 $SAMPLE：缺少 06 PASA GFF"; continue; }

    if [[ -s "$evm_gff" && "$FORCE_RERUN" != "1" ]]; then
        log "已有 EVM GFF，跳过：$evm_gff"
        continue
    fi
    mkdir -p "$outdir"

    use_env braker3
    gene_aug="${outdir}/augustus.gff3"
    gene_gm="${outdir}/genemark.gff3"
    gene_brk="${outdir}/braker.gff3"
    if [[ -n "${aug:-}" ]]; then
        perl "$GTF2GFF" --gff3 --printExon --out="$gene_aug" <"$aug" || gffread "$aug" -o "$gene_aug"
    fi
    if [[ -n "${gm:-}" ]]; then
        perl "$GTF2GFF" --gff3 --printExon --out="$gene_gm" <"$gm" || gffread "$gm" -o "$gene_gm"
    fi
    if [[ "$brk" == *.gtf ]]; then
        perl "$GTF2GFF" --gff3 --printExon --out="$gene_brk" <"$brk" || gffread "$brk" -o "$gene_brk"
    else
        cp "$brk" "$gene_brk"
    fi

    # miniprot 蛋白比对
    mp_gff="${outdir}/miniprot.gff3"
    if [[ ! -s "$mp_gff" || "$FORCE_RERUN" == "1" ]]; then
        use_env aligner
        log "miniprot $PROT_DB -> $mp_gff"
        "$MINIPROT" -t "$THREADS" --gff "$genome" "$PROT_DB" \
            >"$mp_gff" 2>"${outdir}/${sp}.miniprot.log"
        [[ -s "$mp_gff" ]] || die "miniprot 失败：见 ${outdir}/${sp}.miniprot.log"
    fi

    gene_all="${outdir}/gene_predictions.gff3"
    prot_all="${outdir}/protein_alignments.gff3"
    tx_all="${outdir}/transcript_alignments.gff3"
    evm_py=(
        python3 "${SCRIPT_DIR}/prep_evm_evidence.py"
        --gene-out "$gene_all"
        --protein miniprot "$mp_gff" --protein-out "$prot_all"
        --transcript PASA "$pasa_gff" --transcript-out "$tx_all"
        --gene BRAKER "$gene_brk"
    )
    [[ -s "$gene_aug" ]] && evm_py+=(--gene AUGUSTUS "$gene_aug")
    [[ -s "$gene_gm" ]] && evm_py+=(--gene GeneMark "$gene_gm")
    [[ -n "${mikado:-}" && -s "$mikado" ]] && evm_py+=(--gene Mikado "$mikado")
    "${evm_py[@]}"

    [[ -s "$gene_all" ]] || die "gene_predictions 为空"
    [[ -s "$tx_all" ]] || die "transcript_alignments 为空"

    weights="${outdir}/weights.txt"
    : >"$weights"
    [[ -s "$gene_aug" ]] && echo -e "ABINITIO_PREDICTION\tAUGUSTUS\t1" >>"$weights"
    [[ -s "$gene_gm" ]] && echo -e "ABINITIO_PREDICTION\tGeneMark\t1" >>"$weights"
    echo -e "OTHER_PREDICTION\tBRAKER\t4" >>"$weights"
    [[ -n "${mikado:-}" && -s "$mikado" ]] && echo -e "OTHER_PREDICTION\tMikado\t7" >>"$weights"
    echo -e "TRANSCRIPT\tPASA\t10" >>"$weights"
    echo -e "PROTEIN\tminiprot\t8" >>"$weights"

    use_env annot_evm
    export EVM_HOME
    export PATH="${ENV_EVM}/bin:${EVM_HOME}/EvmUtils:${PATH}"
    require_bin EVidenceModeler "annot_evm"

    log "EVidenceModeler sample_id=$sp"
    (
        cd "$outdir"
        EVidenceModeler \
            --sample_id "$sp" \
            --genome "$genome" \
            --weights "$weights" \
            --gene_predictions "$gene_all" \
            --protein_alignments "$prot_all" \
            --transcript_alignments "$tx_all" \
            --segmentSize "$SEGMENT_SIZE" \
            --overlapSize "$OVERLAP_SIZE" \
            --CPU "$THREADS" \
            --exec_dir "$outdir"
    ) >"${outdir}/${sp}.evm.log" 2>&1

    # EVM 2.1 输出名通常是 ${sample_id}.EVM.gff3 或 evm.out.gff3
    found=$(find "$outdir" -maxdepth 1 \( -name "${sp}.EVM.gff3" -o -name '*.EVM.gff3' -o -name 'evm.out.gff3' -o -name '*evm*.gff3' \) | head -1 || true)
    [[ -s "${found:-}" ]] || die "EVM 未产出 GFF3：见 ${outdir}/${sp}.evm.log"
    cp "$found" "$evm_gff"
    log "完成 $SAMPLE  $evm_gff"
done

log "========== 07 EVM 结束 =========="
log "下一步：bash 08_run_pasa_update.sh"
