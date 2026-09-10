#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# 05  Portcullis + TransDecoder + Mikado
#     从 StringTie/Trinity 里挑可靠转录本（给 EVM 当 OTHER_PREDICTION）
#
# 环境：
#   annot_mikado        portcullis, mikado
#   annot_gmap          gmap / gmap_build（Trinity fasta → GFF3）
#   annot_transdecoder  TransDecoder.LongOrfs / Predict
#   annp                diamond
#
# 输入：04 的 stringtie.gtf、Trinity-GG.fasta、BAM、基因组
# 输出：
#   ${sp}_models/05_mikado/${sp}.mikado.loci.gff3
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
MAX_INTRON="${MAX_INTRON:-20000}"

[[ -s "$MIKADO_SCORING" ]] || die "找不到 Mikado plant.yaml：$MIKADO_SCORING"
[[ -s "$PROT_DB" ]] || die "缺少蛋白库：$PROT_DB"

log "THREADS=$THREADS scoring=$MIKADO_SCORING"

for SAMPLE in "${SAMPLES[@]}"; do
    [[ -z "$FILTER_SAMPLE" || "$SAMPLE" == "$FILTER_SAMPLE" ]] || continue
    sp="${SAMPLE%%_*}"
    genome="$(resolve_genome "$SAMPLE")"
    bam="$(resolve_bam "$sp")"
    d04="${SCRIPT_DIR}/${sp}_models/04_transcripts"
    st_gtf="${d04}/${sp}.stringtie.gtf"
    link_trinity_gg_fasta "${d04}/trinity_gg" || true
    tr_fa="${d04}/trinity_gg/Trinity-GG.fasta"
    outdir="${SCRIPT_DIR}/${sp}_models/05_mikado"
    loci="${outdir}/${sp}.mikado.loci.gff3"

    log "============================================================"
    log "样本：$SAMPLE"
    log "============================================================"
    [[ -s "$st_gtf" || -s "$tr_fa" ]] || { log "跳过 $SAMPLE：先跑 04"; continue; }
    [[ -s "$bam" ]] || { log "跳过 $SAMPLE：缺少 BAM"; continue; }

    if [[ -s "$loci" && "$FORCE_RERUN" != "1" ]]; then
        log "已有 Mikado loci，跳过：$loci"
        continue
    fi
    mkdir -p "$outdir"

    # --- Portcullis ---
    pc_dir="${outdir}/portcullis"
    junc=""
    use_env annot_mikado
    require_bin portcullis "annot_mikado"
    if [[ "$FORCE_RERUN" != "1" ]]; then
        junc=$(find "$pc_dir" -name '*pass.junctions.bed' 2>/dev/null | head -1 || true)
    fi
    if [[ -z "${junc:-}" || ! -s "${junc:-}" ]]; then
        log "Portcullis prep/junc/filter --no_ml"
        run_portcullis_no_ml "$genome" "$bam" "$pc_dir" "$THREADS" "$MAX_INTRON" \
            "${outdir}/${sp}.portcullis.log" \
            || die "Portcullis 失败：见 ${outdir}/${sp}.portcullis.log"
        junc=$(find "$pc_dir" -name '*pass.junctions.bed' | head -1 || true)
    fi
    [[ -s "${junc:-}" ]] || die "Portcullis 未产出 junctions bed：见 ${outdir}/${sp}.portcullis.log"

    # --- GMAP Trinity → GFF3 ---
    tr_gff="${outdir}/${sp}.trinity.gmap.gff3"
    if [[ -s "$tr_fa" ]]; then
        if [[ ! -s "$tr_gff" || "$FORCE_RERUN" == "1" ]]; then
            use_env annot_gmap
            require_bin gmap_build "annot_gmap"
            require_bin gmap "annot_gmap"
            idx="${outdir}/gmap_idx"
            mkdir -p "$idx"
            if [[ ! -d "${idx}/${sp}" ]]; then
                log "gmap_build $sp"
                gmap_build -D "$idx" -d "$sp" "$genome" >"${outdir}/${sp}.gmap_build.log" 2>&1
            fi
            log "GMAP Trinity -> $tr_gff"
            gmap -D "$idx" -d "$sp" -f 2 -n 1 -t "$THREADS" --nofails \
                "$tr_fa" >"$tr_gff" 2>"${outdir}/${sp}.gmap.log"
        fi
    fi

    # --- Mikado list ---
    list="${outdir}/mikado.list.tsv"
    : >"$list"
    if [[ -s "$st_gtf" ]]; then
        printf '%s\tStringTie\tFalse\t1\tFalse\tTrue\tFalse\tFalse\n' "$st_gtf" >>"$list"
    fi
    if [[ -s "$tr_gff" ]]; then
        printf '%s\tTrinity\tFalse\t1\tFalse\tTrue\tFalse\tFalse\n' "$tr_gff" >>"$list"
    fi
    [[ -s "$list" ]] || die "Mikado list 为空"

    use_env annot_mikado
    require_bin mikado "annot_mikado"
    cfg="${outdir}/configuration.yaml"
    prep_fa="${outdir}/mikado_prepared.fasta"
    if [[ ! -s "$cfg" || ! -s "$prep_fa" || "$FORCE_RERUN" == "1" ]]; then
        log "mikado configure"
        mikado configure --list "$list" --reference "$genome" --junctions "$junc" \
            --scoring "$MIKADO_SCORING" -bt "$PROT_DB" --mode permissive \
            --max-intron-length "$MAX_INTRON" -t "$THREADS" -od "$outdir" "$cfg" \
            >"${outdir}/${sp}.mikado_configure.log" 2>&1

        log "mikado prepare"
        mikado prepare --configuration "$cfg" --fasta "$genome" -od "$outdir" \
            >"${outdir}/${sp}.mikado_prepare.log" 2>&1
    else
        log "已有 mikado_prepared.fasta，跳过 configure/prepare：$prep_fa"
    fi
    prep_fa=$(find "$outdir" -name 'mikado_prepared.fasta' | head -1 || true)
    [[ -s "${prep_fa:-}" ]] || die "mikado prepare 未产出 fasta"

    # --- TransDecoder on prepared transcripts ---
    use_env annot_transdecoder
    require_bin TransDecoder.LongOrfs "annot_transdecoder"
    td_dir="${outdir}/transdecoder"
    mkdir -p "$td_dir"
    ln -sfn "$prep_fa" "${td_dir}/mikado_prepared.fasta"
    if [[ ! -s "${td_dir}/mikado_prepared.fasta.transdecoder.bed" || "$FORCE_RERUN" == "1" ]]; then
        log "TransDecoder"
        (
            cd "$td_dir"
            TransDecoder.LongOrfs -t mikado_prepared.fasta
            TransDecoder.Predict -t mikado_prepared.fasta
        ) >"${outdir}/${sp}.transdecoder.log" 2>&1
    fi
    orf_bed=$(find "$td_dir" -name '*.transdecoder.bed' | head -1 || true)
    orf_pep=$(find "$td_dir" -name '*.transdecoder.pep' | head -1 || true)
    [[ -s "${orf_bed:-}" ]] || die "TransDecoder 未产出 bed：见 ${outdir}/${sp}.transdecoder.log"

    # --- DIAMOND blastp ---
    use_env annp
    require_bin diamond "annp"
    if [[ ! -s "$PROT_DMND" ]]; then
        log "diamond makedb $PROT_DB"
        "$DIAMOND" makedb --in "$PROT_DB" --db "$PROT_DMND" --threads "$THREADS"
    fi
    blast_xml="${outdir}/diamond.xml"
    if [[ ! -s "$blast_xml" || "$FORCE_RERUN" == "1" ]]; then
        log "diamond blastp -> $blast_xml"
        "$DIAMOND" blastp --db "$PROT_DMND" --query "$orf_pep" --out "$blast_xml" \
            --outfmt 5 --max-target-seqs 5 --evalue 1e-5 --threads "$THREADS" \
            >"${outdir}/${sp}.diamond.log" 2>&1
    fi

    # --- serialise + pick ---
    use_env annot_mikado
    log "mikado serialise"
    mikado serialise --configuration "$cfg" --transcripts "$prep_fa" \
        --orfs "$orf_bed" --xml "$blast_xml" --genome "$genome" \
        --junctions "$junc" -p "$THREADS" --force -od "$outdir" \
        >"${outdir}/${sp}.mikado_serialise.log" 2>&1

    log "mikado pick"
    mikado pick --configuration "$cfg" --fasta "$genome" -p "$THREADS" \
        --source Mikado --loci-out "$(basename "$loci")" -od "$outdir" \
        >"${outdir}/${sp}.mikado_pick.log" 2>&1

    # pick 可能把 loci 写在 outdir 下固定名
    if [[ ! -s "$loci" ]]; then
        alt=$(find "$outdir" -name '*loci.gff*' | head -1 || true)
        [[ -n "${alt:-}" && -s "$alt" ]] && cp "$alt" "$loci"
    fi
    [[ -s "$loci" ]] || die "Mikado pick 未产出 loci GFF3"
    log "完成 $SAMPLE  $loci"
done

log "========== 05 Mikado 结束 =========="
log "下一步：bash 06_run_pasa_align.sh"
