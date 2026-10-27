#!/usr/bin/env bash
# Write the gene functional table for each 13_final set made by 10_final_anno.sh.
#
# final_function_table.py reads 13_final and does not change the GFF or FASTA.
# A final transcript whose protein is identical to a step-09 protein reuses that
# protein's eggNOG, InterProScan, SwissProt, ZH13 v3 and Wm82 a6 results. A
# protein changed in step 10 (extension, merge) is annotated again with the same
# tools, databases and settings as step 09. Those results go to
#   05_diamond/PREFIX.final_novel.*   06_eggnog/PREFIX_final_novel.*
#   07_interpro/PREFIX.final_novel.interproscan.*
# and are reused on a later run when the novel protein FASTA is unchanged.
# TE overlap is recomputed for genes whose coordinates changed. Expression stays
# the step-09 quantification of that gene ID.
#
# Output: 13_final/PREFIX.final.functional_annotation.tsv, one row per gene.
#
# Usage:
#   bash 11_final_function.sh
#   SAMPLE=KD35 bash 11_final_function.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${ROOT:-/tiandata2/zzh/sjy_assembly}"
FILTER_SAMPLE="${SAMPLE:-}"
FUNCTION_THREADS="${FUNCTION_THREADS:-16}"

DOWNSTREAM_DIRS=(
    "${ROOT}/annot/KD35_annotation_downstream"
    "${ROOT}/annot/N2L_annotation_downstream"
    "${ROOT}/batch2/annot/BZ_annotation_downstream"
    "${ROOT}/batch2/annot/DY11_annotation_downstream"
    "${ROOT}/batch2/annot/TH2_annotation_downstream"
    "${ROOT}/batch2/annot/TZX1998_annotation_downstream"
    "${ROOT}/batch2/annot/TZX253_annotation_downstream"
)

timestamp() { date "+%F %T"; }
log() { echo "[$(timestamp)] $*"; }
die() { echo "[$(timestamp)] ERROR: $*" >&2; exit 1; }

downs=()
prefixes=()
for down in "${DOWNSTREAM_DIRS[@]}"; do
    base="$(basename "$down")"
    prefix="${base%%_annotation_downstream}"
    if [[ -n "$FILTER_SAMPLE" && "$prefix" != "$FILTER_SAMPLE" ]]; then
        continue
    fi
    [[ -s "${down}/13_final/${prefix}.final.gff3" ]] || die "missing 13_final for ${prefix}; run 10_final_anno.sh first"
    downs+=("$down")
    prefixes+=("$prefix")
done
[[ ${#downs[@]} -gt 0 ]] || die "no sample matched SAMPLE=${FILTER_SAMPLE}"

# Novel proteins of all selected samples go through eggNOG and InterProScan in
# one run, so the eggNOG database is scanned once.
log "functional tables: ${prefixes[*]}"
python3 "${SCRIPT_DIR}/final_function_table.py" \
    --downstream "${downs[@]}" \
    --prefix "${prefixes[@]}" \
    --threads "$FUNCTION_THREADS" \
    --work "${ROOT}/final_function_work"

log "done"
