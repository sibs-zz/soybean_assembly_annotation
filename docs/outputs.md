# Outputs for a paper

中文见文末。

After assembly 07 and annotation 10, the files that belong in a manuscript / database submission are:

## Genome

| File | Meaning |
|------|---------|
| `07_polish_rename/<sample>/03_final/<sample>_Chr.fa` | Chromosome-level FASTA (Chr01–Chr20 + unplaced if present) |
| `08_fa_qc/` | Merqury QV, BUSCO, QUAST, synteny plots |
| `09_run_hic_heatmap/` | Hi-C cooler + heatmap PDF/PNG |
| `10_run_centro/<sample>/` | Centromere calls |

## Gene annotation (use `13_final/`)

Directory: `<workdir>/annot/<cultivar>_annotation_downstream/13_final/`  
(or `batch2/annot/<cultivar>_annotation_downstream/13_final/`).

`12_renamed/` is the step-09 renamed annotation and drop list. Step 10 reads it and writes the files below. Repeat features stay in `12_renamed/*.renamed.repeat.gff3.gz`; TE IDs are not gene loci.

### Submitted set

| File | Use |
|------|-----|
| `*.final.gff3` | One structural gene set |
| `*.final.transcripts.fa` | mRNA |
| `*.final.cds.fa` | CDS, including the terminal stop |
| `*.final.proteins.fa` | Proteins, without the terminal stop |
| `*.final.genes.list` | One gene ID per line, same order as the GFF |
| `*.final.gene_class.tsv` | `high_confidence`, `expressed_zh13_below_50`, or `homolog_unexpressed` |
| `*.final.criteria.txt` | Filter and repair rules used for this run |
| `*.final.corrections.tsv` | Short-model review; old and new protein when a model was overwritten |
| `*.final.id_map.tsv` | Old gene or transcript ID → kept ID |
| `*.final.redundant_merged.tsv` | Duplicate CDS/protein gene ID removed in step 10 |
| `*.final.extension_conflict.tsv` | Extension restored because the new CDS overlapped another extension |
| `*.final.nested_fragment.tsv` | Same-strand fragment inside a longer gene of the same ZH13 locus |
| `*.final.antisense_pair_qc.tsv` | Every opposite-strand CDS overlap and its evidence tier, written before exclusion |
| `*.final.antisense_model_qc.tsv` | Single-exon antisense model lacking reliable homology |
| `*.final.antisense_excluded.tsv` | Gene removed from the final set by the antisense filter |
| `*.final.orf_dropped.tsv` | ORF models removed while building this set |
| `*.final.busco.short_summary.txt` | BUSCO protein mode, longest protein per gene, embryophyta_odb10 |
| `*.final.busco.full_table.tsv` | One row per BUSCO gene |
| `*.final.busco.missing_busco_list.tsv` | Missing BUSCO IDs |

Do **not** report gene counts from:

- `12_renamed/*.renamed.gff3` (all scaffolds, all loci)
- `12_renamed/*.renamed.nuclear.gff3`
- `12_renamed/*.renamed.high_confidence.*` (step-09 subset, before the step-10 merge and repair)
- `11_publication/` if it still uses an older filter

Functional tables from step 09 remain in `12_renamed/*.renamed.annotation.summary.tsv` and `*.gene_id_mapping.tsv`. They describe the renamed models, not the coordinates changed in step 10.

## Example counts (this project, after step 10)

These numbers document one run. They are not a specification. BUSCO is the longest protein of each final gene, measured before the same-strand duplicate removal and before the antisense exclusion.

| Cultivar | Genes | high_confidence | expressed_zh13_below_50 | homolog_unexpressed | BUSCO complete |
|----------|------:|----------------:|------------------------:|--------------------:|----------------|
| KD35 | 61,657 | 58,902 | 932 | 1,823 | C:97.2% |
| N2L | 58,572 | 56,230 | 800 | 1,542 | C:96.9% |
| BZ | 59,931 | 57,254 | 983 | 1,694 | C:96.8% |
| DY11 | 58,970 | 56,380 | 931 | 1,659 | C:97.2% |
| TH2 | 59,493 | 56,895 | 971 | 1,627 | C:97.5% |
| TZX1998 | 59,222 | 56,739 | 927 | 1,556 | C:96.6% |
| TZX253 | 61,642 | 59,058 | 876 | 1,708 | C:96.4% |

---

## 发表用文件（中文）

基因组用 07 的 `*_Chr.fa`。基因数和投稿注释用 `13_final` 的一套 GFF 与 FASTA。`12_renamed` 是改名和删除记录，其中的 high-confidence 子集是 10 步之前的中间结果。功能表仍在 `12_renamed/*.renamed.annotation.summary.tsv`。
