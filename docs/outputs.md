# Outputs for a paper

中文见文末。

After assembly 07 and annotation 09, the files that belong in a manuscript / database submission are:

## Genome

| File | Meaning |
|------|---------|
| `07_polish_rename/<sample>/03_final/<sample>_Chr.fa` | Chromosome-level FASTA (Chr01–Chr20 + unplaced if present) |
| `08_fa_qc/` | Merqury QV, BUSCO, QUAST, synteny plots |
| `09_run_hic_heatmap/` | Hi-C cooler + heatmap PDF/PNG |
| `10_run_centro/<sample>/` | Centromere calls |

## Gene annotation (use `12_renamed/`)

Directory: `<workdir>/annot/<cultivar>_annotation_downstream/12_renamed/`  
(or the equivalent path if you keep annotation outputs under `ROOT`).

### Main paper set (HC)

| File | Use |
|------|-----|
| `*.renamed.high_confidence.annotated.gff3` | Submit / browse; structural gene models |
| `*.renamed.high_confidence.proteins.fa` | Protein set matching HC genes |
| `*.renamed.high_confidence.transcripts.fa` | mRNA / transcripts |
| `*.renamed.high_confidence.cds.fa` | CDS |
| `*.renamed.high_confidence.genes.list` | One gene ID per line |

Do **not** report gene counts from:

- full `*.renamed.gff3` (all scaffolds, all loci)
- `*.renamed.nuclear.gff3` (Chr01–20, protein ≥ 29 aa, **no** HC quality cut)
- `11_publication/*.high_confidence.*` if generated with an older filter (all chromosomes, no ZH13 identity gate)

### Supplement: lineage-specific

| File | Use |
|------|-----|
| `*.renamed.lineage_specific.genes.list` | Gene IDs |
| `*.renamed.lineage_specific.annotated.gff3` | GFF with Name/product |
| `*.renamed.lineage_specific.transcripts.fa` | mRNA |
| `*.renamed.lineage_specific.cds.fa` | CDS |
| `*.renamed.lineage_specific.proteins.fa` | Proteins |
| `*.renamed.lineage_specific.README.txt` | Filter text and counts |

### Tables

| File | Use |
|------|-----|
| `*.renamed.annotation.summary.tsv` | Per-gene homologs, eggNOG, InterPro, TE, TPM, confidence |
| `*.renamed.gene_functional_annotation.tsv` | Compact functional table |
| `*.gene_id_mapping.tsv` | Old ID → ZH13-style ID |

CDS/mRNA for the **HC** set are written as `*.renamed.high_confidence.transcripts.fa` and `*.renamed.high_confidence.cds.fa`.

## Example counts (KD35 / N2L, this project)

These numbers document one successful run of the HC definition above; they are not a specification.

| Set | KD35 | N2L |
|-----|------|-----|
| HC genes | 59,284 | 56,547 |
| Lineage-specific (expressed, ZH13 < 50) | 1,916 | 1,758 |
| Protein BUSCO embryophyta_odb10 (all predicted proteins) | C:99.3% | C:99.0% |

Published soybean nuclear gene counts typically fall near 5.4–6.0 × 10⁴; ZH13 v3 nuclear is 53,401. HC is intended to sit in that range.

---

## 发表用文件（中文）

基因组用 07 的 `*_Chr.fa`。基因数和投稿注释只用 `12_renamed` 的 **HC 三件套**。全量 / nuclear GFF 是备份和浏览器轨道，不是论文基因数。品种特异另附 `lineage_specific`。功能表用 `*.renamed.annotation.summary.tsv`。
