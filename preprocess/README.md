# Preprocess helpers

Optional, site-specific cleaning. Edit barcodes and paths before use.

| Script | Role |
|--------|------|
| `hifi_bam2fastq.sh` | HiFi BAM → filtered FASTQ/FASTA |
| `hic_fastp.sh` | Hi-C adapter/quality trim (`fastp`) |
| `rna_salt_clean.sh` | RNA-seq cleaning prior to Atlas merge into `rna_clean/` |

Expected cleaned names consumed by the main pipeline:

- `hifi/*_hifi.filt.fastq.gz`
- `ont/<cultivar>_ONT_clean.fq.gz`
- `hic/<cultivar>_clean_{1,2}.fq.gz`
- `rna_clean/<cultivar>_RNA_Atlas_{1,2}.fq.gz`
