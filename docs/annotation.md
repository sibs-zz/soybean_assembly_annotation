# Gene annotation

中文见文末。 / English first; Chinese follows.

Annotation is evidence-based (no Iso-Seq in the default KD35/N2L run). Protein homology uses a Viridiplantae set plus five soybean reference proteomes.

## Order

| Step | Script | Role |
|------|--------|------|
| 01 | `01_EDTA.sh` | EDTA (apptainer) + RepeatMasker soft-mask; TE library and `.out` |
| 02a | `02_LAI.sh` | LTR Assembly Index from EDTA intact LTRs + RM `.out` |
| 02b | `02_run_HISAT2_align.sh` | HISAT2 on masked genome → `${sp}_RNA/${sp}_RNA_sorted.bam` |
| 03 | `03_run_BRAKER3.sh` | BRAKER3 (GeneMark-ETP + ProtHint + RNA BAM + proteins) |
| 04 | `04_run_transcripts.sh` | StringTie + genome-guided Trinity; merge for PASA |
| 05 | `05_run_mikado.sh` | Portcullis (no ML), GMAP Trinity, TransDecoder, DIAMOND, Mikado |
| 06 | `06_run_pasa_align.sh` | PASA SQLite align/assemble (**do not** pass `--trans_gtf`) |
| 07 | `07_run_evm.sh` | miniprot (Viridiplantae) + EVM; weights in `evm.weights.txt` |
| 08 | `08_run_pasa_update.sh` | PASA annot compare `-A -L --annots evm.gff3` (UTRs) |
| 09 | `09_run_polish.sh` | Sequences, homologs, eggNOG, InterPro, TE overlap, Salmon, HC filter, ZH13-style rename, lineage-specific extract |

Shared settings: `annot_common.sh` (sourced by 04–08). Python helpers live in the same folder as the shell scripts.

Typical launch for one sample:

```bash
SAMPLE=KD35_primary THREADS=16 bash 04_run_transcripts.sh
```

`09_run_polish.sh` defaults to 64 threads; InterProScan and eggNOG are the long steps. eggNOG mapper must see Diamond **≥ 2.2** (`--iterate`). Trinity 2.15.2 in this pipeline requires **salmon 1.x**.

## Confidence labels (09)

Assigned per gene from expression (Salmon TPM) and functional hits (SwissProt / soybean refs / eggNOG / InterPro):

- TE overlap ≥ 0.5 **and** no expression **and** no function → `low_confidence_possible_TE`
- expression **and** function → `high_confidence`
- exactly one of expression or function → `medium_confidence`
- neither → `low_confidence`

## High-confidence (HC) vs lineage-specific

HC (main paper set), all required:

1. `Chr01`–`Chr20`
2. confidence ≥ medium
3. TE overlap ratio ≤ 0.5
4. best hit vs ZH13 v3 nuclear proteins identity ≥ 50

Lineage-specific (supplement), all required:

1. `Chr01`–`Chr20`
2. confidence ≥ medium
3. TE overlap ratio ≤ 0.5
4. `expression_supported=yes`
5. ZH13 v3 identity **< 50** (no hit counts as 0)

The two sets are disjoint. IDs are `PREFIX_chrG` + 6 digits, step 100 (e.g. `KD35_01G000100.m1`). Filtering does **not** re-number, so HC IDs may have gaps.

## EVM weights

`annotation/evm.weights.txt` (copy into the EVM workdir if the runner does not already write it):

```
ABINITIO_PREDICTION	AUGUSTUS	1
ABINITIO_PREDICTION	GeneMark	1
OTHER_PREDICTION	BRAKER	4
OTHER_PREDICTION	Mikado	7
TRANSCRIPT	PASA	10
PROTEIN	miniprot	8
```

PASA transcripts outrank protein and ab initio tracks.

## Known operational constraints

- Portcullis is run as prep/junc/filter `--no_ml` (sklearn not required).
- Trinity `--full_cleanup` moves FASTA to `${dir}.Trinity-GG.fasta`; `link_trinity_gg_fasta()` restores the expected path.
- PASA 08 must use `-A -L` with EVM GFF; do not set `GENETIC_CODE=universal` in this soybean run.
- Step 09 PATH should put `eggnog213` (Python) and `annp` Diamond 2.x ahead of any Diamond 0.9.x.

---

## 注释说明（中文）

默认流程对标 ZH13 v3：短读 RNA + 绿色植物蛋白 + 五套大豆参考蛋白，不用 Iso-seq。01 EDTA 必须用容器；03 BRAKER3 需要 GeneMark-ETP 许可。04–08 共用 `annot_common.sh`。EVM 权重见上，PASA 转录本权重最高。

09 步才做功能注释、置信度、HC 过滤和 ZH13 风格改名。论文主表用 HC 三件套；品种特异用 `lineage_specific` 一套。全量 nuclear GFF 仍保留未进 HC 的基因。
