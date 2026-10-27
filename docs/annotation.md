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
| 09 | `09_run_polish.sh` | Sequences, homologs, eggNOG, InterPro, TE overlap, Salmon, confidence labels, ZH13-style rename, drop list |
| 10 | `10_final_anno.sh` | From `12_renamed`: ORF drop, three-class final set, short-model repair, same-strand cleanup, antisense filter, short-fragment filter, BUSCO |
| 11 | `11_final_function.sh` | Functional table for the `13_final` genes; re-annotates proteins changed in step 10 |

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

## Final set (step 10)

`10_final_anno.sh` calls `curate_final_models.py`, `resolve_final_redundancy.py`, `antisense_qc.py`, then `short_fragment_filter.py`. It does not rerun steps 01–09. Input is an existing `12_renamed/` directory. Output is one GFF and one FASTA set in `13_final/`.

Shared gates: `Chr01`–`Chr20`, confidence ≥ medium, TE overlap ≤ 0.5. The three classes are disjoint:

| Class | Extra rule |
|-------|------------|
| `high_confidence` | ZH13 v3 identity ≥ 50; expression not required |
| `expressed_zh13_below_50` | Expression supported and ZH13 v3 identity < 50 (no hit counts as 0). Not rewritten from a full ZH13 protein. |
| `homolog_unexpressed` | No expression, ZH13 v3 identity < 50, and at least one of ZH13 v2, ZH13 v3, Wm82 a2.v1, Wm82 gnm4, Wm82 gnm6 whose own highest-bitscore hit has identity ≥ 50, both coverages ≥ 50%, E-value ≤ 1e-5 |

Before the merge, ORF abnormalities are removed again. A non-multiple-of-3 CDS removes the gene. Missing ATG, an internal stop, or a missing terminal stop removes that transcript. Previously dropped gene IDs stay out. IDs are `PREFIX_chrG` + 6 digits, step 100, and are not renumbered.

Short high-confidence proteins (query coverage ≥ 50% and subject coverage < 50% against ZH13 v3) may be extended or merged from a same-locus miniprot alignment. The alignment must cover ≥ 80% of the ZH13 protein at identity ≥ 80%, with no frameshift and no in-frame stop. The replacement CDS must start with ATG, end with a stop, have length divisible by 3, and have no internal stop. Otherwise the original model is kept. Two extensions whose new same-strand CDS overlap by ≥ 200 bp are both rejected. A merge keeps the lowest locus ID; other IDs go to `*.final.id_map.tsv`.

`resolve_final_redundancy.py` then edits that same directory. Identical CDS and protein sequences on one strand are collapsed to the gene ID that carries the UTR transcripts. A same-strand gene is also removed when every one of its CDS coordinate sets is already a transcript of another gene at the same ZH13 locus; it is not added as an isoform. An extension that created a same-strand CDS overlap of ≥ 200 bp is restored from `12_renamed`. An unexpressed fragment is removed when its CDS and protein sit inside a longer gene of the same ZH13 locus. The three records are `*.final.redundant_merged.tsv`, `*.final.extension_conflict.tsv`, and `*.final.nested_fragment.tsv`.

`antisense_qc.py` keeps an opposite-strand CDS overlap only when both genes have reliable homology (`high_confidence`, or `homolog_unexpressed`). Pairs whose ZH13 v3 homologs also overlap on opposite strands stay, as do pairs whose ZH13 homologs do not. A single-exon model without reliable homology is removed, and so is every other antisense pair outside those two groups. When both genes map to one ZH13 gene, the model with more exons stays. Removed genes are listed in `*.final.antisense_excluded.tsv` and `*.final.id_map.tsv`. Locus IDs are not renumbered.

`short_fragment_filter.py` aligns the longest final protein of each gene to ZH13 v3 and Wm82 a6 with DIAMOND (E-value ≤ 1e-5, five targets), merging HSPs on each reference protein. An unexpressed `high_confidence` gene is a candidate when no reference protein with identity ≥ 50 is covered by half or more, and the highest-bitscore ZH13 v3 gene already has a copy in this genome covered ≥ 80% at identity ≥ 50. Positions of the candidate and of every full-length copy are converted to ZH13 v3 coordinates by interpolating between the nearest anchor genes on each side. An anchor is a gene whose top ZH13 hit has identity ≥ 98 and subject coverage ≥ 0.9 and is not the top hit of any other gene; the anchor of the locus being tested is skipped. Leave-one-out on the anchors gives a median error of 0.3–0.5 kb and a 90th percentile of 2–4 kb, against 0.1–0.25 Mb and 1.1–2.3 Mb for raw genome coordinates (all seven sets). The Chr11–Chr13 translocated segments in N2L, BZ, DY11, TH2 and TZX1998 are therefore compared at their ZH13 locus. A candidate that maps to the chromosome of that ZH13 gene, closer to the ZH13 gene position than every full-length copy, is treated as the same-locus model and kept. Other candidates are removed as partial copies or pseudogene fragments. Candidates and outcomes are in `*.final.short_fragment_qc.tsv`; removed IDs also go to `*.final.id_map.tsv`. Alignments are written to `05_diamond/PREFIX.final_longest.vs_*.tsv`.

```bash
bash 10_final_anno.sh
SAMPLE=KD35 bash 10_final_anno.sh
RUN_BUSCO_FINAL=0 bash 10_final_anno.sh
```

`SAMPLE` here is `KD35`, not `KD35_primary`. BUSCO uses the longest protein of each final gene (`embryophyta_odb10`, protein mode) and copies the summary into `13_final/`.

## Functional table (step 11)

`11_final_function.sh` calls `final_function_table.py` for each `13_final` set. It does not change the GFF or FASTA. A final transcript whose protein is identical to a step-09 protein reuses that protein's eggNOG, InterProScan, SwissProt, ZH13 v3 and Wm82 a6 results. A protein changed in step 10 (extension, merge) is annotated again with the same eggNOG-mapper, InterProScan, SwissProt and reference DIAMOND settings as step 09. Those results are written to `05_diamond/PREFIX.final_novel.*`, `06_eggnog/PREFIX_final_novel.*` and `07_interpro/PREFIX.final_novel.interproscan.*`, and are reused when the novel FASTA is unchanged.

Per gene, eggNOG, SwissProt and reference hits come from the longest protein that has a result; InterPro IDs and GO terms are the union over transcripts. `product` is the eggNOG description, then the SwissProt description, then `hypothetical protein`. TE overlap is recomputed with the step-09 method for genes whose coordinates changed. `max_TPM` and `expression_supported` stay the step-09 quantification of the gene ID. Rerun step 11 after any rerun of step 10.

```bash
bash 11_final_function.sh
SAMPLE=KD35 bash 11_final_function.sh
```

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

09 步做功能注释、置信度、ZH13 风格改名和删除清单，结果留在 `12_renamed/`。10 步读取这一层，再做 ORF 删除、三类基因合并，并在原基因号上修补 high-confidence 短模型，写成 `13_final/` 的一套 GFF 和 FASTA。论文基因数用 `13_final`。`expressed_zh13_below_50` 不按 ZH13 全长蛋白改写。
