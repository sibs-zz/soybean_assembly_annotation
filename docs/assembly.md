# T2T assembly

中文见文末。 / English first; Chinese follows.

## Order

| Step | Script | Input | Output |
|------|--------|--------|--------|
| 01 | `01_run_hifiasm_t2t.sh` | HiFi + ONT ultra-long + Hi-C | `01_assembly_out/<sample>.fasta` (primary contig from `.hic.p_ctg.gfa`) |
| 02 | `02_run_yahs.sh` | Primary FASTA + Hi-C | YaHS scaffolds, BAM with `-5SP` |
| 03 | `03_run_jbat_hic.sh` | YaHS bin/AGP | Juicebox JBAT `.hic` / `.assembly` |
| 04 | `04_t2t.sh` | Reviewed `.assembly` + liftover AGP + contigs | `04_t2t_output/<sample>_T2T.FINAL.fa` |
| 05 | `05_run_gapcloser.sh` | FINAL.fa + ONT | LR_Gapcloser (3 iterations, keep last) |
| 06 | `06_gap_telo.sh` + `06_gap_telo_lib.py` | Gap-closed FA + ONT | Telomere extension + local gap fill |
| 07 | `07_polish_rename.sh` + `07_polish_rename_lib.py` | Filled FA + HiFi | Racon on fill flanks; chromosomes named vs Wm82/ZH13 → `07_polish_rename/<sample>/03_final/<sample>_Chr.fa` |
| 08 | `08_fa_qc.sh` + `08_fa_qc_lib.py` | Chr.fa + HiFi | Merqury QV, BUSCO, QUAST-LG, synteny dotplots |
| 09 | `09_run_hic_heatmap_cooler.sh` | Chr.fa + Hi-C | cooler matrix + publication heatmap |
| 10 | `10_run_centro.sh` | Chr.fa + TE GFF | quarTeT CentroMiner (run after annotation EDTA) |

`check_gap.sh` reports runs of N in a FASTA (1-based closed intervals).

## Design notes

- hifiasm is run in Hi-C + ultra-long mode (`--h1/--h2` and `--ul`) to phase and to traverse centromeric / repetitive knots.
- YaHS mapping uses `bwa mem -5SP` as recommended by the authors, then `samblaster` duplicate marking and name-sorted BAM.
- **Juicebox review is required.** The pipeline does not auto-correct chimeric joins.
- Gap closing uses ONT via LR_Gapcloser; remaining chromosome ends and small gaps are patched from ONT alignments in step 06 (quarTeT TeloExplorer + minimap2 baits).
- Chromosome IDs in the final FASTA follow soybean convention (Chr01–Chr20) by alignment to Wm82 / ZH13, not by contig order from hifiasm.
- QC is publication-oriented: k-mer QV (Merqury, k=21), BUSCO `embryophyta_odb10`, QUAST vs Wm82, and Hi-C contact maps at 500 kb (default).

## Work directory created by the scripts

All numeric prefixes (`01_assembly_out`, `02_yahs_results`, `03_JBAT_output`, … `10_run_centro`) are written under `ROOT`, not inside this git clone, unless you set `ROOT` to the clone (not recommended).

---

## 组装说明（中文）

按上表 01→10 顺序执行。03 之后必须在 Juicebox 里人工审查，再跑 04。07 的 `*_Chr.fa` 才是注释、QC、热图和着丝粒的正式基因组。10 依赖 EDTA/RepeatMasker 的 TE GFF，应在注释 01 完成后再跑。

hifiasm 使用 HiFi + 超长 ONT（`--ul`）+ Hi-C 分型；YaHS 比对必须带 `-5SP`。填缝用 ONT（LR_Gapcloser），端粒与残留 gap 在 06 用 quarTeT + minimap2 补。染色体名称按 Wm82/ZH13 定向，不是按 contig 编号。
