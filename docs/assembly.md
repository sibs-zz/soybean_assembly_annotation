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
| 06 (extra) | `06_gap_telo_final_fix.sh` + `06_gap_telo_final_fix_lib.py` | 06 output + ONT | Buried-telomere trim / tip extension; replaces `06_gap_telo/<sample>/04_final` only when quarTeT finds all 40 ends and no gap (old copy kept as `04_final.bak_*`). Run for DY11 and TH2 |
| 07 | `07_polish_rename.sh` + `07_polish_rename_lib.py` | Filled FA + HiFi | Racon on fill flanks; chromosomes named vs Wm82/ZH13 → `07_polish_rename/<sample>/03_final/<sample>_Chr.fa` |
| 08 | `08_fa_qc.sh` + `08_fa_qc_lib.py` | Chr.fa + HiFi | Merqury QV, BUSCO, QUAST-LG, synteny dotplots |
| 09 | `09_run_hic_heatmap_cooler.sh` | Chr.fa + Hi-C | cooler matrix + publication heatmap |
| 10 | `10_run_centro.sh` | Chr.fa + TE GFF | quarTeT CentroMiner (run after annotation EDTA) |

`check_gap.sh` reports runs of N in a FASTA (1-based closed intervals).

`assembly/` holds the scripts run for KD35 and N2L. `assembly/batch2/` holds the 01–10 drivers run for BZ, DY11, TH2, TZX1998 and TZX253. They differ in sample lists and input names: `01_prep_hifi_hic.sh` converts HiFi BAM to `hifi/<cultivar>_hifi.{fq.gz,fasta}` (the two TH2 HiFi runs are merged) and trims Hi-C to `hic/<cultivar>_hic_{1,2}.fq.gz`; 01 and 02 loop over the samples; 07 and 08 take either HiFi format; 09 reads the `_hic_` names. They call the Python helpers and `check_gap.sh` in `assembly/`.

## Design notes

- hifiasm is run in Hi-C + ultra-long mode (`--h1/--h2` and `--ul`) to phase and to traverse centromeric / repetitive knots.
- YaHS mapping uses `bwa mem -5SP` as recommended by the authors, then `samblaster` duplicate marking and name-sorted BAM.
- **Juicebox review is required.** The pipeline does not auto-correct chimeric joins.
- Gap closing uses ONT via LR_Gapcloser; remaining chromosome ends and small gaps are patched from ONT alignments in step 06 (quarTeT TeloExplorer + minimap2 baits).
- Chromosome IDs in the final FASTA follow soybean convention (Chr01–Chr20) by alignment to Wm82 / ZH13, not by contig order from hifiasm.
- QC is publication-oriented: k-mer QV (Merqury, k=21), BUSCO `embryophyta_odb10`, QUAST vs Wm82, and Hi-C contact maps at 500 kb (default).

## Work directory created by the scripts

All numeric prefixes (`01_assembly_out`, `02_yahs_results`, `03_JBAT_output`, … `10_run_centro`) are written under `ROOT`, not inside this git clone, unless you set `ROOT` to the clone (not recommended).

## Results of this run

From `08_fa_qc/<sample>/05_summary/qc_metrics.tsv` and `annot/<cultivar>_TE/*.LAI`. Every assembly has 20 chromosomes with no N; the rest are unanchored contigs.

| Cultivar | Chr Mb | Unanchored | QV | k-mer completeness | BUSCO C | LAI | Telomeres (strict) |
|---|---:|---:|---:|---:|---:|---:|---:|
| KD35 | 1012.3 | 111 (13.0 Mb) | 53.1 | 99.51% | 99.7% | 9.80 | 40/40 |
| N2L | 1006.5 | 137 (11.3 Mb) | 54.8 | 99.62% | 99.6% | 11.30 | 37/40 |
| BZ | 1006.6 | 99 (8.6 Mb) | 54.1 | 99.64% | 99.6% | 10.57 | 39/40 |
| DY11 | 1008.5 | 91 (8.9 Mb) | 56.3 | 99.54% | 99.6% | 10.56 | 37/40 |
| TH2 | 1006.2 | 64 (6.0 Mb) | 53.9 | 99.58% | 99.7% | 11.99 | 39/40 |
| TZX1998 | 1013.8 | 92 (8.8 Mb) | 54.8 | 99.55% | 99.6% | 11.46 | 38/40 |
| TZX253 | 1009.6 | 114 (9.8 Mb) | 54.7 | 99.54% | 99.7% | 12.50 | 40/40 |

Telomere counts. quarTeT TeloExplorer at `-m 50`, used inside step 06, reports 40/40 for every assembly. The strict count above needs at least 50 TTTAGGG/CCCTAAA copies in the outer 20 kb of the chromosome end; it gives 270 of 280 ends. Of the other ten, the right end of Chr13 in N2L and DY11 has no telomere; the only nearby array is about 40 CCCTAAA copies 20–30 kb inside, in the orientation of a left end. Six ends have the array 19–170 kb inside the tip (BZ Chr13 right, TH2 Chr12 left, TZX1998 Chr09 right and Chr15 left, N2L Chr17 left and Chr19 right). DY11 Chr04 right and Chr18 right have 29–32 copies at the tip. Report the strict count, not the quarTeT count.

45S rDNA. In every assembly the hifiasm contigs of Chr13 break at about 3.2 Mb, where 13–43 copies of the 18S 5′ end (`TACCTGGTTGATCCTGCCAG`) sit; another 446–817 copies are on 41–82 unanchored contigs. The nucleolus organiser region on Chr13 is therefore collapsed, even though the chromosome has no N.

Chr11–Chr13 translocation. N2L, BZ, DY11, TH2 and TZX1998 carry a reciprocal translocation relative to ZH13. Their Chr11 is about 56 Mb and Chr13 about 32 Mb, against 41–42 Mb and 46 Mb in KD35, TZX253 and ZH13. Against ZH13, about 19 Mb of ZH13 Chr13 sits on their Chr11, and 4–5 Mb of ZH13 Chr11 sits on their Chr13. Every junction lies inside one hifiasm contig (for example, TH2 Chr11 is the single 56.24 Mb contig `ptg000018l`), so it is not a Hi-C scaffolding join. Chromosome names follow the larger share of reference support. Gene IDs follow the chromosome a gene sits on, so a gene in a moved segment can be `11G` while its ZH13 orthologue is on Chr13, or the reverse.

---

## 组装说明（中文）

按上表 01→10 顺序执行。03 之后必须在 Juicebox 里人工审查，再跑 04。07 的 `*_Chr.fa` 才是注释、QC、热图和着丝粒的正式基因组。10 依赖 EDTA/RepeatMasker 的 TE GFF，应在注释 01 完成后再跑。

hifiasm 使用 HiFi + 超长 ONT（`--ul`）+ Hi-C 分型；YaHS 比对必须带 `-5SP`。填缝用 ONT（LR_Gapcloser），端粒与残留 gap 在 06 用 quarTeT + minimap2 补。染色体名称按 Wm82/ZH13 定向，不是按 contig 编号。

`assembly/` 是 KD35、N2L 实际运行的脚本；`assembly/batch2/` 是 BZ、DY11、TH2、TZX1998、TZX253 实际运行的 01–10。两者的区别在样品列表和输入命名：`01_prep_hifi_hic.sh` 把 HiFi BAM 转成 `hifi/<品种>_hifi.{fq.gz,fasta}`（TH2 的两次 HiFi 合并），Hi-C 用 fastp 处理成 `hic/<品种>_hic_{1,2}.fq.gz`；01、02 按样品循环；07、08 两种 HiFi 格式都能读；09 读 `_hic_` 命名。Python 辅助脚本和 `check_gap.sh` 共用 `assembly/` 里的。DY11、TH2 在 06 之后又跑了 `06_gap_telo_final_fix.sh`：埋在内侧的端粒按 ONT 证据裁剪或外延，只有 quarTeT 判定 40 端齐全且无 gap 时才替换 `06_gap_telo/<样品>/04_final`，旧结果保留为 `04_final.bak_*`。

### 本次组装结果

七套都是 20 条染色体、染色体上没有 N，其余为未挂载 contig。QV、k-mer 完整度、BUSCO、LAI 和端粒数见上方英文表。

端粒。06 里用的 quarTeT TeloExplorer（`-m 50`）对七套都报 40/40。表里的严格口径要求染色体末端 20 kb 内至少 50 个 TTTAGGG/CCCTAAA，七套合计 270/280。其余 10 个末端中：N2L 和 DY11 的 Chr13 右端没有端粒，附近只有内侧 20–30 kb 处约 40 个 CCCTAAA，方向与左端相同；BZ Chr13 右端、TH2 Chr12 左端、TZX1998 Chr09 右端和 Chr15 左端、N2L Chr17 左端和 Chr19 右端的端粒在末端内侧 19–170 kb；DY11 Chr04 右端和 Chr18 右端末端只有 29–32 个拷贝。论文请用严格口径，不用 quarTeT 的 40/40。

45S rDNA。七套 Chr13 的 hifiasm contig 都在约 3.2 Mb 处断开，那里有 13–43 个 18S 5′ 端序列（`TACCTGGTTGATCCTGCCAG`），另有 446–817 个分布在 41–82 条未挂载 contig 上。Chr13 的核仁组织区是塌缩的，虽然染色体上没有 N。

Chr11–Chr13 易位。N2L、BZ、DY11、TH2、TZX1998 相对 ZH13 带一个相互易位：Chr11 约 56 Mb、Chr13 约 32 Mb，KD35、TZX253 和 ZH13 分别为 41–42 Mb 和 46 Mb。与 ZH13 比对，约 19 Mb 的 ZH13 Chr13 位于这五套的 Chr11 上，4–5 Mb 的 ZH13 Chr11 位于其 Chr13 上。所有接合点都在单条 hifiasm contig 内部（例如 TH2 的 Chr11 就是一条 56.24 Mb 的 contig `ptg000018l`），不是 Hi-C 挂载造成的。染色体按参考支持较多的一侧命名，基因号跟随所在染色体，所以易位片段上的基因可能是 `11G` 而 ZH13 直系同源在 Chr13，反之亦然。
