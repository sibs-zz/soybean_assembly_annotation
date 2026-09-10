# Software environments

Do not merge these conda envs. Version conflicts are real (salmon 1 vs 2; Diamond 0.9 vs 2.2).

## Assembly / QC

| Env name | Typical use |
|----------|-------------|
| `qc_merqury` | Merqury / meryl k-mer QV |
| `qc_busco` | BUSCO 5.7.x, `embryophyta_odb10` |
| `qc_quast` | QUAST-LG |
| `qc_synteny` | minimap2, racon (polish) |
| `hic_tools` | bwa, pairtools, cooler, heatmap |
| `quarTeTdependencies` | quarTeT TeloExplorer / CentroMiner |
| `py37` | seqkit, matplotlib helpers |
| `nanopore` | seqkit (shared with annotation) |

System or local builds (set paths in `config/config.env`):

- [hifiasm](https://github.com/chhylp123/hifiasm)
- [YaHS](https://github.com/c-zhou/yahs) + `juicer` helper
- [juicer_tools](https://github.com/aidenlab/juicer) 1.9.9 jar
- [LR_Gapcloser](https://github.com/CAFS-bioinformatics/LR_Gapcloser)
- [quarTeT](https://github.com/zy-zhang/quarTeT)
- bwa, samtools, samblaster, minimap2

## Annotation

Set `ANACONDA` to the directory that contains these env folders (`annot_common.sh`).

| Env name | Typical use |
|----------|-------------|
| `EDTA` | RepeatMasker, LAI, AGAT, cd-hit-est |
| `aligner` | HISAT2, StringTie, miniprot, samtools |
| `braker3` | BRAKER3, `gtf2gff.pl` |
| `annot_gmap` | GMAP (Mikado / Trinity transcripts) |
| `annot_trinity` | Trinity 2.15.2 + **salmon 1.x** wrapper |
| `annot_transdecoder` | TransDecoder |
| `annot_pasa` | PASA 2.5.3 (`PASAHOME`) |
| `annot_mikado` | Mikado (pandas 2.x tested), Portcullis |
| `annot_evm` | EVM 2.1.0 (`EVM_HOME`) |
| `annp` | gffread, **Diamond 2.2.x**, bedtools, salmon (09 polish) |
| `eggnog213` | eggNOG-mapper 2.1.3 (Python); put `annp` Diamond on PATH |
| `qc_busco` | protein BUSCO in 09 |

Also required on disk:

- GeneMark-ETP and ProtHint (BRAKER3)
- InterProScan
- EDTA apptainer/singularity image (step 01)
- Diamond DB: SwissProt; Viridiplantae; five soybean proteomes

## Creating envs

Export from a working machine (example):

```bash
conda env export -n annot_pasa --from-history > env/annot_pasa.from-history.yml
```

`--from-history` is preferred for GitHub (no local build strings). Recreate:

```bash
mamba env create -f env/annot_pasa.from-history.yml
```

If you cannot export, install the named tools into identically named envs and keep `ANACONDA=/path/to/anaconda3/envs`.

## Diamond / salmon pitfalls

- System Diamond 0.9.x **cannot** run eggNOG `--iterate`. 09 polish must see Diamond 2.2+ (`annp`).
- Trinity 2.15.2 still passes `--minAssignedFrags` to salmon; salmon 2 removed that flag. Keep salmon 1.x on PATH for step 04.

## 中文摘要

组装 QC、Hi-C、注释不要塞进同一个 conda。Trinity 用 salmon 1.x，eggNOG 用 Diamond 2.2+。`ANACONDA` 指 env 的父目录。GeneMark-ETP、InterProScan、EDTA 容器、YaHS 等按 `config/config.env.example` 单独指定路径。
