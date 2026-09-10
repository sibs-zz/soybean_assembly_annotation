# 大豆 T2T 基因组组装与基因注释流程

**English version: [README.md](README.md).**

本仓库提供大豆（*Glycine max*）端粒到端粒（T2T）组装与证据驱动基因注释的可重复流程，在品种基因组（KD35、N2L 及后续材料）上用 PacBio HiFi、超长 ONT、Hi-C 以及短读 RNA-seq 跑通。

仓库只包含**脚本、辅助代码、环境说明和流程文档**，不含原始测序、组装序列或注释 GFF/FASTA。数据放在独立工作目录，用 `ROOT` 指向即可。

---

## 流程概览

```mermaid
flowchart TD
  subgraph reads["输入数据"]
    HiFi["PacBio HiFi"]
    ONT["超长 ONT"]
    HiC["Hi-C"]
    RNA["RNA-seq Atlas"]
  end

  subgraph asm["T2T 组装"]
    A01["01 hifiasm --ul + Hi-C"]
    A02["02 YaHS 挂载"]
    A03["03 Juicebox 人工纠错"]
    A04["04 juicer post"]
    A05["05 LR_Gapcloser"]
    A06["06 端粒 / gap 填补"]
    A07["07 Racon 抛光 + 染色体改名"]
    A08["08 Merqury / BUSCO / QUAST"]
    A09["09 Hi-C 热图"]
    A10["10 着丝粒 quarTeT"]
  end

  subgraph ann["基因注释"]
    N01["01 EDTA + RepeatMasker"]
    N02["02 HISAT2 + LAI"]
    N03["03 BRAKER3"]
    N04["04 StringTie + Trinity-GG"]
    N05["05 Mikado"]
    N06["06 PASA 比对"]
    N07["07 EVM"]
    N08["08 PASA UTR 更新"]
    N09["09 功能注释 + HC + 改名"]
  end

  HiFi --> A01
  ONT --> A01
  HiC --> A01
  A01 --> A02 --> A03 --> A04 --> A05 --> A06 --> A07
  A07 --> A08
  A07 --> A09
  N01 --> A10
  A07 --> N01 --> N02 --> N03
  RNA --> N02
  N02 --> N04 --> N05 --> N06 --> N07 --> N08 --> N09
  N03 --> N07
```

**发表用高置信（HC）基因集**（09 步）：

| 过滤 | 阈值 |
|------|------|
| 染色体 | 仅 Chr01–Chr20 |
| 置信度 | ≥ medium |
| TE 重叠比例 | ≤ 0.5 |
| 相对 ZH13 v3 蛋白 identity | ≥ 50% |

另输出互斥的**品种特异集**（有表达、TE 干净、ZH13 identity < 50%），供补充材料和差异分析。全量核染色体 GFF 仍保留这些位点，只是不计入 HC 基因数。

---

## 目录结构

```
soybean_assembly_annotation/
├── README.md / README.zh-CN.md
├── LICENSE / CITATION.cff / .gitignore
├── config/config.env.example      # 路径与 HC 过滤参数
├── env/                           # conda 环境
├── docs/                          # 方法说明
├── assembly/                      # 组装 01–10
├── annotation/                    # EDTA → BRAKER → EVM → polish
└── preprocess/                    # 可选：HiFi / Hi-C / RNA 清洗
```

每个脚本顶部有**配置区**。其中的 `/tiandata2/`、`/home/zzh/` **只是示例路径，必须改成你自己的目录**，或先 `export ROOT=...`、`export ANACONDA=...`（其余变量见 `config/config.env.example`）。基因组、reads、注释 GFF/FASTA 不要提交进本仓库。

---

## 依赖

- Linux 服务器（各步常用 16–64 线程）
- Conda / Mamba（见 [env/README.md](env/README.md)）
- 部分工具不在 conda 里：YaHS、juicer_tools、LR_Gapcloser、quarTeT、GeneMark-ETP、ProtHint、InterProScan
- 参考：Wm82 / ZH13 基因组、Viridiplantae 蛋白库、SwissProt Diamond 库、BUSCO `embryophyta_odb10`、五套大豆参考蛋白（ZH13 v2/v3，Wm82 a2/gnm4/gnm6）

GeneMark-ETP 需学术许可，运行 BRAKER3 前请自行申请。

---

## 快速开始

```bash
git clone https://github.com/sibs-zz/soybean_assembly_annotation.git
cd soybean_assembly_annotation
cp config/config.env.example config/config.env
# 必须自己改：ROOT、conda、YaHS/juicer/EDTA、参考基因组等
# 脚本里仍有 /tiandata2/、/home/zzh/ 示例路径，请替换，或 export 下列变量覆盖 ${VAR:-default}

export ROOT=/path/to/workdir
export ANACONDA=/path/to/anaconda3/envs
export THREADS=32
```

工作目录最低需要：

```
$ROOT/hifi/          HiFi FASTQ/FASTA
$ROOT/ont/           清洗后的超长 ONT
$ROOT/hic/           清洗后的 Hi-C 双端
$ROOT/rna_clean/     Atlas RNA-seq 双端
```

按编号顺序跑：先组装 01–07 得到 `*_Chr.fa`，再注释。组装 10 着丝粒需要 EDTA 的 TE GFF，应在注释 01 之后。

单样本：`SAMPLE=N2L_primary THREADS=16 bash 09_run_polish.sh`。

分步说明：[docs/assembly.md](docs/assembly.md)、[docs/annotation.md](docs/annotation.md)、[docs/outputs.md](docs/outputs.md)。

---

## 论文主注释用哪几个文件

09 步会写成 ZH13 风格 ID（如 `KD35_01G000100.m1`）。**对外报基因数、投注释**用 `12_renamed/` 下：

- `*.renamed.high_confidence.annotated.gff3`
- `*.renamed.high_confidence.proteins.fa`
- `*.renamed.high_confidence.genes.list`

不要用全量 nuclear GFF（约 8 万位点）报基因数。`11_publication/` 若仍是旧过滤，也不要当最终 HC。

特异集文件名为 `*.renamed.lineage_specific.*`（list / GFF / mRNA / CDS / 蛋白）。

---

## 引用

请引用后续发表的大豆基因组论文（如有）以及本仓库。Methods 中请同时引用 hifiasm、YaHS、BRAKER3、PASA、EVM、EDTA、eggNOG、InterProScan、BUSCO 等工具原文。

---

## 许可

本仓库代码为 [MIT License](LICENSE)。第三方软件仍按其原许可使用（尤其是 GeneMark-ETP）。

---

## 复现时注意

1. Juicebox **人工纠错**是流程的一部分：03 生成 JBAT，04 读取 review 后的 `.assembly`。
2. 多个 conda 环境是有意拆开的（Trinity 需要 salmon 1.x；eggNOG 需要 Diamond ≥2.2 的 `--iterate`）。见 [env/README.md](env/README.md)。
3. `09_run_polish.sh` 默认 64 线程；注释 04–08 建议 `THREADS=16`。
4. **跑之前请自己配路径**：在脚本里搜索 `/tiandata2/`、`/home/zzh/` 并替换，或按 `config/config.env.example` `export ROOT` / `ANACONDA`。基因组和注释结果放在工作目录，不要放进 git。
