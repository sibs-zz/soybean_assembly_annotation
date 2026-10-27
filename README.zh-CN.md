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
    N09["09 功能注释 + 改名"]
    N10["10 最终注释：修补、去冗余、反义筛选、短片段、BUSCO"]
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
  N02 --> N04 --> N05 --> N06 --> N07 --> N08 --> N09 --> N10
  N03 --> N07
```

**发表用基因集**是 10 步写成的 `13_final/`，不是 09 步的 high-confidence 子集。三类基因互不重叠，共同条件是 Chr01–Chr20、置信度 ≥ medium、TE 重叠 ≤ 0.5：

| 类别 | 附加条件 |
|------|----------|
| `high_confidence` | 相对 ZH13 v3 的相似度 ≥ 50%，不要求表达 |
| `expressed_zh13_below_50` | 有表达，相对 ZH13 v3 的相似度 < 50%；不用 ZH13 全长蛋白改写 |
| `homolog_unexpressed` | 没有表达，相对 ZH13 v3 的相似度 < 50%，并且至少一套大豆参考蛋白的最高 bitscore 命中达到相似度 ≥ 50%、双向覆盖 ≥ 50% |

10 步会修补 high-confidence 短模型，去掉同链重复和碎片，反义 CDS 重叠只保留两侧都有可靠同源的基因对，去掉已有完整拷贝的无表达短片段，最后做蛋白 BUSCO。位点编号不重排。`12_renamed/` 保留改名后的模型和删除清单。

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
├── annotation/                    # EDTA → BRAKER → EVM → polish → 最终注释
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

注释 09 的单样本写法是 `SAMPLE=N2L_primary THREADS=16 bash 09_run_polish.sh`。10 用品种前缀：`SAMPLE=KD35 bash 10_final_anno.sh`。不跑 BUSCO 时加 `RUN_BUSCO_FINAL=0`。

分步说明：[docs/assembly.md](docs/assembly.md)、[docs/annotation.md](docs/annotation.md)、[docs/outputs.md](docs/outputs.md)。

---

## 论文主注释用哪几个文件

09 步把基因改成 ZH13 风格 ID（如 `KD35_01G000100.m1`），结果在 `12_renamed/`。这一层保留改名后的全量注释和删除清单，用来追溯，不作为投稿基因集。

10 步（`10_final_anno.sh`，辅助脚本 `curate_final_models.py`、`resolve_final_redundancy.py`、`antisense_qc.py` 和 `short_fragment_filter.py`）读取 `12_renamed`，在 `13_final/` 写成一套注释：

1. 再查一遍 ORF。CDS 长度不能被 3 整除就删除整个基因；没有 CDS、没有 ATG、有内部终止或没有终止密码子就删除该转录本。09 步已经删除的基因不会回来。位点编号不重排。
2. 保留三类互不重叠的核基因（Chr01–Chr20，置信度至少 medium，转座元件重叠不超过 0.5）：
   - `high_confidence`：相对 ZH13 v3 的相似度 ≥ 50，不要求表达
   - `expressed_zh13_below_50`：有表达，相对 ZH13 v3 的相似度 < 50
   - `homolog_unexpressed`：没有表达，相对 ZH13 v3 的相似度 < 50，并且五套大豆参考蛋白里至少一套的最高 bitscore 命中同时满足相似度 ≥ 50、查询覆盖 ≥ 50%、参考覆盖 ≥ 50%、E-value ≤ 1e-5
3. 只修补 high-confidence 短模型。同位点 miniprot 能恢复至少 80% 的 ZH13 蛋白、没有移码和框内终止时，才提出新 CDS。新 CDS 必须有 ATG、有终止密码子、长度能被 3 整除、没有内部终止。任一条件不满足就保留原模型。同一次修订里两个新 CDS 在同链重叠达到 200 bp 时，这两个延长都不写入。`expressed_zh13_below_50` 不用 ZH13 全长蛋白改写。
4. 仍在这一套 `13_final` 里收尾：CDS 和蛋白完全相同的重复基因号合并，保留带 UTR 转录本的基因号。同链、同一 ZH13 基因，且一个基因的每条 CDS 坐标都已经是另一个基因的转录本时，删除多余基因号，不加成异构体。延长后才发生同链 CDS 重叠的模型恢复原 CDS。没有表达、短 CDS 和短蛋白都落在同一 ZH13 基因的长模型里的同链碎片从主注释去掉。
5. 反义 CDS 重叠只保留两侧都有可靠同源的基因对。ZH13 v3 上也有同样反义 CDS 重叠的记为较高；ZH13 没有同样重叠的记为中等至较高。缺少可靠同源的单外显子反义模型，以及其余未归入这两档的反义对，从这一套注释排除，基因号记在 `*.final.antisense_excluded.tsv`。位点编号不重排。
6. 每个基因取最终最长蛋白，比对 ZH13 v3 和 Wm82 a6。没有表达的 `high_confidence` 基因，如果对相似度 ≥ 50 的参考蛋白覆盖都不到一半，而同一个 ZH13 基因在本基因组已有覆盖 ≥ 80% 的拷贝，就从主注释去掉。它位于该 ZH13 基因所在染色体、且比所有完整拷贝更靠近 ZH13 位置时，视为同位点模型保留。全部候选及处理结果记在 `*.final.short_fragment_qc.tsv`。

**对外报基因数、投注释**用 `13_final/`：

- `*.final.gff3`
- `*.final.transcripts.fa`
- `*.final.cds.fa`
- `*.final.proteins.fa`
- `*.final.genes.list`
- `*.final.gene_class.tsv`

`*.final.corrections.tsv` 记录每个被审查的短模型；改过的行同时保存旧蛋白和新蛋白。`*.final.id_map.tsv` 只记录发生替换的基因号或转录本号。筛选规则在 `*.final.criteria.txt`。每个基因取最长蛋白做的 BUSCO 摘要在 `*.final.busco.short_summary.txt`。

重复序列仍在 `12_renamed/*.renamed.repeat.gff3.gz`。TE ID 不改成基因号。不要用 `12_renamed` 的全量 GFF 或 `11_publication/` 报基因数。

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
