# KML-CYP2D6 分析流程

## 项目概述

基于靶向 NGS 数据的 CYP2D6 / CYP2C19 / CYP2C9 药物基因组学 (PGx) 基因分型流程。使用 PharmGKB 等位基因定义表构建参考库，不使用 PharmVAR/PharmCAT。

## 流程调用方式

```bash
python -m src.kml_cyp2d6 \
  --input-tab template/input.tsv \
  --output-dir /path/to/results \
  --threads 32
```

输入格式 (TSV/XLSX，三列，无表头):

| sample_name | R1_fastq_path | R2_fastq_path |
|-------------|---------------|---------------|

CLI 入口: `src/kml_cyp2d6/cli.py:19` -> `prepare_fastq_by_samptab()` + `run_snakemake()`

流程包含两步:

1. FASTQ 准备 — 将输入 FASTQ 拷贝/压缩到 `.rawdata/{sample}_{1,2}.fastq.gz`
2. Snakemake 工作流 — 生成配置文件并启动 `src/wf-cyp2d6/Snakefile`

---

## 完整分析流程 (Snakemake DAG)

### Phase 0: 配置与准备

生成 `snakemake.yaml` 配置文件，聚合以下模块:

- `src/config/database.py` — 参考基因组、GATK known sites、BED、PharmGKB 表、CoNVaDING 控件路径
- `src/config/env.py` — Conda 环境名 (basic/basic2/python3.12/convading)
- `src/config/software.py` — Activate 脚本路径、CoNVaDING Perl 路径

配置文件生成: `src/kml_cyp2d6/snakemake.py:14`

### Phase 1: 预处理

| 步骤 | 规则文件 | 规则 | 工具/脚本 | 说明 |
|------|---------|------|-----------|------|
| BED 排序 | `bed.smk` | `bedtools_sort` | bedtools sort | 对探针 BED 排序 |
| 原始 QC | `fastqc.smk` | `fastqc` | FastQC | 每个样本原始 FASTQ QC |
| 聚合 QC | `fastqc.smk` | `multiqc` | MultiQC | 合并 FastQC 报告 |
| 质控过滤 | `fastp.smk` | `fastp` | fastp | `-q 15 -u 40 -l 15 --correction`，接头修剪+质量过滤 |
| 汇总 | `fastp.smk` | `fq_stats_summary` | `fq_all_samples_qc.py` | 聚合 fastp JSON 为汇总 TSV |

### Phase 2: 比对

| 步骤 | 规则 | 工具/脚本 | 说明 |
|------|------|-----------|------|
| `bwa_mem` | bwa mem + samtools sort | 带 RG 标签比对到 hg38 |
| `samtools_index` | samtools index | BAM 建索引 |
| `samtools_stats` | samtools stats | 靶向区域统计 |
| `samtools_stats_all` | samtools stats | 全基因组统计 |
| `samtools_depth` | samtools depth | 靶向区域深度 |
| `bam_stats` | `bam_stats.py` | 合并 stats+depth 为 per-sample CSV |
| `bam_stats_summary` | `bam_stats_summary.py` | 聚合所有样本 BAM 统计 |
| `samtools_bedcov` | samtools bedcov | 靶向区域覆盖度 |

### Phase 3: BQSR (碱基质量重校准)

| 步骤 | 规则 | 工具/脚本 | 说明 |
|------|------|-----------|------|
| `mark_duplicates` | Picard MarkDuplicates | `REMOVE_DUPLICATES=true` |
| `recalibrate_base_qualities` | GATK BaseRecalibrator | 已知 SNP/INDEL 位点重校准 |
| `apply_base_quality_recalibration` | GATK ApplyBQSR | BQSR 应用到 BAM |

### Phase 4: 变异检测

| 步骤 | 规则 | 工具/脚本 | 说明 |
|------|------|-----------|------|
| `haplotype_caller` | GATK HaplotypeCaller | 靶向 BED 区域检测变异 |

### Phase 5: 变异过滤

| 步骤 | 规则 | 工具/脚本 | 说明 |
|------|------|-----------|------|
| `gatk_select_snps` | GATK SelectVariants | 提取 SNPs |
| `gatk_select_indels` | GATK SelectVariants | 提取 INDELs |
| `gatk_filter_snps` | GATK VariantFiltration | `QD < 2.0 \|\| FS > 60.0 \|\| MQ < 40.0` |
| `gatk_filter_indels` | GATK VariantFiltration | `QD < 2.0 \|\| FS > 200.0` |
| `merge_vcfs` | Picard MergeVcfs | 合并过滤后的 SNP+INDEL |
| `bcftools_view` | bcftools view | 仅保留 PASS 变异 |
| `bcftools_filter_het` | bcftools filter | HET: VAF 0.25~0.75 |
| `bcftools_filter_hom` | bcftools filter | HOM: VAF > 0.9 |
| `merge_het_hom_filter_vcfs` | Picard MergeVcfs | 合并 HET+HOM VCF |

### Phase 6: CNV 检测

使用 CoNVaDING 进行拷贝数变异检测:

| 步骤 | 规则 | 说明 |
|------|------|------|
| `start_with_bam` | 归一化，使用预构建 control 数据集 |
| `start_with_match_score` | 选择最佳 control 样本 (controlSamples=30) |
| `start_with_best_score` | CNV 检测 (ratioCutOffLow=0.65, ratioCutOffHigh=1.4) |
| `create_final_list` | 生成最终 CNV 列表 |
| `parse_cnv` | `parse_cnv.py` | 解析 CNV 结果为 `cnv.txt` |

### Phase 7: 等位基因鉴定 (CYP2D6)

| 步骤 | 规则 | 脚本 | 说明 |
|------|------|------|------|
| `analysis_allele_snps` | `analysis_allele_snps.py` | 与 PharmGKB 定义比较，计算 PRESENT RATIO |
| `all_snp_allele` | `all_snp_allele.py` | 列出所有检测到的 SNP 等位基因 |
| `call_raw_allele` | `call_raw_allele.py` | 结合 SNP stats + 突变注释，推断原始 star-allele |
| `paste_allele_cnv` | paste 命令 | 合并等位基因 + CNV 结果 |

### Phase 8: 报告生成 (CYP2D6)

| 步骤 | 规则 | 脚本/工具 | 说明 |
|------|------|-----------|------|
| `report_loci_info` | `report_loci_info.py` | 位点级报告 (含 CNV ratio) |
| `csv2xlsx_allele_snp` | csvtk csv2xlsx | 转为 XLSX |
| `summary_allele` | `summary_allele.py` | 等位基因汇总 |
| `summary_diplotype_phenotype` | `summary_diplo_pheno.py` | 推断双倍型+表型 |
| `all_summary` | `all_summary.py` | FASTQ + BAM QC 汇总 -> `panel-qc-summary.xlsx` |

### Phase 9: CYP2C19 分型

规则文件: `cyp2c19.smk`

- 复用 CYP2D6 的 allele calling 规则，使用 CYP2C19 参考表
- 默认等位基因: `*38` (而非 `*1`)
- `*36` = 全基因缺失, `*37` = 部分外显子缺失
- CNV 解析: `parse_cnv_c19.py` (ratioCutOffLow=0.65)
- 双倍型/表型: `summary_diplo_pheno_c19.py`

### Phase 10: CYP2C9 分型

规则文件: `cyp2c9.smk`

- 复用 CYP2C19 的 call_raw_allele 规则，使用 CYP2C9 参考表
- 默认等位基因: `*1`
- 无双倍型相关的 DEL/CNV

---

## 基因分型核心逻辑

### PRESENT RATIO 计算 (`analysis_allele_snps.py`)

```
PRESENT RATIO = 分型测到的位点数 / 分型中理论的位点数
- 理论位点数量固定(不含简并碱基位点)
- 简并碱基位点替换成变异碱基参与计算
- 数值超过 1 则取 1
```

### VAF 过滤阈值 (`filter.smk`)

| 基因型 | VAF 条件 |
|--------|----------|
| HOM | VAF > 0.9 |
| HET | 0.25 < VAF < 0.75 |

### GATK 硬过滤阈值

| 变异类型 | 过滤条件 |
|----------|----------|
| SNP | `QD < 2.0 \|\| FS > 60.0 \|\| MQ < 40.0` |
| INDEL | `QD < 2.0 \|\| FS > 200.0` |

### CoNVaDING ratio cutoffs

| 参数 | CYP2D6 | CYP2C19 |
|------|--------|---------|
| ratioCutOffLow | 0.65 | 0.65 |
| ratioCutOffHigh | 1.4 | — |

---

## 输出目录结构

```
output-dir/
├── .rawdata/                    # 准备的 FASTQ
├── .temp/                       # 临时文件 + snakemake.yaml + snakemake.log
├── multiqc/                     # MultiQC 报告
├── fastp/                       # 质控后 FASTQ + JSON
├── fastqc/                      # FastQC 报告
├── align/                       # BAM + stats + depth + bedcov
├── bqsr/                        # 去重 BAM + BQSR 重校准
├── calls/                       # GATK HaplotypeCaller VCF
├── filter/                      # 过滤后 VCF
├── convading/{sample}/          # CoNVaDING CNV 结果
│   ├── StartWithBam/
│   ├── StartWithMatchScore/
│   ├── StartWithBestScore/
│   ├── CreateFinalList/
│   └── cnv.txt
├── allele/                      # 等位基因明细
├── upload/                      # 最终报告
│   ├── panel-qc-summary.xlsx
│   ├── all.allele.summary.tsv
│   ├── all.diplotype_phenotype.tsv
│   ├── all.diplotype_phenotype.xlsx
│   └── {sample}.summary.xlsx
├── cyp2c19/                     # CYP2C19 结果 (结构与 upload/ 类似)
│   └── all.diplotype_phenotype.tsv
├── cyp2c9/                      # CYP2C9 结果 (结构与 upload/ 类似)
│   └── all.diplotype_phenotype.tsv
└── .log/                        # 日志 + benchmark
```

## 数据库配置

硬编码路径位于 `src/config/database.py`，指向服务器路径:

- 参考基因组: hg38
- Known sites: 1000G / dbSNP / Mills (GATK 资源)
- PharmGKB 定义表: `assets/CYP2D6/allele_definition_table_reportable.csv`

## Conda 环境

| 名称 | 用途 |
|------|------|
| basic | 基础生物信息工具 (bwa, samtools, bcftools, GATK 等) |
| basic2 | 额外工具 (fastp, Picard, bedtools, csvtk 等) |
| python3.12 | Python 分析脚本 (pandas, numpy, openpyxl, biopython 等) |
| convading | CoNVaDING Perl 环境 |

---

## 开发注意

1. 不使用 Poetry — 在 snakemake 中环境冲突 (numpy/pandas 版本冲突)
2. 不使用 PharmVAR — 使用 PharmGKB 对照表
3. 不使用 Docker/Singularity — 依赖 Conda 环境
4. 脚本格式化: `snakefmt src/wf-cyp2d6`
5. 测试: `python -m pytest tests/` (生成 snakemake config 验证)
