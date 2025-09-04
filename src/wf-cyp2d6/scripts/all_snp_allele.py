from pathlib import Path
import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")


# 样本变异表
vcf = snakemake.input[0]
name = Path(vcf).stem.split('.')[0]
# 跳过 VCF ## 的表头
lines = [line for line in open(vcf).readlines() if not line.startswith('##')]
matrix = [line.strip().split('\t') for line in lines]
# 读 '#CHROM', 'POS', 'ID', 'REF', 'ALT', 和突变信息列
vcf_df = pd.DataFrame(matrix[1:], columns=matrix[0])[['#CHROM', 'POS', 'ID', 'REF', 'ALT', name]]
# 样本名统一转成 Detail
vcf_df.rename(columns={name: 'DETAIL'}, inplace=True)
vcf_df['GENOTYPE'] = vcf_df['DETAIL'].str.split(':').str[0]
# pharmgkb 参考表
pharmgkb = pd.read_csv(snakemake.input[1], dtype=str)
# * 输出
# 基于样本变异合并
outdf = pd.merge(vcf_df, pharmgkb, on=['#CHROM', 'POS', 'REF', 'ALT'], how='left').fillna('-')
outdf.to_csv(snakemake.output[0], sep='\t', index=False)
