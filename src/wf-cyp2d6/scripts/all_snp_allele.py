from pathlib import Path
import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")


# * 样本变异表
vcf = snakemake.input[0]
name = Path(vcf).stem.split('.')[0]
# 跳过 VCF ## 的表头
lines = [line for line in open(vcf).readlines() if not line.startswith('##')]
matrix = [line.strip().split('\t') for line in lines]
# 读 '#CHROM', 'POS', 'ID', 'REF', 'ALT', 和突变信息列
curdf = pd.DataFrame(matrix[1:], columns=matrix[0])[['#CHROM', 'POS', 'ID', 'REF', 'ALT', name]]
# 样本名统一转成 Detail
curdf.rename(columns={name: 'DETAIL'}, inplace=True)

# * pharmvar 参考表
pharmvar = pd.read_csv(snakemake.input[1], dtype=str)

# * 基于样本变异合并
outdf = pd.merge(curdf, pharmvar, on=['#CHROM', 'POS', 'REF', 'ALT'], how='left').fillna('-')
# outdf['CLEAN_ALLELE'] = outdf['ALLELE'].apply(lambda x: x.split('.')[0])
# outdf = outdf.groupby(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'DETAIL'])['CLEAN_ALLELE'].apply(lambda x: ','.join(set(x))).reset_index()
# outdf.rename(columns={'CLEAN_ALLELE': 'ALLELE'}, inplace=True)

# * 输出
outdf.to_csv(snakemake.output[0], sep='\t', index=False)
