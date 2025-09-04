"""
根据 CNV 和 allele1/allele2 推断 CYP2D6  Diplotype, 并注释出 Phenotype
"""

import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")

# * 输入
# 分型信息汇总
df = pd.read_csv(snakemake.input[0], sep='\t')
# deplotype phenotype 对照字典
deplo_pheno_dict = pd.read_csv(
    snakemake.input[1],
    usecols=['CYP2D6 Diplotype', 'Coded Diplotype/Phenotype Summary']
).set_index('CYP2D6 Diplotype').to_dict()['Coded Diplotype/Phenotype Summary']

# 添加 deplotype 双倍型
for idx, row in df.iterrows():
    if (row.CNV_TYPE == 'DEL') and (row.CNV_RATIO < float(snakemake.params.ratio_cutoff_low)):
        sorted_alleles = [f'*{i}' for i in sorted([5, row.ALLELE1])]
    else:
        sorted_alleles = [f'*{i}' for i in sorted([row.ALLELE1, row.ALLELE2])]
    # df.at[idx, 'DEPLOTYPE'] = '/'.join(sorted_alleles)
    deplotype = '/'.join(sorted_alleles)
    # phenotype 表型
    if deplotype in deplo_pheno_dict:
        phenotype = deplo_pheno_dict[deplotype]
    else:
        phenotype = 'Unknown'
    df.at[idx, 'DEPLOTYPE'] = deplotype
    df.at[idx, 'PHENOTYPE'] = phenotype

# * 输出倍型和表型结果
outdf = df[['SAMPLE', 'DEPLOTYPE', 'PHENOTYPE',
            'CNV_TYPE', 'EXON_NUMBER', 'CNV_RATIO', 'ALLELE1', 'ALLELE2', 'TOP_UNIQUE']]
outdf.to_csv(snakemake.output[0], sep='\t', index=False)
outdf.to_excel(snakemake.output[1], index=False)
