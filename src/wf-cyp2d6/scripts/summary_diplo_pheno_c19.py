import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")

# * 输入
# 分型信息汇总
df = pd.read_csv(snakemake.input[0], sep='\t')
# deplotype phenotype 对照字典
deplo_pheno_dict = pd.read_csv(
    snakemake.input[1],
    usecols=['Diplotype', 'Coded Diplotype/Phenotype Summary']
).set_index('Diplotype').to_dict()['Coded Diplotype/Phenotype Summary']

# 添加 deplotype 双倍型
for idx, row in df.iterrows():
    # CYP2C19 为 9 个外显子, 但靶区域为 12 个区段, 只能按照区段去看是部分外显子缺失还是全部缺失
    # CYP2C19*36 全基因缺失, *37 部分外显子缺失
    if row['DELETE_EXON_COUNT'] == 12:
        sorted_alleles = [f'*{i}' for i in sorted([row.ALLELE1, 36])]
    elif 0 < row['DELETE_EXON_COUNT'] < 12:
        sorted_alleles = [f'*{i}' for i in sorted([row.ALLELE1, 37])]
    else:
        sorted_alleles = [f'*{i}' for i in sorted([row.ALLELE1, row.ALLELE2])]
    diplotype = '/'.join(sorted_alleles)
    # phenotype 表型
    if diplotype in deplo_pheno_dict:
        phenotype = deplo_pheno_dict[diplotype]
    else:
        phenotype = 'Unknown'
    df.at[idx, 'DIPLOTYPE'] = diplotype
    df.at[idx, 'PHENOTYPE'] = phenotype

# * 输出双倍型和表型结果
outdf = df[['SAMPLE', 'DIPLOTYPE', 'PHENOTYPE',
            'DELETE_EXON_COUNT', 'ALLELE1', 'ALLELE2', 'TOP_UNIQUE']]
outdf.to_csv(snakemake.output[0], sep='\t', index=False)
outdf.to_excel(snakemake.output[1], index=False)
