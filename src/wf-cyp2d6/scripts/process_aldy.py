import pandas as pd
import json
import sys

sys.stderr = open(snakemake.log[0], "w")

# * IO
aldy_res = snakemake.input[0]
out_snp_tsv = snakemake.output[0]
out_genotype = snakemake.output[1]

# * genotype-rsid 参考
with open(snakemake.params[0]) as f:
    genotype_rsid_dict = json.load(f)

# * MAIN
# 只读第一个 solution
records = []
with open(aldy_res) as f:
    for line in f:
        if line.startswith('#Solution 1'):
            continue
        elif line.startswith('#Solution 2'):
            break
        else:
            records.append(line.split('\t')[:13])
df = pd.DataFrame(records[1:], columns=records[0])[
    ['Copy', 'Allele', 'Location', 'Type', 'Coverage', 'dbSNP']]
# coverage 有空值的过滤掉
df = df[df['Coverage'] != '']
df['Coverage'] = df['Coverage'].astype(int)
# ! Coverage 过滤阈值 10X
df = df[df['Coverage'] >= 10]
df['NewAllele'] = '*' + df['Allele'].str.split('.').str[0]
# * 过滤掉不在参考中的分型
df = df[df['NewAllele'].isin(genotype_rsid_dict)]
df['Allele'] = 'CYP2D6' + df['NewAllele']
# 输出表格
detected_alleles = df['NewAllele'].unique().tolist()
detected_alleles_snps = [
    rsid for alleles in detected_alleles for rsid in genotype_rsid_dict[alleles]]
# 去重
outdf = df[df['dbSNP'].isin(detected_alleles_snps)][['Allele', 'dbSNP',
                                                     'Location', 'Type', 'Coverage']].drop_duplicates()
outdf.to_csv(out_snp_tsv, index=False, sep='\t')
# 输出分型
genotype = 'CYP2D6 (*1/*1)'
if len(outdf) == 0:
    # *1 野生型
    genotype = 'CYP2D6 (*1/*1)'
else:
    unique_copies = df['Copy'].unique()
    unique_alleles = df['NewAllele'].unique()
    n_copies = len(unique_copies)
    n_alleles = len(unique_alleles)
    # *5 只有一个拷贝
    if n_copies == 1:
        genotype = f'CYP2D6 ({unique_alleles[0]}/*5)'
    # 两个拷贝
    elif n_copies == 2:
        if n_alleles == 1:
            # 纯合
            genotype = f'CYP2D6 ({unique_alleles[0]}/{unique_alleles[0]})'
        elif n_alleles == 2:
            # 杂合
            genotype = f'CYP2D6 ({unique_alleles[0]}/{unique_alleles[1]})'
        else:
            raise Exception('分型错误，查看 aldy 文件')
    # 多拷贝
    elif n_copies > 2:
        val_count = df[['Copy', 'NewAllele']].drop_duplicates()['NewAllele'].value_counts()
        alleles = [
            f'{allele}x{count}' if count > 1 else str(allele)
            for allele, count in val_count.items()
        ]
        genotype = 'CYP2D6 (' + '/'.join(alleles) + ')'

# 分型和 SNP 明细分开写
with open(out_genotype, 'w') as f:
    f.write(genotype)
