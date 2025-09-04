import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")


def call_allele(allele_snp_df) -> tuple:
    """
    调用分型
    :param allele_snp_df: allele snp present 统计表
    :return: 分型结果
    """
    # allele snp present 为空表时, 是 *1/*1
    if len(allele_snp_df) == 0:
        return '1', '1', 'Y'
    # * 输入的应该是排过序的
    allele1, allele2 = allele_snp_df.iloc[0]['ALLELE1'], allele_snp_df.iloc[0]['ALLELE2']
    # 判断最优组合是否非唯一
    top_ratio = allele_snp_df['PRESENT_RATIO'].max()
    top_count = allele_snp_df['PRESENT_COUNT'].max()
    if len(allele_snp_df[(allele_snp_df['PRESENT_RATIO'] == top_ratio)
                         & (allele_snp_df['PRESENT_COUNT'] == top_count)]) > 1:
        return allele1, allele2, 'N'
    else:
        return allele1, allele2, 'Y'


# * 输入
# allele snp present 统计表
allele_snp_df = pd.read_csv(snakemake.input[0], dtype=str, sep='\t')
# 所有 snp 和 allele 信息表
all_snp_df = pd.read_csv(snakemake.input[1], dtype=str, sep='\t')
# NM_000106 注释
nmdf = pd.read_csv(snakemake.input[2], dtype=str).fillna('-')
nmdf.drop('ID', axis=1, inplace=True)

# * 输出
# 分型结果 (不考虑 *5 的情况)
allele1, allele2, unique_info = call_allele(allele_snp_df)
with open(snakemake.output[0], 'w') as f:
    f.write('ALLELE1\tALLELE2\tTOP_UNIQUE\n')
    f.write(f'{allele1}\t{allele2}\t{unique_info}\n')

# 分型 SNP 加 NM_000106.6 注释
outdf = pd.merge(all_snp_df, nmdf, on=['#CHROM', 'POS', 'REF', 'ALT'], how='left')
outdf[['ALLELE', 'ID', 'VARIANT', '#CHROM', 'POS', 'REF', 'ALT', 'GENOTYPE',
       'DETAIL', 'CORE_SNP']].to_csv(snakemake.output[1], index=False, sep='\t')
