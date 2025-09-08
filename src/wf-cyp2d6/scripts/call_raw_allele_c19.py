import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")


def call_allele(allele_snp_df, default_allele) -> tuple:
    """
    调用分型
    :param allele_snp_df: allele snp present 统计表
    :return: 分型结果
    """
    # allele snp present 为空表时, CYP2C19 默认是 *38/*38
    if len(allele_snp_df) == 0:
        return default_allele, default_allele, 'Y'
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
default_allele = snakemake.params.get('default_allele', '38')
# allele snp present 统计表
allele_snp_df = pd.read_csv(snakemake.input[0], dtype=str, sep='\t')

# * 输出
# 分型结果 (不考虑 *5 的情况)
allele1, allele2, unique_info = call_allele(allele_snp_df, default_allele)
with open(snakemake.output[0], 'w') as f:
    f.write('ALLELE1\tALLELE2\tTOP_UNIQUE\n')
    f.write(f'{allele1}\t{allele2}\t{unique_info}\n')
