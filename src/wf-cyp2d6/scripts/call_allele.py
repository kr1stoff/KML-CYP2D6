import pandas as pd
from pathlib import Path
import sys

sys.stderr = open(snakemake.log[0], "w")

# * 输入
vcf = snakemake.input[0]
convading_res = snakemake.input[1] + '/recal.best.score.totallist.txt'

# * PharmVAR 和 CYP2D6 对照库, NM_000106 注释
pharmvar = pd.read_csv(snakemake.input[2], dtype=str)
# 每个分型的 SNP 数量
allele_snp_num = pharmvar['ALLELE'].value_counts().to_dict()
# 关联 NM_000106
nmdf = pd.read_csv(snakemake.input[3], dtype=str).fillna('-')
nmdf.drop('ID', axis=1, inplace=True)

# * 分型结果
# 在 SNP 中纯合 COPY 算 2,
# 第一轮匹配第一个 Allele 然后对应删除使用到的 SNP 一个 COPY;
# 排除 COPY 为 0 的在进行下一轮匹配, 得到第二个 Allele
name = Path(vcf).stem.split('.')[0]
# 跳过 VCF ## 的表头
lines = [line for line in open(vcf).readlines() if not line.startswith('##')]
matrix = [line.strip().split('\t') for line in lines]
# 读 '#CHROM', 'POS', 'ID', 'REF', 'ALT', 和突变信息列
curdf = pd.DataFrame(matrix[1:], columns=matrix[0])[['#CHROM', 'POS', 'ID', 'REF', 'ALT', name]]
# 样本名统一转成 Detail
curdf.rename(columns={name: 'DETAIL'}, inplace=True)
# 提取 Genotype 和 Copy 信息
curdf['GENOTYPE'] = curdf['DETAIL'].str.split(':').str[0]
curdf['COPY'] = curdf['GENOTYPE'].apply(lambda x: 1 if x == '0/1' else 2)
# 仅保留 chr22 上的变异, CYP2D6 在 chr22 上
curdf = curdf[curdf['#CHROM'] == 'chr22']

# 合并当前样本变异信息和 pharmvar 对应分型
key4cols = ['#CHROM', 'POS', 'REF', 'ALT']
merged_df = pd.merge(curdf, pharmvar, how='left', on=key4cols)


def get_allele(curdf):
    """使用当前的变异信息表获取分型"""
    # 所有备选的 allele + snp 统计
    allele_present_list = []
    merged_df = pd.merge(curdf, pharmvar, how='left', on=key4cols)
    # 初始化数值, 如果没有变异就是 *1
    allele = 'CYP2D6_1'
    snp_allele_df = pd.DataFrame(columns=['ALLELE', 'Sample_count', 'PharmVar_count', 'Present'])
    # 不是 *1 的情况
    if merged_df.shape[0] != 0:
        allele_count_ser = merged_df['ALLELE'].value_counts()
        for allele in allele_count_ser.index:
            allele_present_list.append(
                [allele, allele_count_ser[allele],
                 allele_snp_num[allele], allele_count_ser[allele]/allele_snp_num[allele]])
        snp_allele_df = pd.DataFrame(
            allele_present_list, columns=['ALLELE', 'Sample_count', 'PharmVar_count', 'Present']
        ).sort_values('Present', ascending=False)
        snp_allele_df = snp_allele_df.sort_values(
            by=["Present", "Sample_count"], ascending=[False, False])
        allele = snp_allele_df.iloc[0]["ALLELE"]
    return allele, snp_allele_df


# 第一个 Allele
allele1, snp_allele_df = get_allele(curdf)
# ! 分型位点不是 100% 出现就写"不是 100% 覆盖分型位点"
if (allele1 != 'CYP2D6_1') and (snp_allele_df.iloc[0]["Present"] < 1):
    allele2 = "不是 100% 覆盖分型位点"
else:
    # 第二个 Allele
    # 挑出 allele1 的 snp 从 curdf 中删除一个 copy
    allele1_snp_df = pharmvar[pharmvar['ALLELE'] == allele1][key4cols]
    allele1_snp_idx = allele1_snp_df.set_index(key4cols).index
    curdf2 = curdf.set_index(key4cols).copy()
    curdf2.loc[allele1_snp_idx, 'COPY'] -= 1
    # COPY 为 0 的跳过
    curdf2 = curdf2[curdf2['COPY'] > 0]
    allele2, _ = get_allele(curdf2)

# * 拷贝数结果
# data/mengxf/Project/KML250731-cyp2d6-pipeline/results/250818/convading/PGXt1/StartWithBestScore/recal.best.score.totallist.txt
df = pd.read_csv(convading_res, sep='\t', usecols=['GENE', 'AUTO_RATIO', 'ABBERATION'])
cyp2d6 = df[df['GENE'] == 'CYP2D6'].reset_index(drop=True)
ratio = cyp2d6['AUTO_RATIO'].mean()
cnv_value_count_ser = cyp2d6[cyp2d6['ABBERATION'] != '.']['ABBERATION'].value_counts()
# 野生型
if len(cnv_value_count_ser) == 0:
    cnv_type, exon_number = 'WT', 0
else:
    if cnv_value_count_ser.shape[0] > 1:
        cnv_type, exon_number = 'CNV突变类型不唯一', 0
    else:
        cnv_type, exon_number = cnv_value_count_ser.index[0], cnv_value_count_ser.values[0]

# * 输出
with open(snakemake.output[0], 'w') as f:
    f.write('CNV-TYPE\tEXON-NUMBER\tRATIO\tALLELE1\tALLELE2\n')
    f.write(f'{cnv_type}\t{exon_number}\t{ratio}\t{allele1}\t{allele2}\n')
# 输出 SNP 明细
snp_detail_df = merged_df[merged_df['ALLELE'].isin(set([allele1, allele2]))].reset_index(drop=True)
outdf = pd.merge(snp_detail_df, nmdf, on=['#CHROM', 'POS', 'REF', 'ALT'], how='left')[
    ['ALLELE', 'ID', 'VARIANT', '#CHROM', 'POS', 'REF', 'ALT', 'GENOTYPE', 'DETAIL']]
outdf.to_csv(snakemake.output[1], index=False, sep='\t')
# 输出 SNP 统计
snp_allele_df.to_csv(snakemake.output[2], index=False, sep='\t')
