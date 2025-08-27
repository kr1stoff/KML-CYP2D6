import pandas as pd
from pathlib import Path
import sys

sys.stderr = open(snakemake.log[0], "w")


def get_inputs() -> tuple:
    """
    获取输入文件
    :return:
    :vcf: snp vcf 文件
    :pharmvar_df: 参考表格
    :pharmvar_allele_snp_count: 每个分型的标记位点数量字典
    :nmdf: 关联 NM_000106 注释表格
    """
    # * 输入 SNP 和 CNV 结果文件
    vcf_file = snakemake.input[0]
    # * 输入PharmVAR 和 CYP2D6 对照库
    pharmvar_df = pd.read_csv(snakemake.input[1], dtype=str)
    # 每个分型的 SNP 数量
    pharmvar_allele_snp_count = pharmvar_df['ALLELE'].value_counts().to_dict()
    # * 输入 NM_000106 注释
    nmdf = pd.read_csv(snakemake.input[2], dtype=str).fillna('-')
    nmdf.drop('ID', axis=1, inplace=True)
    return vcf_file, pharmvar_df, pharmvar_allele_snp_count, nmdf


def parse_vcf(vcf_file: str) -> pd.DataFrame:
    """解析 VCF 文件, 输出 chr22 的表格"""
    # 在 SNP 中纯合 COPY 算 2,
    # 第一轮匹配第一个 Allele 然后对应删除使用到的 SNP 一个 COPY;
    # 排除 COPY 为 0 的在进行下一轮匹配, 得到第二个 Allele
    name = Path(vcf_file).stem.split('.')[0]
    # 跳过 VCF ## 的表头
    lines = [line for line in open(vcf_file).readlines() if not line.startswith('##')]
    matrix = [line.strip().split('\t') for line in lines]
    # 读 '#CHROM', 'POS', 'ID', 'REF', 'ALT', 和突变信息列
    df = pd.DataFrame(matrix[1:], columns=matrix[0])[['#CHROM', 'POS', 'ID', 'REF', 'ALT', name]]
    # 样本名统一转成 Detail
    df.rename(columns={name: 'DETAIL'}, inplace=True)
    # 提取 Genotype 和 Copy 信息
    df['GENOTYPE'] = df['DETAIL'].str.split(':').str[0]
    df['COPY'] = df['GENOTYPE'].apply(lambda x: 1 if x == '0/1' else 2)
    # 仅保留 chr22 上的变异, CYP2D6 在 chr22 上
    return df[df['#CHROM'] == 'chr22']


def get_snp_allele_df(df, pharmvar_allele_snp_count, present_threshold=0.6):
    """
    统计每个 allele 的 SNP 检出情况, 返回过滤后的 DataFrame
    :param df: 当前样本的 SNP 信息
    :param pharmvar_allele_snp_count: 每个分型的标记位点数量
    :param present_threshold: 过滤阈值, 当前为 >0.6
    :return: 过滤后的 SNP 统计 DataFrame
    """
    allele_count_ser = df['ALLELE'].value_counts()
    allele_present_list = [
        [a, allele_count_ser[a], pharmvar_allele_snp_count[a],
            allele_count_ser[a]/pharmvar_allele_snp_count[a]]
        for a in allele_count_ser.index
    ]
    snp_allele_df = pd.DataFrame(
        allele_present_list, columns=['ALLELE', 'SAMPLE-COUNT', 'PHARMVAR-COUNT', 'PRESENT']
    )
    snp_allele_df = snp_allele_df.sort_values(
        by=["PRESENT", "SAMPLE-COUNT"], ascending=[False, False]
    )
    return snp_allele_df[snp_allele_df['PRESENT'] > present_threshold]


def get_allele1_allele2_present_df(curdf, pharmvar_df, key4cols, pharmvar_allele_snp_count):
    # 所有备选的 allele + snp 统计
    merged_df = pd.merge(curdf, pharmvar_df, how='left', on=key4cols)
    # 初始化数值, 如果没有变异就是 *1
    allele1_allele2_present_df = pd.DataFrame(
        columns=['ALLELE1', 'ALLELE1-PRESENT-COUNT', 'ALLELE1-PHARMVAR-COUNT',
                 'ALLELE2', 'ALLELE2-PRESENT-COUNT', 'ALLELE2-PHARMVAR-COUNT'])
    # 删除 ALLELE 为 NA 的条目, 避免引发报错
    merged_df = merged_df[~merged_df['ALLELE'].isna()]
    # 不是 *1 的情况
    if not merged_df.empty:
        snp_allele_df = get_snp_allele_df(merged_df, pharmvar_allele_snp_count)
        # 可能经过 present 阈值 dataframe 就空了
        if snp_allele_df.shape[0] > 1:
            # 枚举所有 allele1 + allele2 组合
            # * call 第一个 allele
            allele1_allele2_counts = []
            for allele in snp_allele_df["ALLELE"]:
                # 如果没有完全匹配, 则标记检测到的 SNP和 PHARMVAR 的 SNP 数量
                cur_rec = snp_allele_df.loc[snp_allele_df['ALLELE'] == allele].iloc[0]
                # 输出当前分型检出和 PHARMVAR 的 SNP 数量
                present_snp_count = cur_rec['SAMPLE-COUNT']
                pharmvar_snp_count = cur_rec['PHARMVAR-COUNT']
                # 获取第二个 Allele 的 DataFrame
                curdf2 = curdf.copy()
                curdf2['KEY'] = curdf2['#CHROM'] + '-' + \
                    curdf2['POS'].astype(str) + '-' + curdf2['REF'] + '-' + curdf2['ALT']
                # 第一个 allele 过的 snp 拷贝 - 1
                remove_sites = [
                    '-'.join(r.tolist())
                    for _, r in pharmvar_df[pharmvar_df['ALLELE'] == allele.split('(')[0]][key4cols].iterrows()
                ]
                curdf2.loc[curdf2['KEY'].isin(remove_sites), 'COPY'] -= 1
                curdf2 = curdf2[curdf2['COPY'] > 0]
                # * call 第二个 allele
                # 初始化 present/pharmvar snp count 数值
                allele2, present_snp_count2, pharmvar_snp_count2 = 'CYP2D6_1', 0, 0
                merged_df2 = pd.merge(curdf2, pharmvar_df, how='left', on=key4cols)
                # 删除 ALLELE 为 NA 的条目, 避免引发报错
                merged_df2 = merged_df2[~merged_df2['ALLELE'].isna()]
                # 不是 *1 的情况
                if not merged_df2.empty:
                    snp_allele_df2 = get_snp_allele_df(merged_df2, pharmvar_allele_snp_count)
                    if snp_allele_df2.shape[0] > 1:
                        allele2 = snp_allele_df2.iloc[0]["ALLELE"]
                        # 如果没有完全匹配, 则标记检测到的 SNP和 PHARMVAR 的 SNP 数量
                        top1record2 = snp_allele_df2.iloc[0]
                        present_snp_count2 = top1record2['SAMPLE-COUNT']
                        pharmvar_snp_count2 = top1record2['PHARMVAR-COUNT']
                allele1_allele2_counts.append(
                    [allele, present_snp_count, pharmvar_snp_count, allele2, present_snp_count2, pharmvar_snp_count2])
            # allele1 + allele2 位点检出信息表
            allele1_allele2_present_df = pd.DataFrame(
                allele1_allele2_counts,
                columns=['ALLELE1', 'ALLELE1-PRESENT-COUNT', 'ALLELE1-PHARMVAR-COUNT',
                         'ALLELE2', 'ALLELE2-PRESENT-COUNT', 'ALLELE2-PHARMVAR-COUNT'])
            # allele1 + allele2 位点检出数量
            allele1_allele2_present_df['PRESENT-COUNT'] = allele1_allele2_present_df.apply(
                lambda x: x['ALLELE1-PRESENT-COUNT'] + x['ALLELE2-PRESENT-COUNT'], axis=1)
            # allele1 + allele2 在 PharmVar 中的 SNP 数量
            allele1_allele2_present_df['PHARMVAR-COUNT'] = allele1_allele2_present_df.apply(
                lambda x: x['ALLELE1-PHARMVAR-COUNT'] + x['ALLELE2-PHARMVAR-COUNT'], axis=1)
            # allele1 + allele2 检出比例
            allele1_allele2_present_df['PRESENT-RATIO'] = round(
                allele1_allele2_present_df['PRESENT-COUNT'] / allele1_allele2_present_df['PHARMVAR-COUNT'], 3)
            # 删除 allele1 和 allele2 重复组合
            allele1_allele2_present_df['remove-key'] = [
                set([r['ALLELE1'], r['ALLELE2']])
                for _, r in allele1_allele2_present_df.iterrows()
            ]
            allele1_allele2_present_df = allele1_allele2_present_df.drop_duplicates(
                'remove-key').drop(columns=['remove-key'])
            allele1_allele2_present_df = allele1_allele2_present_df.sort_values(
                by=['PRESENT-RATIO', 'PRESENT-COUNT'], ascending=[False, False]).reset_index(drop=True)
    return allele1_allele2_present_df


def main():
    vcf_file, pharmvar_df, pharmvar_allele_snp_count, nmdf = get_inputs()
    # 合并当前样本变异信息和 pharmvar_df 对应分型
    key4cols = ['#CHROM', 'POS', 'REF', 'ALT']
    curdf = parse_vcf(vcf_file)
    #
    allele1_allele2_present_df = get_allele1_allele2_present_df(
        curdf, pharmvar_df, key4cols, pharmvar_allele_snp_count)

    # * 输出 SNP 统计
    allele1_allele2_present_df.to_csv(snakemake.output[2], index=False, sep='\t')
    # * 输出明确分型的 SNP 明细
    allele1, allele2 = 'CYP2D6_1', 'CYP2D6_1'
    if len(allele1_allele2_present_df) > 0:
        allele1, allele2 = allele1_allele2_present_df.iloc[0]['ALLELE1'], allele1_allele2_present_df.iloc[0]['ALLELE2']
    merged_df = pd.merge(curdf, pharmvar_df, how='left', on=key4cols)
    snp_detail_df = merged_df[
        merged_df['ALLELE'].isin(set([allele1, allele2]))
    ].reset_index(drop=True)
    # 加上 NM_000106.6 注释
    outdf = pd.merge(snp_detail_df, nmdf, on=['#CHROM', 'POS', 'REF', 'ALT'], how='left')
    outdf = outdf[['ALLELE', 'ID', 'VARIANT', '#CHROM', 'POS', 'REF', 'ALT', 'GENOTYPE', 'DETAIL']]
    outdf.to_csv(snakemake.output[1], index=False, sep='\t')
    # * 输出分型结果 (不考虑 *5 的情况)
    with open(snakemake.output[0], 'w') as f:
        f.write('ALLELE1\tALLELE2\n')
        f.write(f'{allele1}\t{allele2}\n')


if __name__ == "__main__":
    main()
