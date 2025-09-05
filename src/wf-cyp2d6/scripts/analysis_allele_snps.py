import pandas as pd
from pathlib import Path
import sys
import json

sys.stderr = open(snakemake.log[0], "w")


def get_inputs() -> tuple:
    """
    获取输入文件
    :return:
    :vcf: snp vcf 文件
    :pharmgkb_df: 参考表格
    :pharmgkb_allele_snp_count: 每个分型的标记位点数量字典
    :nmdf: 关联 NM_000106 注释表格
    """
    # 输入 SNP 和 CNV 结果文件
    vcf_file = snakemake.input[0]
    # 合并当前样本变异信息和 pharmgkb_df 对应分型
    vcf_df = parse_vcf(vcf_file)
    # 输入 PharmGKB 和 CYP2D6 对照库
    pharmgkb_df = pd.read_csv(snakemake.input[1], dtype=str)
    # 分型核心 SNP 数量
    allele_core_snp_count = pharmgkb_df[~pharmgkb_df['CORE_SNP'].isna()
                                        ][['ALLELE']].groupby('ALLELE').value_counts().to_dict()
    # 每个分型的 SNP 数量
    with open(snakemake.input[2]) as f:
        pharmgkb_allele_snp_count = json.load(f)
    return vcf_df, pharmgkb_df, pharmgkb_allele_snp_count, allele_core_snp_count


def parse_vcf(vcf_file: str) -> pd.DataFrame:
    """解析 VCF 文件, 输出 chr22 的表格"""
    # ! 纯合/杂合 SNP 拷贝数量分型逻辑
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


def get_snp_allele_df(df, pharmgkb_allele_snp_count, allele_core_snp_count, present_threshold=0.6) -> pd.DataFrame:
    """
    统计每个 allele 的 SNP 检出情况, 返回过滤后的 DataFrame
    :param df: 当前样本的 SNP 信息
    :param pharmgkb_allele_snp_count: 每个 allele 的标记位点数量
    :param allele_core_snp_count: 每个 allele 的核心 SNP 数量
    :param present_threshold: 过滤阈值, 当前为 >0.6
    :return: 过滤后的 SNP 统计 DataFrame
    """
    allele_count_ser = df['ALLELE'].value_counts()
    # CYP2D6*1 型 pharmgkb 和 present 都是取 1
    allele_present_list = [
        [a, allele_count_ser[a], pharmgkb_allele_snp_count[a],
            allele_count_ser[a]/pharmgkb_allele_snp_count[a]]
        if a != '1' else [a, allele_count_ser[a], 1, 1]
        for a in allele_count_ser.index
    ]
    snp_allele_df = pd.DataFrame(
        allele_present_list,
        columns=['ALLELE', 'SAMPLE_COUNT', 'PHARMGKB_COUNT', 'PRESENT']
    )
    # ! 核心 SNP 筛选逻辑
    # 样本中检测到 allele 核心 snp 数量统计, 是否与 pharmgkb 理论的核心数量一致
    core_pass_alleles = []
    for allele, count in df[df['CORE_SNP'] == 'CORE'].value_counts('ALLELE').to_dict().items():
        if count == allele_core_snp_count[allele]:
            core_pass_alleles.append(allele)
    # present 大于阈值, 且分型核心 snp 都存在
    return snp_allele_df[(snp_allele_df['ALLELE'].isin(core_pass_alleles)) & (snp_allele_df['PRESENT'] > present_threshold)]


def get_allele1_allele2_present_df(vcf_df, pharmgkb_df, pharmgkb_allele_snp_count, allele_core_snp_count) -> pd.DataFrame:
    key4cols = ['#CHROM', 'POS', 'REF', 'ALT']
    # 所有备选的 allele + snp 统计, 删除 ALLELE 为 NA 的条目, 避免引发报错
    merged_df = pd.merge(vcf_df, pharmgkb_df, how='left', on=key4cols)
    merged_df = merged_df[~merged_df['ALLELE'].isna()]
    # 初始化结果表格
    out_columns = ['ALLELE1', 'ALLELE1_PRESENT_COUNT', 'ALLELE1_PHARMGKB_COUNT',
                   'ALLELE2', 'ALLELE2_PRESENT_COUNT', 'ALLELE2_PHARMGKB_COUNT']
    df = pd.DataFrame(columns=out_columns)
    # 没有变异直接返回空表
    if merged_df.empty:
        return df
    # 获取 allele snp 检出情况, 如果不满足阈值, 直接返回空表
    snp_allele_df = get_snp_allele_df(merged_df, pharmgkb_allele_snp_count, allele_core_snp_count)
    if snp_allele_df.empty:
        return df
    # allele1, allele2 检出 snp 统计
    allele1_allele2_counts = []
    # 先建好 key 避免循环反复创建
    vcf_df['KEY'] = vcf_df['#CHROM'] + '-' + \
        vcf_df['POS'].astype(str) + '-' + vcf_df['REF'] + '-' + vcf_df['ALT']
    # 枚举所有 allele1 + allele2 组合
    for allele1 in snp_allele_df['ALLELE']:
        # * call 第一个 allele
        rec1 = snp_allele_df.loc[snp_allele_df['ALLELE'] == allele1].iloc[0]
        persent1, pharmgkb1 = rec1['SAMPLE_COUNT'], rec1['PHARMGKB_COUNT']
        # 第一个 allele 过的 snp 拷贝 - 1
        remove_sites = set(
            '-'.join(r.tolist())
            for _, r in pharmgkb_df[pharmgkb_df['ALLELE'] == allele1][key4cols].iterrows()
        )
        # * call 第二个 allele
        curdf2 = vcf_df.copy()
        curdf2.loc[curdf2['KEY'].isin(remove_sites), 'COPY'] -= 1
        curdf2 = curdf2[curdf2['COPY'] > 0]
        merged_df2 = pd.merge(curdf2, pharmgkb_df, how='left', on=key4cols)
        merged_df2 = merged_df2[~merged_df2['ALLELE'].isna()]
        if not merged_df2.empty:
            snp_allele_df2 = get_snp_allele_df(
                merged_df2, pharmgkb_allele_snp_count, allele_core_snp_count)
            if not snp_allele_df2.empty:
                for allele2 in snp_allele_df2['ALLELE'].unique():
                    rec2 = snp_allele_df2.loc[snp_allele_df2['ALLELE'] == allele2].iloc[0]
                    present2, pharmgkb2 = rec2['SAMPLE_COUNT'], rec2['PHARMGKB_COUNT']
                    allele1_allele2_counts.append(
                        [allele1, persent1, pharmgkb1, allele2, present2, pharmgkb2])
                # 已添加所有 allele2, 跳过后续
                continue
        # 如果 allele2 没有匹配到分型, 补充默认值
        allele1_allele2_counts.append([allele1, persent1, pharmgkb1, '1', 0, 0])
    # 汇总总表, 统计 allele1 + allele2 的 SNP 检出情况
    return pd.DataFrame(allele1_allele2_counts, columns=out_columns)


def get_all_allele_present_df(indf: pd.DataFrame) -> pd.DataFrame:
    """
    计算每个 allele 检出的比例, 每个 allele 为 0.5, 总 allele 为加和
    :param df: 输入的 DataFrame, 包含 ALLELE1, ALLELE1_PRESENT_COUNT, ALLELE1_PHARMGKB_COUNT,
               ALLELE2, ALLELE2_PRESENT_COUNT, ALLELE2_PHARMGKB_COUNT 列
    :return: 输出的 DataFrame, 包含 ALLELE1_PRESENT_RATIO, ALLELE2_PRESENT_RATIO, PRESENT_RATIO, PRESENT_COUNT 列
    """
    # 不能合并 allele1 和 allele2 算 RATIO, 单个 ALLELE 可能不是 100%
    # 每个 allele 为 0.5, 总 allele 为加和
    df = indf.copy()
    df['ALLELE1_PRESENT_RATIO'] = (df['ALLELE1_PRESENT_COUNT'] / df['ALLELE1_PHARMGKB_COUNT'])
    df['ALLELE1_PRESENT_RATIO'] = df['ALLELE1_PRESENT_RATIO'].apply(lambda x: 1 if x > 1 else x)
    df['ALLELE2_PRESENT_RATIO'] = (df['ALLELE2_PRESENT_COUNT'] / df['ALLELE2_PHARMGKB_COUNT'])
    df['ALLELE2_PRESENT_RATIO'] = df['ALLELE2_PRESENT_RATIO'].apply(lambda x: 1 if x > 1 else x)
    # CYP2D6*1 型 ratio ALLELE2_PRESENT_RATIO 是 0/0
    df.loc[df['ALLELE2_PRESENT_RATIO'].isna(), 'ALLELE2_PRESENT_RATIO'] = 1
    df['PRESENT_RATIO'] = (df['ALLELE1_PRESENT_RATIO'] + df['ALLELE2_PRESENT_RATIO']).round(3)
    df['PRESENT_COUNT'] = df['ALLELE1_PRESENT_COUNT'] + df['ALLELE2_PRESENT_COUNT']
    # 删除 allele1 和 allele2 重复组合
    df['remove_key'] = [frozenset([r['ALLELE1'], r['ALLELE2']]) for _, r in df.iterrows()]
    df = df.drop_duplicates('remove_key').drop(columns=['remove_key'])
    # ! 排序逻辑, present 1.ratio -> 2.count
    allele1_allele2_present_df = df.sort_values(
        by=['PRESENT_RATIO', 'PRESENT_COUNT'], ascending=[False, False]
    ).reset_index(drop=True)
    return allele1_allele2_present_df


def main():
    vcf_df, pharmgkb_df, pharmgkb_allele_snp_count, allele_core_snp_count = get_inputs()
    # 输出 allele1 和 allele2 检出 SNP 统计
    allele1_allele2_present_df = get_allele1_allele2_present_df(
        vcf_df, pharmgkb_df, pharmgkb_allele_snp_count, allele_core_snp_count)
    # 添加总 allele present 信息
    allele_present_df = get_all_allele_present_df(allele1_allele2_present_df)
    # 输出排序一下表头
    allele_present_df[[
        'ALLELE1', 'ALLELE2', 'PRESENT_RATIO', 'PRESENT_COUNT',
        'ALLELE1_PRESENT_COUNT', 'ALLELE1_PHARMGKB_COUNT',
        'ALLELE2_PRESENT_COUNT', 'ALLELE2_PHARMGKB_COUNT',
        'ALLELE1_PRESENT_RATIO', 'ALLELE2_PRESENT_RATIO'
    ]].to_csv(snakemake.output[0], index=False, sep='\t')


if __name__ == "__main__":
    main()
