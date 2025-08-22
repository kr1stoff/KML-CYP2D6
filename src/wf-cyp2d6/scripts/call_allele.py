import pandas as pd
from pathlib import Path
import sys

sys.stderr = open(snakemake.log[0], "w")


def get_inputs() -> tuple:
    """
    获取输入文件
    :return:
    :vcf: snp vcf 文件
    :convading: convading 结果文件
    :pharmvar: 参考表格
    :allele_snp_num: 每个分型的标记位点数量字典
    :nmdf: 关联 NM_000106 注释表格
    """
    # 输入 SNP 和 CNV 结果文件
    vcf_file = snakemake.input[0]
    convading_res = snakemake.input[1] + '/recal.best.score.totallist.txt'
    # PharmVAR 和 CYP2D6 对照库, NM_000106 注释
    pharmvar = pd.read_csv(snakemake.input[2], dtype=str)
    # 每个分型的 SNP 数量
    allele_snp_num = pharmvar['ALLELE'].value_counts().to_dict()
    # 关联 NM_000106
    nmdf = pd.read_csv(snakemake.input[3], dtype=str).fillna('-')
    nmdf.drop('ID', axis=1, inplace=True)
    return vcf_file, convading_res, pharmvar, allele_snp_num, nmdf


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


def get_allele(curdf: pd.DataFrame, pharmvar: pd.DataFrame, key4cols: list, allele_snp_num: dict):
    """
    使用当前的变异信息表获取分型, 如果不是 100% 覆盖, 则标记检测到的 SNP 和 PHARMVAR 的 SNP 数量
    :param curdf: 当前样本的 SNP 信息
    :param pharmvar: pharmvar 参考表
    :param key4cols: 合并的 key 列, ['#CHROM', 'POS', 'REF', 'ALT']
    :param allele_snp_num: 每个分型的标记位点数量
    :return:
    :allele: 分型
    :snp_allele_df: 每个分型的 SNP 数量, PHARMVAR 的 SNP 数量, 存在的分型位点比例
    """
    # 所有备选的 allele + snp 统计
    merged_df = pd.merge(curdf, pharmvar, how='left', on=key4cols)
    # 初始化数值, 如果没有变异就是 *1
    allele = 'CYP2D6_1'
    snp_allele_df = pd.DataFrame(columns=['ALLELE', 'SAMPLE-COUNT', 'PHARMVAR-COUNT', 'PRESENT'])
    # * 删除 ALLELE 为 NA 的条目, 避免引发报错
    merged_df = merged_df[~merged_df['ALLELE'].isna()]
    # 不是 *1 的情况
    if not merged_df.empty:
        allele_count_ser = merged_df['ALLELE'].value_counts()
        allele_present_list = [
            [a, allele_count_ser[a], allele_snp_num[a], allele_count_ser[a]/allele_snp_num[a]]
            for a in allele_count_ser.index]
        snp_allele_df = pd.DataFrame(
            allele_present_list, columns=['ALLELE', 'SAMPLE-COUNT', 'PHARMVAR-COUNT', 'PRESENT']
        ).sort_values('PRESENT', ascending=False)
        snp_allele_df = snp_allele_df.sort_values(
            by=["PRESENT", "SAMPLE-COUNT"], ascending=[False, False])
        allele = snp_allele_df.iloc[0]["ALLELE"]
        # 如果没有完全匹配, 则标记检测到的 SNP和 PHARMVAR 的 SNP 数量
        if snp_allele_df.iloc[0]["PRESENT"] < 1:
            top1record = snp_allele_df.iloc[0]
            allele = f"{top1record["ALLELE"]}({top1record['SAMPLE-COUNT']}/{top1record['PHARMVAR-COUNT']})"
    return allele, snp_allele_df


def get_allele2_dataframe(curdf: pd.DataFrame, pharmvar: pd.DataFrame, allele1: str, key4cols: list) -> pd.DataFrame:
    """
    获取第二个 Allele 的 DataFrame
    :param curdf:
    :param pharmvar:
    :param key4cols:
    :return: 第二个 Allele 的 DataFrame, 含位置和 COPY 信息
    """
    df = curdf.copy()
    df['KEY'] = df['#CHROM'] + '-' + df['POS'].astype(str) + '-' + df['REF'] + '-' + df['ALT']
    remove_sites = ['-'.join(r.tolist()) for _, r in pharmvar[pharmvar['ALLELE']
                                                              == allele1.split('(')[0]][key4cols].iterrows()]
    df.loc[df['KEY'].isin(remove_sites), 'COPY'] -= 1
    return df[df['COPY'] > 0]


def parse_cnv(convading_res: str):
    """
    解析 CNV 结果
    :param convading_res:
    :return:
    :cnv_type: CNV 类型
    :exon_number: 变异的 exon 数量
    :ratio: CNV 比例
    """
    df = pd.read_csv(convading_res, sep='\t', usecols=['GENE', 'AUTO_RATIO', 'ABBERATION'])
    cyp2d6 = df[df['GENE'] == 'CYP2D6'].reset_index(drop=True)
    ratio = cyp2d6['AUTO_RATIO'].mean()
    cnv_value_count_ser = cyp2d6[cyp2d6['ABBERATION'] != '.']['ABBERATION'].value_counts()
    # 野生型
    if len(cnv_value_count_ser) == 0:
        return 'WT', 0, ratio
    else:
        if cnv_value_count_ser.shape[0] > 1:
            return 'CNV突变类型不唯一', 0, ratio
        else:
            return cnv_value_count_ser.index[0], cnv_value_count_ser.values[0], ratio


def main():
    vcf_file, convading_res, pharmvar, allele_snp_num, nmdf = get_inputs()
    # 合并当前样本变异信息和 pharmvar 对应分型
    key4cols = ['#CHROM', 'POS', 'REF', 'ALT']
    curdf = parse_vcf(vcf_file)
    # * 第一个 Allele
    allele1, snp_allele_df = get_allele(curdf, pharmvar, key4cols, allele_snp_num)
    # * 第二个 Allele
    curdf2 = get_allele2_dataframe(curdf, pharmvar, allele1, key4cols)
    allele2, _ = get_allele(curdf2, pharmvar, key4cols, allele_snp_num)
    # CNV 结果
    cnv_type, exon_number, ratio = parse_cnv(convading_res)
    # 输出结果
    with open(snakemake.output[0], 'w') as f:
        f.write('CNV-TYPE\tEXON-NUMBER\tCNV-RATIO\tALLELE1\tALLELE2\n')
        f.write(f'{cnv_type}\t{exon_number}\t{ratio}\t{allele1}\t{allele2}\n')
    # 输出 SNP 明细
    merged_df = pd.merge(curdf, pharmvar, how='left', on=key4cols)
    snp_detail_df = merged_df[
        merged_df['ALLELE'].isin(set([allele1.split('(')[0], allele2.split('(')[0]]))
    ].reset_index(drop=True)
    outdf = pd.merge(
        snp_detail_df, nmdf, on=['#CHROM', 'POS', 'REF', 'ALT'], how='left'
    )[['ALLELE', 'ID', 'VARIANT', '#CHROM', 'POS', 'REF', 'ALT', 'GENOTYPE', 'DETAIL']]
    outdf.to_csv(snakemake.output[1], index=False, sep='\t')
    # 输出 SNP 统计
    snp_allele_df.to_csv(snakemake.output[2], index=False, sep='\t')


if __name__ == "__main__":
    main()
