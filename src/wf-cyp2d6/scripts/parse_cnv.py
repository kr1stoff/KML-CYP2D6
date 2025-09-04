import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")


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
    ratio = round(cyp2d6['AUTO_RATIO'].mean(), 4)
    cnv_value_count_ser = cyp2d6[cyp2d6['ABBERATION'] != '.']['ABBERATION'].value_counts()
    # 野生型
    if len(cnv_value_count_ser) == 0:
        return 'WT', 0, ratio
    else:
        if cnv_value_count_ser.shape[0] > 1:
            return 'CNV突变类型不唯一', 0, ratio
        else:
            return cnv_value_count_ser.index[0], cnv_value_count_ser.values[0], ratio


convading_res = snakemake.input[0] + '/recal.best.score.totallist.txt'
cnv_type, exon_number, ratio = parse_cnv(convading_res)
# 输出结果
with open(snakemake.output[0], 'w') as f:
    f.write('CNV_TYPE\tEXON_NUMBER\tCNV_RATIO\n')
    f.write(f'{cnv_type}\t{exon_number}\t{ratio}\n')
