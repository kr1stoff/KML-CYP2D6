import pandas as pd
import json
from openpyxl import Workbook
import sys


# * IO
# aldy_res = '/data/mengxf/Project/KML250707-cyp2d6-igtk/work/250707-3samples-anlysis/V12/V12.aldy.txt'
if len(sys.argv) < 2:
    print("Usage: \n\tpython process_aldy.py <aldy-result> <output-excel>")
    exit(1)
aldy_res = sys.argv[1]
output_excel = sys.argv[2]

# * genotype-rsid 参考
with open('/data/mengxf/Project/KML250707-cyp2d6-igtk/database/genotype-rsid-dict.json') as f:
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


# 写入结果文件
wb = Workbook()
ws = wb.active
if ws:
    ws.append([f"Genotype: {genotype}"])
    # 空一行
    ws.append([])
    # 表头
    ws.append(outdf.columns.tolist())
    # 内容
    for _, row in outdf.iterrows():
        ws.append(row.tolist())
    wb.save(output_excel)
else:
    raise Exception('ws is None!')
