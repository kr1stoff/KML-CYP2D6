import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")


df = pd.read_csv(snakemake.input[0] + '/recal.best.score.totallist.txt',
                 sep='\t', usecols=['GENE', 'AUTO_RATIO', 'ABBERATION'])
cyp2c19 = df[df['GENE'] == "CYP2C19"].reset_index(drop=True)
# 统计 (DEL) 和 (AUTO_RATIO < 0.65) 的外显子数
del_exon_count = len(cyp2c19[(cyp2c19['ABBERATION'] == 'DEL') & (cyp2c19['AUTO_RATIO'] < 0.65)])

with open(snakemake.output[0], 'w') as f:
    f.write('DELETE_EXON_COUNT\n')
    f.write(f'{del_exon_count}\n')
