import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")

# 检出的分型 SNP 位点
all_snp_df = pd.read_csv(snakemake.input[0], sep='\t', usecols=['ID', 'GENOTYPE'])
resent_snp_genotype = all_snp_df.set_index('ID').to_dict()['GENOTYPE']
# 报告位点参考
report_site_df = pd.read_csv(snakemake.input[1])
# 输出导入报告的结果
# SNP
results = []
for it in report_site_df.itertuples():
    if it.rsID in resent_snp_genotype:
        if resent_snp_genotype == '0/1':
            results.append([it.ReportSiteInfo, it.HET])
        else:
            results.append([it.ReportSiteInfo, it.HOM])
    else:
        results.append([it.ReportSiteInfo, it.WT])
# todo CYP2D6*5 CNV
pd.DataFrame(results, columns=['Test', 'Result']).to_csv(snakemake.output[0], index=False, sep='\t')
