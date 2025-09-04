import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")

# 检出的分型 SNP 位点
all_snp_df = pd.read_csv(snakemake.input[0], sep='\t', usecols=['ID', 'GENOTYPE'])
present_snp_genotype = all_snp_df.set_index('ID').to_dict()['GENOTYPE']
# 报告位点参考
report_loci_df = pd.read_csv(snakemake.input[1])
# SNP
results = []
for it in report_loci_df.itertuples():
    if str(it.rsID) in present_snp_genotype:
        if present_snp_genotype[it.rsID] == '0/1':
            results.append([it.ReportSiteInfo, it.HET])
        else:
            results.append([it.ReportSiteInfo, it.HOM])
    else:
        results.append([it.ReportSiteInfo, it.WT])
# CYP2D6*5 CNV
df = pd.read_csv(snakemake.input[2], sep='\t')
cnv = df.iloc[0].to_dict()
if (cnv['CNV_TYPE'] == 'DEL') and (cnv['CNV_RATIO'] < float(snakemake.params.ratio_cutoff_low)):
    results.append(['full gene deletion', 'Positive'])
    results.append(['CYP2D6 copy number', str(round(cnv['CNV_RATIO']*2))])
elif (cnv['CNV_TYPE'] == 'DUP') and (cnv['CNV_RATIO'] > float(snakemake.params.ratio_cutoff_high)):
    results.append(['full gene duplication', 'Negative'])
    results.append(['CYP2D6 copy number', str(round(cnv['CNV_RATIO']*2))])
else:
    results.append(['full gene duplication', 'Negative'])
    results.append(['CYP2D6 copy number', '2'])
# * 输出导入报告的结果
pd.DataFrame(results, columns=['Test', 'Result']).to_csv(snakemake.output[0], index=False, sep='\t')
