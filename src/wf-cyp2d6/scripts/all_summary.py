import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")

dfqc = pd.read_csv(snakemake.input[0], sep="\t")
dfbam = pd.read_csv(snakemake.input[1], sep="\t")
dfmerged = pd.merge(dfqc, dfbam, on='Sample', how='left')
dfmerged.to_csv(snakemake.output[0], index=False, sep='\t')
dfmerged.to_excel(snakemake.output[1], index=False)
