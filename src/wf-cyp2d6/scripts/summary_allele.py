from pathlib import Path
import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")

allele_files = snakemake.input

dfs = []
for afile in allele_files:
    name = Path(afile).stem.split('.')[0]
    df = pd.read_csv(afile, sep='\t')
    df.insert(0, 'SAMPLE', name)
    dfs.append(df)
outdf = pd.concat(dfs, axis=0)
outdf.sort_values('SAMPLE', inplace=True)

outdf.to_csv(snakemake.output[0], sep='\t', index=False)
