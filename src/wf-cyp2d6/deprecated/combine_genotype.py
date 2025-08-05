import sys
from pathlib import Path
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

out_recs = []
for infile in snakemake.input:
    sample = Path(infile).stem.split('.')[0]
    with open(infile, "r") as f:
        out_recs.append([sample, f.read().strip()])

outdf = pd.DataFrame(out_recs, columns=["sample", "genotype"])
outdf.to_csv(snakemake.output[0], index=False, sep="\t")
outdf.to_excel(snakemake.output[1], index=False)
