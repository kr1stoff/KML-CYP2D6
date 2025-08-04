from src.kml_cyp2d6.snakemake import create_snakemake_configfile, run_snakemake
from pathlib import Path


input_tab = '/data/mengxf/GitHub/KML-CYP2D6/template/input.tsv'
workdir = '/data/mengxf/Project/KML250731-cyp2d6-pipeline/results/250731'
reference = '/data/mengxf/Database/reference/hg38/hg38.fa'
bed = '/data/mengxf/GitHub/KML-CYP2D6/assets/probeCov.gene.bed'
threads = 32

create_snakemake_configfile(input_tab, Path(workdir), threads, reference, bed)
