from src.kml_cyp2d6.snakemake import create_snakemake_configfile, run_snakemake
from pathlib import Path


input_tab = '/data/mengxf/Project/KML250804-cyp2d6-40baseline/work/250804-input-tsv/input.pgx.tsv'
workdir = '/data/mengxf/Project/KML250804-cyp2d6-40baseline/results/250804'
reference = '/data/mengxf/Database/reference/hg38/hg38.fa'
bed = '/data/mengxf/GitHub/KML-CYP2D6/assets/probeCov.gene.bed'
threads = 32

create_snakemake_configfile(input_tab, Path(workdir), threads, reference, bed)
