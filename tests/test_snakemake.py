from src.kml_cyp2d6.snakemake import create_snakemake_configfile, run_snakemake
from pathlib import Path


input_tab = '/data/mengxf/GitHub/KML-CYP2D6/template/input2.tsv'
reference = '/data/mengxf/Database/reference/hg38/hg38.fa'
bed = '/data/mengxf/GitHub/KML-CYP2D6/assets/probeCov.gene.bed'
output_dir = '/data/mengxf/Project/KML250731-cyp2d6-pipeline/results/250804'
threads = 32
reference = str(Path(reference).resolve())
bed = str(Path(bed).resolve())
output_dir = Path(output_dir).resolve()

create_snakemake_configfile(input_tab, output_dir, threads, reference, bed)
