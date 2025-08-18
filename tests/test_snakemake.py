from src.kml_cyp2d6.snakemake import create_snakemake_configfile, run_snakemake
from pathlib import Path


input_tab = '/data/mengxf/GitHub/KML-CYP2D6/template/input2.tsv'
output_dir = '/data/mengxf/Project/KML250731-cyp2d6-pipeline/results/250811'
output_dir = Path(output_dir).resolve()
threads = 32

create_snakemake_configfile(input_tab, output_dir, threads)
