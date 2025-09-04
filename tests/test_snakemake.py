from src.kml_cyp2d6.snakemake import create_snakemake_configfile, run_snakemake
from pathlib import Path


input_tab = '/data/mengxf/Project/KML250804-cyp2d6-40baseline/input.pgx.tsv'
output_dir = '/data/mengxf/Project/KML250804-cyp2d6-40baseline/results/250811'
# input_tab = '/data/mengxf/Project/KML250813-CYP2D6-YANZHENG/input/input.tsv'
# output_dir = '/data/mengxf/Project/KML250813-CYP2D6-YANZHENG/results/250819'
# input_tab = '/data/mengxf/Project/KML250731-cyp2d6-pipeline/input.tsv'
# output_dir = '/data/mengxf/Project/KML250731-cyp2d6-pipeline/results/250818'
output_dir = Path(output_dir).resolve()
threads = 32

create_snakemake_configfile(input_tab, output_dir, threads)
