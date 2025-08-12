from pathlib import Path
from src.kml_cyp2d6.fastq import prepare_fastq_by_samptab
from src.kml_cyp2d6.snakemake import create_snakemake_configfile


input_tab = '/data/mengxf/Project/KML250804-cyp2d6-40baseline/work/250804-input-tsv/input.pgx.tsv'
output_dir = '/data/mengxf/Project/KML250804-cyp2d6-40baseline/results/250811'
threads = 32


output_dir = Path(output_dir).resolve()
# fastq
prepare_fastq_by_samptab(output_dir, input_tab, threads)
# snakemake
create_snakemake_configfile(input_tab, output_dir, threads)
