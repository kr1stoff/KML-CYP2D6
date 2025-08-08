from pathlib import Path
from src.kml_cyp2d6.fastq import prepare_fastq_by_samptab
from src.kml_cyp2d6.snakemake import create_snakemake_configfile


input_tab = '/data/mengxf/GitHub/KML-CYP2D6/template/input2.tsv'
output_dir = '/data/mengxf/Project/KML250731-cyp2d6-pipeline/results/250804'
threads = 32


output_dir = Path(output_dir).resolve()
# fastq
prepare_fastq_by_samptab(output_dir, input_tab, threads)
# snakemake
create_snakemake_configfile(input_tab, output_dir, threads)
