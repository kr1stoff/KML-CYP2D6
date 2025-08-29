import click
import logging
from pathlib import Path

from src.kml_cyp2d6.fastq import prepare_fastq_by_samptab
from src.kml_cyp2d6.snakemake import run_snakemake

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')


@click.command()
@click.option('--input-tab', required=True, help='输入样本表, 包含样本名, read1 和 read2 的路径')
@click.option('--output-dir', default='kml-cyp2d6-result', show_default=True, help='输出文件夹')
@click.option('--threads', default=8, type=int, show_default=True, help='线程数')
@click.help_option(help='获取帮助信息')
def main(input_tab, output_dir, threads):
    """KML CYP2D6 分型分析流程"""
    logging.info("开始 KML CYP2D6 分型分析流程")
    output_dir = Path(output_dir).resolve()
    # fastq
    prepare_fastq_by_samptab(output_dir, input_tab, threads)
    # snakemake
    run_snakemake(input_tab, output_dir, threads)
