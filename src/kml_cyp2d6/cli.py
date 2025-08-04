import click
import logging
from pathlib import Path

from src.kml_cyp2d6.fastq import prepare_fastq_by_samptab
from src.kml_cyp2d6.check import check_ref_and_bed
from src.kml_cyp2d6.snakemake import run_snakemake

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')


@click.command()
@click.option('--input-tab', required=True, help='输入样本表, 包含样本名, read1 和 read2 的路径')
@click.option('--output-dir', default='kml-cyp2d6-result', show_default=True, help='输出文件夹')
@click.option('--reference', required=True, help='参考基因组路径, 提前创建(BWA/Samtools)索引')
@click.option('--bed', required=True, help='BED 文件路径, 用于指定试剂盒靶区域, 基因组版本需与参考基因组一致')
@click.option('--threads', default=8, type=int, show_default=True, help='线程数')
@click.help_option(help='获取帮助信息')
def main(input_tab, output_dir, reference, bed, threads):
    """KML Panel 快速质量评估流程"""
    logging.info("开始 KML Panel 快速质量评估流程")
    reference = str(Path(reference).resolve())
    bed = str(Path(bed).resolve())
    output_dir = Path(output_dir).resolve()
    # fastq
    prepare_fastq_by_samptab(output_dir, input_tab, threads)
    # 检查 reference 和 bed
    check_ref_and_bed(reference, bed)
    # snakemake
    run_snakemake(input_tab, output_dir, threads, reference, bed)
