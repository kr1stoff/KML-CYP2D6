import yaml
from pathlib import Path
from subprocess import run
import logging

from src.kml_cyp2d6.fastq import get_sample_names_by_samptab
from src.kml_cyp2d6.config import get_thread_dict
from src.config.env import CONDA_ENV_DICT
from src.config.software import ACTIVATE


def create_snakemake_configfile(input_tab: str, workdir: Path, threads: int, reference: str, bed: str) -> str:
    """创建 snakemake 配置文件"""
    logging.info('创建 snakemake 配置文件')
    samples = get_sample_names_by_samptab(input_tab)
    tempdir = workdir.joinpath('.temp')
    tempdir.mkdir(exist_ok=True, parents=True)
    dict_smk = {
        'workdir': str(workdir),
        'samples': samples,
        'threads': get_thread_dict(threads),
        'conda': CONDA_ENV_DICT,
        'reference': reference,
        'bed': bed,
    }
    configfile = f'{workdir}/.temp/snakemake.yaml'
    with open(configfile, 'w') as f:
        yaml.dump(dict_smk, f)
    return configfile


def run_snakemake(input_tab: str, workdir: Path, threads: int, reference: str, bed: str) -> None:
    """
    运行 snakemake 工作流
    :param input_tab:       样本信息表, 样本名/read1/read2
    :param workdir:         分析结果目录
    :param threads:         最大线程数
    :param reference:       参考基因组路径
    :param bed:             BED 文件路径
    :return:
    """
    logging.info('运行 snakemake')
    # 创建 snakemake 配置文件
    configfile = create_snakemake_configfile(input_tab, workdir, threads, reference, bed)
    # 运行 snakemake 流程
    snakefile = Path(__file__).resolve().parents[1].joinpath('wf-cyp2d6/Snakefile')
    logfile = f'{workdir}/.temp/snakemake.log'
    cml = f"""
    source {ACTIVATE} snakemake
    # use-conda
    snakemake -c {threads} --use-conda -s {snakefile} --configfile {configfile} --ignore-incomplete --scheduler greedy
    conda deactivate
    """
    logging.debug(cml)
    proc = run(cml, shell=True, executable='/bin/bash', capture_output=True, encoding='utf-8')
    # 输出出来这段日志
    with open(logfile, 'w') as f:
        f.write(f'[STDOUT]\n{proc.stdout}\n[STDERR]\n{proc.stderr}')
