from subprocess import run
import pandas as pd


def check_ref_and_bed(reference: str, bed: str) -> None:
    """
    检查参考基因组和 BED 文件
    :reference: 参考基因组路径
    :bed: BED 文件路径
    """
    # reference id 列表
    # samtools faidx reference
    ref_idx = reference + '.fai'
    ref_ids = pd.read_table(ref_idx, header=None, usecols=[0]).iloc[:, 0].tolist()
    # bed 第一列
    bed_ids = pd.read_table(bed, header=None, usecols=[0]).iloc[:, 0].unique().tolist()
    not_in_ref_ids = [bed_id for bed_id in bed_ids if bed_id not in ref_ids]
    assert len(not_in_ref_ids) == 0, f"BED 文件中存在不在参考基因组中的 ID: {', '.join(not_in_ref_ids)}"
