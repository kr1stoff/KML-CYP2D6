from pathlib import Path


def get_thread_dict(max_threads: int) -> dict:
    """获线程数字典, 最大/高/低 线程数"""
    return {
        'high': max(1, max_threads // 2),
        'low': max(1, max_threads // 8),
        'max': max_threads
    }


def get_asset_dict() -> dict:
    """获取资源字典"""
    return {
        'genotype_rsid': f'{Path(__file__).parents[2]}/assets/genotype-rsid-dict.json',
    }
