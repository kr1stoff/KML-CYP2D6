def get_thread_dict(max_threads: int) -> dict:
    """获线程数字典, 最大/高/低 线程数"""
    return {
        'high': max(1, max_threads // 2),
        'low': max(1, max_threads // 8),
        'max': max_threads
    }
