import os


def get_default_threads() -> int:
    cpu_count = os.cpu_count()
    return cpu_count if cpu_count is not None else 32
