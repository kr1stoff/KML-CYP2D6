rule bedtools_sort:
    input:
        in_file=config["database"]["bed"],
    output:
        ".temp/target.sorted.bed",
    log:
        ".log/prepare/bedtools_sort.log",
    benchmark:
        ".log/prepare/bedtools_sort.bm"
    conda:
        config["conda"]["basic2"]
    wrapper:
        f"file:{workflow.basedir}/wrappers/bedtools/sort"
