rule bedtools_sort:
    output:
        ".temp/target.sorted.bed",
    benchmark:
        ".log/prepare/bedtools_sort.bm"
    log:
        ".log/prepare/bedtools_sort.log",
    conda:
        config["conda"]["basic2"]
    shell:
        "bedtools sort -i {config[bed]} > {output} 2> {log}"
