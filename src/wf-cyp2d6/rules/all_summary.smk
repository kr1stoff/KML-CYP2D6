rule all_summary:
    input:
        rules.fq_stats_summary.output,
        rules.bam_stats_summary.output,
    output:
        "panel-qc-summary.tsv",
        "panel-qc-summary.xlsx",
    benchmark:
        ".log/all_summary.bm"
    log:
        ".log/all_summary.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/all_summary.py"
