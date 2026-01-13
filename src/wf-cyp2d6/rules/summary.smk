rule all_summary:
    input:
        rules.fq_stats_summary.output,
        rules.bam_stats_summary.output,
    output:
        "upload/panel-qc-summary.tsv",
        "upload/panel-qc-summary.xlsx",
    benchmark:
        ".log/upload/all_summary.bm"
    log:
        ".log/upload/all_summary.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/all_summary.py"
