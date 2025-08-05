rule fastp_pe:
    input:
        sample=[
            ".rawdata/{sample}_1.fastq.gz",
            ".rawdata/{sample}_2.fastq.gz",
        ],
    output:
        trimmed=["fastp/{sample}.1.fastq", "fastp/{sample}.2.fastq"],
        html="fastp/{sample}.html",
        json="fastp/{sample}.json",
    log:
        ".logs/fastp/{sample}.fastp_pe.log",
    benchmark:
        ".log/fastp/{sample}.fastp_pe.bm"
    conda:
        config["conda"]["basic2"]
    threads: config["threads"]["low"]
    params:
        extra="-q 15 -u 40 -l 15 --cut_right --cut_window_size 4 --cut_mean_quality 20 --correction",
    wrapper:
        f"file:{workflow.basedir}/wrappers/fastp"


rule fq_stats_summary:
    input:
        expand("fastp/{sample}.json", sample=config["samples"]),
    output:
        "fastp/fq_summary.tsv",
    benchmark:
        ".log/fastp/fq_all_samples_qc.bm"
    log:
        ".log/fastp/fq_all_samples_qc.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/fq_all_samples_qc.py"
