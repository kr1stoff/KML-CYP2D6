rule fastp_pe:
    input:
        ".rawdata/{sample}_1.fastq.gz",
        ".rawdata/{sample}_2.fastq.gz",
    output:
        j="fastp/{sample}.json",
        h="fastp/{sample}.html",
        o="fastp/{sample}.1.fastq.gz",
        O="fastp/{sample}.2.fastq.gz",
    benchmark:
        ".log/fastp/{sample}.fastp_pe.bm"
    log:
        ".log/fastp/{sample}.fastp_pe.log",
    conda:
        config["conda"]["basic2"]
    threads: config["threads"]["low"]
    shell:
        "fastp -w {threads} -j {output.j} -h {output.h} -o {output.o} -O {output.O} -i {input[0]} -I {input[1]} &> {log}"


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
