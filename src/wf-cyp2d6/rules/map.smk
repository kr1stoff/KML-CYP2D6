rule bwa_mem:
    input:
        ".rawdata/{sample}_1.fastq.gz",
        ".rawdata/{sample}_2.fastq.gz",
    output:
        "align/{sample}.bam",
    benchmark:
        ".log/align/{sample}.bwa_mem.bm"
    log:
        ".log/align/{sample}.bwa_mem.log",
    conda:
        config["conda"]["basic"]
    params:
        bwa="-M -Y -R '@RG\\tID:{sample}\\tSM:{sample}'",
        view="-hbS",
    threads: config["threads"]["high"]
    shell:
        """
        bwa mem -t {threads} {params.bwa} {config[reference]} {input} 2> {log} | \
            samtools view -@ {threads} {params.view} - 2>> {log} | \
            samtools sort -@ {threads} -o {output} - 2>> {log}
        samtools index {output} 2>> {log}
        """


rule samtools_stat:
    input:
        rules.bwa_mem.output,
        rules.bedtools_sort.output,
    output:
        all_stat="align/{sample}.bam.stat",
        target_stat="align/{sample}.bam.target.stat",
    benchmark:
        ".log/align/{sample}.samtools_stat.bm"
    log:
        ".log/align/{sample}.samtools_stat.log",
    conda:
        config["conda"]["basic"]
    threads: config["threads"]["low"]
    shell:
        """
        samtools stat {input[0]} | grep ^SN | cut -f 2- > {output.all_stat} 2> {log}
        samtools stat -t {input[1]} {input[0]} | grep ^SN | cut -f 2- > {output.target_stat} 2>> {log}
        """


rule samtools_depth:
    input:
        rules.bwa_mem.output,
        ".temp/target.sorted.bed",
    output:
        "align/{sample}.bam.target.depth",
    benchmark:
        ".log/align/{sample}.samtools_depth.bm"
    log:
        ".log/align/{sample}.samtools_depth.log",
    conda:
        config["conda"]["basic"]
    shell:
        "samtools depth -b {input[1]} -a {input[0]} -o {output} 2> {log}"


rule bam_stats:
    input:
        rules.samtools_stat.output.all_stat,
        rules.samtools_stat.output.target_stat,
        rules.samtools_depth.output,
    output:
        "align/{sample}.stats.csv",
    benchmark:
        ".log/align/{sample}.bam_stats.bm"
    log:
        ".log/align/{sample}.bam_stats.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/bam_stats.py"


rule bam_stats_summary:
    input:
        expand("align/{sample}.stats.csv", sample=config["samples"]),
    output:
        "align/bam_summary.tsv",
    benchmark:
        ".log/align/bam_stats_summary.bm"
    log:
        ".log/align/bam_stats_summary.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/bam_stats_summary.py"


rule samtools_bedcov:
    input:
        rules.bwa_mem.output,
        rules.bedtools_sort.output,
    output:
        "align/{sample}.bam.target.bedcov",
    benchmark:
        ".log/align/{sample}.samtools_bedcov.bm"
    log:
        ".log/align/{sample}.samtools_bedcov.log",
    conda:
        config["conda"]["basic"]
    shell:
        "samtools bedcov -c {input[1]} {input[0]} > {output} 2> {log}"
