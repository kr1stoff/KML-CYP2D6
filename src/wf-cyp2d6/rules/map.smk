rule bwa_mem:
    input:
        reads=rules.fastp.output.trimmed,
        idx=multiext(
            config["database"]["reference"], ".amb", ".ann", ".bwt", ".pac", ".sa"
        ),
    output:
        "align/{sample}.bam",
    benchmark:
        ".log/align/{sample}.bwa_mem.bm"
    log:
        ".log/align/{sample}.bwa_mem.log",
    conda:
        config["conda"]["basic"]
    params:
        extra=r"-M -Y -R '@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA'",
        sorting="samtools",  # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
        sort_extra="",  # Extra args for samtools/picard.
        tmp_dir="/tmp/",  # Path to temp dir. (optional)
    threads: config["threads"]["high"]
    wrapper:
        f"file:{workflow.basedir}/wrappers/bwa/mem"


rule samtools_index:
    input:
        rules.bwa_mem.output,
    output:
        "align/{sample}.bam.bai",
    benchmark:
        ".log/align/{sample}.samtools_index.bm"
    log:
        ".log/align/{sample}.samtools_index.log",
    params:
        extra="",  # optional params string
    threads: config["threads"]["low"]  # This value - 1 will be sent to -@
    conda:
        config["conda"]["basic"]
    wrapper:
        f"file:{workflow.basedir}/wrappers/samtools/index"


rule samtools_stats:
    input:
        bam=rules.bwa_mem.output,
        bed=rules.bedtools_sort.output,  # Optional input, specify target regions
    output:
        "align/{sample}.bam.target.stat",
    benchmark:
        ".log/align/{sample}.samtools_stats.bm"
    log:
        ".log/align/{sample}.samtools_stats.log",
    conda:
        config["conda"]["basic"]
    threads: config["threads"]["low"]
    wrapper:
        f"file:{workflow.basedir}/wrappers/samtools/stats"


use rule samtools_stats as samtools_stats_all with:
    input:
        rules.bwa_mem.output,
    output:
        "align/{sample}.bam.stat",
    benchmark:
        ".log/align/{sample}.samtools_stats_all.bm"
    log:
        ".log/align/{sample}.samtools_stats_all.log",


rule samtools_depth:
    input:
        bams=[
            rules.bwa_mem.output,
        ],
        bed=".temp/target.sorted.bed",  # optional
    output:
        "align/{sample}.bam.target.depth",
    benchmark:
        ".log/align/{sample}.samtools_depth.bm"
    log:
        ".log/align/{sample}.samtools_depth.log",
    conda:
        config["conda"]["basic"]
    params:
        # optional bed file passed to -b
        extra="",  # optional additional parameters as string
    wrapper:
        f"file:{workflow.basedir}/wrappers/samtools/depth"


rule bam_stats:
    input:
        rules.samtools_stats_all.output,
        rules.samtools_stats.output,
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
        rules.samtools_index.output,
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
