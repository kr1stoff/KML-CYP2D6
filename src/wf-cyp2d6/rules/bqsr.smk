# rule mark_duplicates:
#     input:
#         "results/mapped/{sample}-{unit}.sorted.bam",
#     output:
#         bam=temp("results/dedup/{sample}-{unit}.bam"),
#         metrics="results/qc/dedup/{sample}-{unit}.metrics.txt",
#     log:
#         "logs/picard/dedup/{sample}-{unit}.log",
#     params:
#         config["params"]["picard"]["MarkDuplicates"],
#     wrapper:
#         "0.74.0/bio/picard/markduplicates"


rule mark_duplicates:
    input:
        rules.bwa_mem.output
    output:
        bam="align/{sample}.dedup.bam",
        metrics="align/{sample}.dedup.metrics.txt",
    benchmark:
        ".log/align/{sample}.bwa_mem.bm"
    log:
        ".log/align/{sample}.bwa_mem.log",
    conda:
        config["conda"]["basic"]
    params:
        config["params"]["picard"]["MarkDuplicates"],
    wrapper:
        "0.74.0/bio/picard/markduplicates"


rule recalibrate_base_qualities:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref="resources/genome.fasta",
        dict="resources/genome.dict",
        known="resources/variation.noiupac.vcf.gz",
        known_idx="resources/variation.noiupac.vcf.gz.tbi",
    output:
        recal_table="results/recal/{sample}-{unit}.grp",
    log:
        "logs/gatk/bqsr/{sample}-{unit}.log",
    params:
        extra=get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"],
    resources:
        mem_mb=1024,
    wrapper:
        "0.74.0/bio/gatk/baserecalibrator"


rule apply_base_quality_recalibration:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref="resources/genome.fasta",
        dict="resources/genome.dict",
        recal_table="results/recal/{sample}-{unit}.grp",
    output:
        bam=protected("results/recal/{sample}-{unit}.bam"),
    log:
        "logs/gatk/apply-bqsr/{sample}-{unit}.log",
    params:
        extra=get_regions_param(),
    resources:
        mem_mb=1024,
    wrapper:
        "0.74.0/bio/gatk/applybqsr"


rule samtools_index:
    input:
        "{prefix}.bam",
    output:
        "{prefix}.bam.bai",
    log:
        "logs/samtools/index/{prefix}.log",
    wrapper:
        "0.74.0/bio/samtools/index"
