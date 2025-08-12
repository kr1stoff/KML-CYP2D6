rule mark_duplicates:
    input:
        bams=rules.bwa_mem.output,
    output:
        bam="bqsr/{sample}.dedup.bam",
        metrics="bqsr/{sample}.metrics.txt",
    benchmark:
        ".log/bqsr/{sample}.mark_duplicates.bm"
    log:
        ".log/bqsr/{sample}.mark_duplicates.log",
    conda:
        config["conda"]["basic2"]
    params:
        "REMOVE_DUPLICATES=true",
    wrapper:
        f"file:{workflow.basedir}/wrappers/picard/markduplicates"


use rule samtools_index as mark_duplicates_index with:
    input:
        rules.mark_duplicates.output.bam,
    output:
        "bqsr/{sample}.dedup.bam.bai",
    benchmark:
        ".log/bqsr/{sample}.mark_duplicates_index.bm"
    log:
        ".log/bqsr/{sample}.mark_duplicates_index.log",


rule recalibrate_base_qualities:
    input:
        bam=rules.mark_duplicates.output.bam,
        bai=rules.mark_duplicates_index.output,
        ref=config["database"]["reference"],
        dict=config["database"]["dict"],
        known=[
            config["database"]["known_site_1000g"],
            config["database"]["known_site_dbsnp"],
            config["database"]["known_site_mills"],
        ],
        known_idx=[
            config["database"]["known_site_1000g_idx"],
            config["database"]["known_site_dbsnp_idx"],
            config["database"]["known_site_mills_idx"],
        ],
    output:
        recal_table="bqsr/{sample}-recal.grp",
    benchmark:
        ".log/bqsr/{sample}.recalibrate_base_qualities.bm"
    log:
        ".log/bqsr/{sample}.recalibrate_base_qualities.log",
    conda:
        config["conda"]["basic"]
    params:
        extra="--intervals " + config["database"]["bed"],
    wrapper:
        f"file:{workflow.basedir}/wrappers/gatk/baserecalibrator"


rule apply_base_quality_recalibration:
    input:
        bam=rules.mark_duplicates.output.bam,
        bai=rules.mark_duplicates_index.output,
        ref=config["database"]["reference"],
        dict=config["database"]["dict"],
        recal_table=rules.recalibrate_base_qualities.output,
    output:
        bam="bqsr/{sample}/recal.bam",
    benchmark:
        ".log/bqsr/{sample}.apply_base_quality_recalibration.bm"
    log:
        ".log/bqsr/{sample}.apply_base_quality_recalibration.log",
    conda:
        config["conda"]["basic"]
    params:
        extra="--intervals " + config["database"]["bed"],
    wrapper:
        f"file:{workflow.basedir}/wrappers/gatk/applybqsr"


use rule samtools_index as recal_index with:
    input:
        rules.apply_base_quality_recalibration.output.bam,
    output:
        "bqsr/{sample}/recal.bam.bai",
    benchmark:
        ".log/bqsr/{sample}.recal_index.bm"
    log:
        ".log/bqsr/{sample}.recal_index.log",
