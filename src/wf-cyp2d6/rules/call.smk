rule haplotype_caller:
    input:
        # single or list of bam files
        bam=rules.apply_base_quality_recalibration.output.bam,
        ref=config["database"]["reference"],
        known=config["database"]["known_site_dbsnp"],  # optional
    output:
        vcf="calls/{sample}.vcf",
        # bam="{sample}.assemb_haplo.bam",
    benchmark:
        ".log/calls/{sample}.haplotype_caller.bm"
    log:
        ".log/calls/{sample}.haplotype_caller.log",
    conda:
        config["conda"]["basic"]
    params:
        extra="--intervals " + config["database"]["bed"],  # optional
        java_opts="",  # optional
    threads: config["threads"]["low"]
    resources:
        mem_mb=16384,
    wrapper:
        f"file:{workflow.basedir}/wrappers/gatk/haplotypecaller"
