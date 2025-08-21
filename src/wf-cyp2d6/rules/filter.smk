rule gatk_select_snps:
    input:
        vcf=rules.haplotype_caller.output.vcf,
        ref=config["database"]["reference"],
    output:
        vcf="filter/{sample}.snps.vcf",
    benchmark:
        ".log/filter/{sample}.gatk_select_snps.bm"
    log:
        ".log/filter/{sample}.gatk_select_snps.log",
    conda:
        config["conda"]["basic"]
    params:
        extra="--select-type-to-include SNP",  # optional filter arguments, see GATK docs
        java_opts="",  # optional
    resources:
        mem_mb=16384,
    wrapper:
        f"file:{workflow.basedir}/wrappers/gatk/selectvariants"


use rule gatk_select_snps as gatk_select_indels with:
    output:
        vcf="filter/{sample}.indels.vcf",
    benchmark:
        ".log/filter/{sample}.gatk_select_indels.bm"
    log:
        ".log/filter/{sample}.gatk_select_indels.log",
    params:
        extra="--select-type-to-include INDEL",  # optional filter arguments, see GATK docs
        java_opts="",  # optional


rule gatk_filter_snps:
    input:
        vcf=rules.gatk_select_snps.output.vcf,
        ref=config["database"]["reference"],
        intervals=config["database"]["bed"],
    output:
        vcf="filter/{sample}.snps.filtered.vcf",
    log:
        ".log/filter/{sample}.gatk_filter_snps.log",
    benchmark:
        ".log/filter/{sample}.gatk_filter_snps.bm"
    conda:
        config["conda"]["basic"]
    params:
        # MQRankSum < -12.5 || ReadPosRankSum < -8.0, 当前 Panel 过小, 不适合使用该参数
        filters={"snp-hard-filter": "QD < 2.0 || FS > 60.0 || MQ < 40.0"},
        extra="",  # optional arguments, see GATK docs
        java_opts="",  # optional
    resources:
        mem_mb=16384,
    wrapper:
        f"file:{workflow.basedir}/wrappers/gatk/variantfiltration"


use rule gatk_filter_snps as gatk_filter_indels with:
    input:
        vcf=rules.gatk_select_indels.output.vcf,
        ref=config["database"]["reference"],
        intervals=config["database"]["bed"],
    output:
        vcf="filter/{sample}.indels.filtered.vcf",
    log:
        ".log/filter/{sample}.gatk_filter_indels.log",
    benchmark:
        ".log/filter/{sample}.gatk_filter_indels.bm"
    params:
        # ReadPosRankSum < -20.0, 当前 Panel 过小, 不适合使用该参数
        filters={"indel-hard-filter": "QD < 2.0 || FS > 200.0"},
        extra="",  # optional arguments, see GATK docs
        java_opts="",  # optional


rule merge_vcfs:
    input:
        vcfs=[
            rules.gatk_filter_snps.output.vcf,
            rules.gatk_filter_indels.output.vcf,
        ],
    output:
        "filter/{sample}.merged.vcf",
    log:
        ".log/filter/{sample}.merge_vcfs.log",
    benchmark:
        ".log/filter/{sample}.merge_vcfs.bm"
    conda:
        config["conda"]["basic2"]
    params:
        extra="",
        # optional specification of memory usage of the JVM that snakemake will respect with global
        # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
        # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
        # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=16384,
    wrapper:
        f"file:{workflow.basedir}/wrappers/picard/mergevcfs"


rule bcftools_view:
    input:
        rules.merge_vcfs.output,
    output:
        "filter/{sample}.filtered.vcf",
    log:
        ".log/filter/{sample}.bcftools_view.log",
    benchmark:
        ".log/filter/{sample}.bcftools_view.bm"
    conda:
        config["conda"]["basic"]
    params:
        extra="-f 'PASS'",
    wrapper:
        f"file:{workflow.basedir}/wrappers/bcftools/view"
