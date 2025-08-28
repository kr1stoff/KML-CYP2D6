rule analysis_allele_snps:
    input:
        rules.bcftools_view.output,
        config["database"]["pharmvar"],
        config["database"]["annotation"],
    output:
        "allele/{sample}.raw.allele.txt",
        "allele/{sample}.allele.snp.detail.tsv",
        "allele/{sample}.allele.snp.stats.tsv",
    log:
        ".log/allele/{sample}.analysis_allele_snps.log",
    benchmark:
        ".log/allele/{sample}.analysis_allele_snps.bm"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/analysis_allele_snps.py"


rule paste_allele_cnv:
    input:
        rules.parse_cnv.output,
        rules.analysis_allele_snps.output[0],
    output:
        "allele/{sample}.allele.txt",
    log:
        ".log/allele/{sample}.paste_allele_cnv.log",
    benchmark:
        ".log/allele/{sample}.paste_allele_cnv.bm"
    shell:
        "paste {input} > {output} 2> {log}"


rule all_snp_allele:
    input:
        rules.bcftools_view.output,
        config["database"]["pharmvar"],
    output:
        "allele/{sample}.all.snp.detail.tsv",
    log:
        ".log/allele/{sample}.all_snp_allele.log",
    benchmark:
        ".log/allele/{sample}.all_snp_allele.bm"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/all_snp_allele.py"
