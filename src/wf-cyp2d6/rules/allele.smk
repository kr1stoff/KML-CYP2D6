rule call_allele:
    input:
        rules.bcftools_view.output,
        rules.start_with_best_score.output,
        config["database"]["pharmvar"],
        config["database"]["annotation"],
    output:
        "allele/{sample}.allele.txt",
        "allele/{sample}.allele.snp.detail.tsv",
        "allele/{sample}.allele.snp.stats.tsv",
    log:
        ".log/allele/{sample}.call_allele.log",
    benchmark:
        ".log/allele/{sample}.call_allele.bm"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/call_allele.py"


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
