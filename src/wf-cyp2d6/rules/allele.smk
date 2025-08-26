rule analysis_allele_snps:
    input:
        rules.bcftools_view.output,
        # rules.start_with_best_score.output,
        config["database"]["pharmvar"],
        config["database"]["annotation"],
    output:
        "allele/{sample}.allele.txt",
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


rule summary_allele:
    input:
        expand("allele/{sample}.allele.txt", sample=config["samples"]),
    output:
        "allele/allele.summary.tsv",
        "allele/allele.summary.xlsx",
    log:
        ".log/allele/summary_allele.log",
    benchmark:
        ".log/allele/summary_allele.bm"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/summary_allele.py"


rule csv2xlsx_allele_and_all_snp:
    input:
        rules.analysis_allele_snps.output[1],
        rules.analysis_allele_snps.output[2],
        rules.all_snp_allele.output[0],
    output:
        "allele/{sample}.summary.xlsx",
    log:
        ".log/allele/{sample}.csv2xlsx_allele_and_all_snp.log",
    benchmark:
        ".log/allele/{sample}.csv2xlsx_allele_and_all_snp.bm"
    conda:
        config["conda"]["basic"]
    params:
        "--comment-char '' --tabs --format-numbers",
    shell:
        "csvtk csv2xlsx {params} {input} -o {output} 2> {log}"
