use rule analysis_allele_snps as cyp2c9_analysis_allele_snps with:
    input:
        rules.bcftools_view.output,
        config["database"]["c9_pharmgkb"],
        config["database"]["c9_pharmgkb_locus_count"],
    output:
        "cyp2c9/{sample}.allele.snp.stats.tsv",
    log:
        ".log/cyp2c9/{sample}.analysis_allele_snps.log",
    benchmark:
        ".log/cyp2c9/{sample}.analysis_allele_snps.bm"
    params:
        chrom="chr10",


use rule all_snp_allele as cyp2c9_all_snp_allele with:
    input:
        rules.bcftools_view.output,
        config["database"]["c9_pharmgkb"],
    output:
        "cyp2c9/{sample}.all.snp.detail.tsv",
    log:
        ".log/cyp2c9/{sample}.all_snp_allele.log",
    benchmark:
        ".log/cyp2c9/{sample}.all_snp_allele.bm"


use rule cyp2c19_call_raw_allele as cyp2c9_call_raw_allele with:
    input:
        rules.cyp2c9_analysis_allele_snps.output,
    output:
        "cyp2c9/{sample}.raw.allele.txt",
    log:
        ".log/cyp2c9/{sample}.call_raw_allele.log",
    benchmark:
        ".log/cyp2c9/{sample}.call_raw_allele.bm"
    params:
        default_allele="1",


use rule csv2xlsx_allele_snp as cyp2c9_csv2xlsx_allele_snp with:
    input:
        rules.cyp2c9_analysis_allele_snps.output,
        rules.cyp2c9_all_snp_allele.output,
    output:
        "cyp2c9/{sample}.summary.xlsx",
    log:
        ".log/cyp2c9/{sample}.csv2xlsx_allele_snp.log",
    benchmark:
        ".log/cyp2c9/{sample}.csv2xlsx_allele_snp.bm"


use rule summary_allele as cyp2c9_summary_allele with:
    input:
        expand("cyp2c9/{sample}.raw.allele.txt", sample=config["samples"]),
    output:
        "cyp2c9/all.allele.summary.tsv",
    log:
        ".log/cyp2c9/summary_allele.log",
    benchmark:
        ".log/cyp2c9/summary_allele.bm"


# CYP2C9 没有 DEL 相关的分型和表型
rule summary_diplotype_phenotype_cyp2c9:
    input:
        rules.cyp2c9_summary_allele.output,
        config["database"]["c9_diplotype_phenotype"],
    output:
        "cyp2c9/all.diplotype_phenotype.tsv",
        "cyp2c9/all.diplotype_phenotype.xlsx",
    log:
        ".log/cyp2c9/summary_diplotype_phenotype.log",
    benchmark:
        ".log/cyp2c9/summary_diplotype_phenotype.bm"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/summary_diplo_pheno_c9.py"
