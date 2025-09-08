use rule analysis_allele_snps as cyp2c19_analysis_allele_snps with:
    input:
        rules.bcftools_view.output,
        config["database"]["c19_pharmgkb"],
        config["database"]["c19_pharmgkb_locus_count"],
    output:
        "cyp2c19/{sample}.allele.snp.stats.tsv",
    log:
        ".log/cyp2c19/{sample}.analysis_allele_snps.log",
    benchmark:
        ".log/cyp2c19/{sample}.analysis_allele_snps.bm"
    params:
        default_allele="38",
        chrom="chr10"


use rule all_snp_allele as cyp2c19_all_snp_allele with:
    input:
        rules.bcftools_view.output,
        config["database"]["c19_pharmgkb"],
    output:
        "cyp2c19/{sample}.all.snp.detail.tsv",
    log:
        ".log/cyp2c19/{sample}.all_snp_allele.log",
    benchmark:
        ".log/cyp2c19/{sample}.all_snp_allele.bm"


rule cyp2c19_call_raw_allele:
    input:
        rules.cyp2c19_analysis_allele_snps.output,
    output:
        "cyp2c19/{sample}.raw.allele.txt",
    log:
        ".log/cyp2c19/{sample}.call_raw_allele.log",
    benchmark:
        ".log/cyp2c19/{sample}.call_raw_allele.bm"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/call_raw_allele_c19.py"


rule cyp2c19_parse_cnv:
    input:
        rules.start_with_best_score.output,
    output:
        "cyp2c19/{sample}.cnv.txt",
    log:
        ".log/cyp2c19/{sample}.parse_cnv.log",
    benchmark:
        ".log/cyp2c19/{sample}.parse_cnv.bm"
    conda:
        config["conda"]["python"]
    params:
        ratio_cutoff_low=0.65,
    script:
        "../scripts/parse_cnv_c19.py"


use rule paste_allele_cnv as cyp2c19_paste_allele_cnv with:
    input:
        rules.cyp2c19_parse_cnv.output,
        rules.cyp2c19_call_raw_allele.output,
    output:
        "cyp2c19/{sample}.allele.txt",
    log:
        ".log/cyp2c19/{sample}.paste_allele_cnv.log",
    benchmark:
        ".log/cyp2c19/{sample}.paste_allele_cnv.bm"


use rule csv2xlsx_allele_snp as cyp2c19_csv2xlsx_allele_snp with:
    input:
        rules.cyp2c19_analysis_allele_snps.output,
        rules.cyp2c19_all_snp_allele.output,
    output:
        "cyp2c19/{sample}.summary.xlsx",
    log:
        ".log/cyp2c19/{sample}.csv2xlsx_allele_snp.log",
    benchmark:
        ".log/cyp2c19/{sample}.csv2xlsx_allele_snp.bm"


use rule summary_allele as cyp2c19_summary_allele with:
    input:
        expand("cyp2c19/{sample}.allele.txt", sample=config["samples"]),
    output:
        "cyp2c19/all.allele.summary.tsv",
    log:
        ".log/cyp2c19/all.summary_allele.log",
    benchmark:
        ".log/cyp2c19/all.summary_allele.bm"


# TODO *36 全基因缺失, *37 部分外显子缺失
rule summary_diplotype_phenotype_cyp2c19:
    input:
        rules.cyp2c19_summary_allele.output,
        config["database"]["c19_diplotype_phenotype"],
    output:
        "cyp2c19/all.diplotype_phenotype.tsv",
        "cyp2c19/all.diplotype_phenotype.xlsx",
    log:
        ".log/cyp2c19/all.summary_diplotype_phenotype.log",
    benchmark:
        ".log/cyp2c19/all.summary_diplotype_phenotype.bm"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/summary_diplo_pheno_c19.py"
