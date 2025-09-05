# ---------------------------------------
# 对于有简并碱基的分型, 理论的位点数量固定, 简并碱基位点可以任意多
# - 理论的分型位点的数量不包含简并碱基
# - 分型时简并都替换成变异碱基
# 计算 PRESENT RATIO (分型测到的位点数 / 分型中理论的位点数)
# - 如果数值超过 1, 则取 1
# ---------------------------------------
rule analysis_allele_snps:
    input:
        rules.bcftools_view.output,
        config["database"]["pharmgkb"],
        config["database"]["pharmgkb_locus_count"],
    output:
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
        config["database"]["pharmgkb"],
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


rule call_raw_allele:
    input:
        rules.analysis_allele_snps.output,
        rules.all_snp_allele.output,
        config["database"]["annotation"],
    output:
        "allele/{sample}.raw.allele.txt",
        "allele/{sample}.allele.snp.detail.tsv",
    log:
        ".log/allele/{sample}.call_raw_allele.log",
    benchmark:
        ".log/allele/{sample}.call_raw_allele.bm"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/call_raw_allele.py"


rule paste_allele_cnv:
    input:
        rules.parse_cnv.output,
        rules.call_raw_allele.output[0],
    output:
        "allele/{sample}.allele.txt",
    log:
        ".log/allele/{sample}.paste_allele_cnv.log",
    benchmark:
        ".log/allele/{sample}.paste_allele_cnv.bm"
    shell:
        "paste {input} > {output} 2> {log}"
