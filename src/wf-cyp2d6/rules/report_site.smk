rule report_site_info:
    input:
        rules.analysis_allele_snps.output[1],
        config["database"]["report_sites"],
        rules.parse_cnv.output,
    output:
        "report/{sample}.report.site.info.tsv",
    log:
        ".log/report/{sample}.report_site_info.log",
    benchmark:
        ".log/report/{sample}.report_site_info.bm"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/report_site_info.py"


rule csv2xlsx_allele_snp:
    input:
        rules.report_site_info.output,
        rules.analysis_allele_snps.output[2],
        rules.analysis_allele_snps.output[1],
        rules.all_snp_allele.output[0],
    output:
        "report/{sample}.summary.xlsx",
    log:
        ".log/report/{sample}.csv2xlsx_allele_and_all_snp.log",
    benchmark:
        ".log/report/{sample}.csv2xlsx_allele_and_all_snp.bm"
    conda:
        config["conda"]["basic"]
    params:
        "--comment-char '' --tabs --format-numbers",
    shell:
        "csvtk csv2xlsx {params} {input} -o {output} 2> {log}"


rule summary_allele:
    input:
        expand("allele/{sample}.allele.txt", sample=config["samples"]),
    output:
        "report/allele.summary.tsv",
        "report/allele.summary.xlsx",
    log:
        ".log/report/summary_allele.log",
    benchmark:
        ".log/report/summary_allele.bm"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/summary_allele.py"
