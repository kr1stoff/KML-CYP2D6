rule report_loci_info:
    input:
        rules.all_snp_allele.output,
        config["database"]["report_locus"],
        rules.parse_cnv.output,
    output:
        "report/{sample}.report.loci.info.tsv",
    log:
        ".log/report/{sample}.report_loci_info.log",
    benchmark:
        ".log/report/{sample}.report_loci_info.bm"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/report_loci_info.py"


rule csv2xlsx_allele_snp:
    input:
        rules.report_loci_info.output,
        rules.analysis_allele_snps.output,
        rules.call_raw_allele.output[1],
        rules.all_snp_allele.output,
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
