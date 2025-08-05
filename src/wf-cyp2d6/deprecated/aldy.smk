rule run_aldy:
    input:
        rules.bwa_mem.output,
    output:
        "aldy/{sample}.aldy.txt",
    benchmark:
        ".log/aldy/{sample}.aldy.bm"
    log:
        ".log/aldy/{sample}.aldy.log",
    conda:
        config["conda"]["aldy"]
    params:
        "--profile pgx1 --gene cyp2d6",
    resources:
        parallel_tasks=2,
    shell:
        "aldy genotype {params} {input} --output {output} --log {log}"


rule parse_aldy:
    input:
        rules.run_aldy.output,
    output:
        "aldy/{sample}.snp.tsv",
        "aldy/{sample}.genotype.txt",
    benchmark:
        ".log/aldy/{sample}.parse_aldy.bm"
    log:
        ".log/aldy/{sample}.parse_aldy.log",
    params:
        config["assets"]["genotype_rsid"],
    conda:
        config["conda"]["python"]
    script:
        "../scripts/process_aldy.py"


rule combine_genotype:
    input:
        expand("aldy/{sample}.genotype.txt", sample=config["samples"]),
    output:
        "aldy/all-sample-genotype.tsv",
        "aldy/all-sample-genotype.xlsx",
    benchmark:
        ".log/aldy/combine_genotype.bm"
    log:
        ".log/aldy/combine_genotype.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/combine_genotype.py"


rule combine_snp_results:
    input:
        expand("aldy/{sample}.snp.tsv", sample=config["samples"]),
    output:
        "aldy/all-sample-snp.xlsx",
    benchmark:
        ".log/aldy/combine_snp_results.bm"
    log:
        ".log/aldy/combine_snp_results.log",
    conda:
        config["conda"]["basic"]
    shell:
        "csvtk csv2xlsx -t {input} -o {output} 2> {log}"
