rule fastqc:
    input:
        ".rawdata/{sample}_1.fastq.gz",
        ".rawdata/{sample}_2.fastq.gz",
    output:
        directory("fastqc/{sample}"),
    benchmark:
        ".log/fastqc/{sample}.bm"
    log:
        ".log/fastqc/{sample}.log",
    conda:
        config["conda"]["basic"]
    threads: config["threads"]["low"]
    shell:
        "mkdir {output} && fastqc {input} -o {output} -t {threads} --extract &> {log}"


rule multiqc:
    input:
        expand("fastqc/{sample}", sample=config["samples"]),
    output:
        directory("multiqc"),
    benchmark:
        ".log/multiqc/multiqc.bm"
    log:
        ".log/multiqc/multiqc.log",
    conda:
        config["conda"]["basic"]
    shell:
        "multiqc {input} --outdir {output} 2> {log}"
