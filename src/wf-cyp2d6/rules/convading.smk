"""
# 创建 controls 数据集
# 1. 生成均一化数据集. 使用 BQSR 流程已经去重, 无需 -rmDup 参数
perl /data/mengxf/Software/GitHub/CoNVaDING/CoNVaDING.pl \
    -mode StartWithBam \
    -inputDir all \
    -bed probeCov.gene.bed \
    -useSampleAsControl \
    -controlsDir controls \
    -outputDir .tmp/controls-out
# 2. 生成 TargetQcList
perl /data/mengxf/Software/GitHub/CoNVaDING/CoNVaDING.pl \
    -mode GenerateTargetQcList \
    -controlsDir controls \
    -inputDir controls \
    -outputDir GenerateTargetQcList
"""


rule start_with_bam:
    input:
        bam=rules.apply_base_quality_recalibration.output.bam,
        bai=rules.recal_index.output,
        bed=config["database"]["bed"],
    output:
        directory("convading/{sample}/StartWithBam"),
    benchmark:
        ".log/convading/{sample}.start_with_bam.bm"
    log:
        ".log/convading/{sample}.start_with_bam.log",
    conda:
        config["conda"]["convading"]
    params:
        "-mode StartWithBam -controlsDir " + config["database"]["convading_controls"],
    shell:
        """
        indir=$(dirname {input.bam})
        perl {config[software][convading]} {params} \
            -inputDir $indir \
            -bed {input.bed} \
            -outputDir {output} &> {log}
        """


rule start_with_match_score:
    input:
        rules.start_with_bam.output,
    output:
        directory("convading/{sample}/StartWithMatchScore"),
    benchmark:
        ".log/convading/{sample}.start_with_match_score.bm"
    log:
        ".log/convading/{sample}.start_with_match_score.log",
    conda:
        config["conda"]["convading"]
    params:
        # ! 质控样本数根据实际情况修改
        "-mode StartWithMatchScore -controlSamples 30 -controlsDir "
        + config["database"]["convading_controls"],
    shell:
        """
        perl {config[software][convading]} {params} \
            -inputDir {input} \
            -outputDir {output} &> {log}
        """


rule start_with_best_score:
    input:
        rules.start_with_match_score.output,
    output:
        directory("convading/{sample}/StartWithBestScore"),
    benchmark:
        ".log/convading/{sample}.start_with_best_score.bm"
    log:
        ".log/convading/{sample}.start_with_best_score.log",
    conda:
        config["conda"]["convading"]
    params:
        "-mode StartWithBestScore -ratioCutOffLow 0.65 -ratioCutOffHigh 1.4 -controlsDir "
        + config["database"]["convading_controls"],
    shell:
        """
        perl {config[software][convading]} {params} \
            -inputDir {input} \
            -outputDir {output} &> {log}
        """


rule create_final_list:
    input:
        rules.start_with_best_score.output,
    output:
        directory("convading/{sample}/CreateFinalList"),
    benchmark:
        ".log/convading/{sample}.create_final_list.bm"
    log:
        ".log/convading/{sample}.create_final_list.log",
    conda:
        config["conda"]["convading"]
    params:
        "-mode CreateFinalList -targetQcList "
        + config["database"]["convading_target_qc_list"],
    shell:
        """
        perl {config[software][convading]} {params} \
            -inputDir {input} \
            -outputDir {output} &> {log}
        """


rule parse_cnv:
    input:
        rules.start_with_best_score.output,
    output:
        "convading/{sample}/cnv.txt",
    benchmark:
        ".log/convading/{sample}.parse_cnv.bm"
    log:
        ".log/convading/{sample}.parse_cnv.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/parse_cnv.py"
