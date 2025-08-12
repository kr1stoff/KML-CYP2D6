cd /data/mengxf/Project/KML250707-cyp2d6-igtk/work/250710-tgcnv/
mkdir V12-out

source /home/mengxf/miniforge3/bin/activate convading

# * bam 文件需要建索引
samtools index -M */*.bam

# * 0. 创建 controls 数据集, 使用 BAM 文件作为输入
# 创建 controls 数据集
perl /data/mengxf/Software/GitHub/CoNVaDING/CoNVaDING.pl \
    -mode StartWithBam \
    -inputDir all \
    -bed probeCov.gene.bed \
    -useSampleAsControl \
    -controlsDir controls \
    -outputDir .tmp/controls-out

# * 1. 使用 BAM 文件作为输入
# BAM 作为输入
perl /data/mengxf/Software/GitHub/CoNVaDING/CoNVaDING.pl \
    -mode StartWithBam \
    -inputDir V12 \
    -bed probeCov.gene.bed \
    -controlsDir controls \
    -outputDir V12-out/StartWithBam

# * 2. 选择最佳 control 样本
# ! controlSamples 后续需要修改, 默认 30, 按照实际数量修改
perl /data/mengxf/Software/GitHub/CoNVaDING/CoNVaDING.pl \
    -mode StartWithMatchScore \
    -controlSamples 30 \
    -inputDir V12-out/StartWithBam \
    -controlsDir controls \
    -outputDir V12-out/StartWithMatchScore

# * 3. CNV 检测
# ! 调整拷贝数阈值, 后续改回默认值. ratioCutOffLow 0.65, ratioCutOffHigh 1.4
perl /data/mengxf/Software/GitHub/CoNVaDING/CoNVaDING.pl \
    -mode StartWithBestScore \
    -ratioCutOffLow 0.65 \
    -ratioCutOffHigh 1.4 \
    -controlsDir controls \
    -inputDir V12-out/StartWithMatchScore \
    -outputDir V12-out/StartWithBestScore

# * 4. QC 列表
perl /data/mengxf/Software/GitHub/CoNVaDING/CoNVaDING.pl \
    -mode GenerateTargetQcList \
    -controlsDir controls \
    -inputDir controls \
    -outputDir V12-out/GenerateTargetQcList

# * 5. 最终列表
perl /data/mengxf/Software/GitHub/CoNVaDING/CoNVaDING.pl \
    -mode CreateFinalList \
    -inputDir V12-out/StartWithBestScore \
    -targetQcList V12-out/GenerateTargetQcList/targetQcList.txt \
    -outputDir V12-out/CreateFinalList

conda deactivate
