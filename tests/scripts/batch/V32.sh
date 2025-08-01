NAME="V32"
READ1="/data/share/samba/public/bioinfo/KML250707-cyp2d6-igtk/BX9K256021-T2685V2hg38-20250707123000/Raw/V32_R1.fastq.gz"
READ2="/data/share/samba/public/bioinfo/KML250707-cyp2d6-igtk/BX9K256021-T2685V2hg38-20250707123000/Raw/V32_R2.fastq.gz"
OUTDIR="/data/mengxf/Project/KML250707-cyp2d6-igtk/work/250707-3samples-anlysis/V32"
# NAME="V12"
# READ1="/data/share/samba/public/bioinfo/KML250707-cyp2d6-igtk/BX9K256021-T2685V2hg38-20250707123000/Raw/V12_R1.fastq.gz"
# READ2="/data/share/samba/public/bioinfo/KML250707-cyp2d6-igtk/BX9K256021-T2685V2hg38-20250707123000/Raw/V12_R2.fastq.gz"
# OUTDIR="/data/mengxf/Project/KML250707-cyp2d6-igtk/work/250707-3samples-anlysis/V12"

mkdir -p $OUTDIR

# 比对
source /home/mengxf/miniforge3/bin/activate basic

bwa mem -t 32 -M -Y -R '@RG\tID:'$NAME'\tSM:'$NAME \
    /data/mengxf/Database/reference/hg38/hg38.fa \
    $READ1 \
    $READ2 |
    samtools view -@ 32 -hbS - |
    samtools sort -@ 32 -o $OUTDIR/$NAME.sorted.bam -

samtools index $OUTDIR/$NAME.sorted.bam
# 覆盖度 & 深度统计
# bed 区域的 base 数统计
samtools bedcov -c \
    /data/mengxf/Project/KML250707-cyp2d6-igtk/target/probeCov.gene.bed \
    $OUTDIR/$NAME.sorted.bam \
    >$OUTDIR/$NAME.sorted.bam.bedcov
# CYP2D6 & CYP2D7 & CYP2D8P 区域 base 数统计
samtools bedcov -c \
    /data/mengxf/Project/KML250707-cyp2d6-igtk/target/CYP2D6-7-8P.bed \
    $OUTDIR/$NAME.sorted.bam \
    >$OUTDIR/$NAME.sorted.bam.bedcov.cyp2d678p
# bed 区域内的各个位置 base 数统计
samtools depth \
    -b /data/mengxf/Project/KML250707-cyp2d6-igtk/target/probeCov.gene.bed \
    $OUTDIR/$NAME.sorted.bam \
    >$OUTDIR/$NAME.sorted.bam.depth
# 所有比对位置 depth 统计
samtools depth \
    $OUTDIR/$NAME.sorted.bam \
    >$OUTDIR/$NAME.sorted.bam.all.depth
# 比对统计
samtools stat $OUTDIR/$NAME.sorted.bam | grep ^SN | cut -f 2- >$OUTDIR/$NAME.sorted.bam.stat

conda deactivate

# 分型和拷贝数 aldy pgx1
source /home/mengxf/miniforge3/bin/activate aldy

aldy genotype \
    --profile pgx1 \
    --gene cyp2d6 \
    $OUTDIR/$NAME.sorted.bam \
    --output $OUTDIR/$NAME.aldy.txt \
    --log $OUTDIR/$NAME.aldy.log

conda deactivate

# 解析 aldy 结果, 输出分型
/home/mengxf/miniforge3/envs/python3.12/bin/python \
    /data/mengxf/Project/KML250707-cyp2d6-igtk/work/250707-3samples-anlysis/jupyter/process_aldy.py \
    $OUTDIR/$NAME.aldy.txt \
    $OUTDIR/$NAME.genotype.xlsx
