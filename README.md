# KML-CYP2D6

CYP2D6 分析流程

## 命令行

- poetry 运行

  ```bash
    /home/mengxf/miniforge3/envs/python3.12/bin/python -m src.kml_cyp2d6 \
      --input-tab template/input.tsv \
      --output-dir /data/mengxf/Project/KML250731-cyp2d6-pipeline/results/250731 \
      --threads 32
  ```

- 同步分析结果目录

  ```bash
  rsync -auvP --delete --exclude '**.bam' --exclude '**.gz' \
    /data/mengxf/Project/KML250813-CYP2D6-YANZHENG/results/ \
    /data/share/samba/public/bioinfo/KML250813-CYP2D6-YANZHENG/results/
  ```

- 同步 SAV 文件

    ```bash
    rsync -auvP --delete \
      --include 'Images/***' --include 'InterOp/***' --include 'Thumbnail_Images/***' \
      --include 'RunInfo.xml' --include 'RunParameters.xml' \
      --exclude '*' \
      /data/rawdata/illumina/NEXTseq500/250707_NB501947_0941_AHKKYYBGXW/ \
      /data/share/samba/public/bioinfo/KML250709-lvis-jiance-run1-2/250707_NB501947_0941_AHKKYYBGXW/
    ```

## 开发

1. PharmCAT 分析不了, 需要自建
2. 不使用 PharmVAR, 使用 PharmGKB 对照表生成参考库

## 注意

1. 不在使用 poetry, 在 snakemake 中总是影响环境, 报 numpy 和 pandas 版本冲突, 直接使用 python 解释器
