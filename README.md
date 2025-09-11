# KML-CYP2D6

![GitHub followers](https://img.shields.io/github/followers/kr1stoff)
![GitHub Created At](https://img.shields.io/github/created-at/kr1stoff/KML-CYP2D6)
![GitHub commit activity](https://img.shields.io/github/commit-activity/w/kr1stoff/KML-CYP2D6)
![GitHub Release](https://img.shields.io/github/v/release/kr1stoff/KML-CYP2D6)
![GitHub License](https://img.shields.io/github/license/kr1stoff/KML-CYP2D6)

CYP2D6 分析流程

## 开发

1. PharmCAT 分析不了, 需要自建
2. 不使用 PharmVAR, 使用 PharmGKB 对照表生成参考库
3. 新增 CYP2C19 分析流程, 其中
    - C19 默认分型是 \*38 并非 \*1
    - C19\*36 为全基因缺失, 类似 D6\*5
    - C19\*37 为 C19 基因缺失部分外显子

## 命令行

- 分析

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

## 更新

## 注意

1. 不在使用 poetry, 在 snakemake 中总是影响环境, 报 numpy 和 pandas 版本冲突, 直接使用 python 解释器
