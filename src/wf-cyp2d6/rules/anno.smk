rule annotate_variants:
    input:
        calls=rules.merge_vcfs.output,  # .vcf, .vcf.gz or .bcf
        cache=config["database"]["vep_cache"],  # can be omitted if fasta and gff are specified
        plugins=config["database"]["vep_plugins"],
        # optionally add reference genome fasta
        # fasta="genome.fasta",
        # fai="genome.fasta.fai", # fasta index
        # gff="annotation.gff",
        # csi="annotation.gff.csi", # tabix index
        # add mandatory aux-files required by some plugins if not present in the VEP plugin directory specified above.
        # aux files must be defined as following: "<plugin> = /path/to/file" where plugin must be in lowercase
        # revel = path/to/revel_scores.tsv.gz
    output:
        calls="anno/{sample}.vcf",  # .vcf, .vcf.gz or .bcf
        stats="anno/{sample}.html",
    params:
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        plugins=["LoFtool"],
        extra=" --everything",
        # optional: extra arguments
    log:
        ".log/anno/{sample}.annotate_variants.log",
    benchmark:
        ".log/anno/{sample}.annotate_variants.bm"
    conda:
        config["conda"]["basic"]
    threads: config["threads"]["low"]
    wrapper:
        f"file:{workflow.basedir}/wrappers/vep/annotate"
 