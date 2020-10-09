# vim: set filetype=sh :

rule star_map_reads:
    """
    Map the samples to the genome (STAR).
    """
    input:
        index=config["index"],
        fastq=["%s/FASTQ/CONCATENATED/{samples}_R1_001.merged.fastq.gz" % (config["project-folder"]),
	             "%s/FASTQ/CONCATENATED/{samples}_R2_001.merged.fastq.gz" % (config["project-folder"])]
    output:
        file="%s/BAM/STAR/{samples}.bam" % (config["project-folder"]),
        dir=directory("%s/BAM/STAR/{samples}" % (config["project-folder"]))
    log:
        "%s/logs/STAR/star_map.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/STAR/star_map.{samples}.benchmark.tsv" % (config["project-folder"])
    threads: 20
    conda: "envs/star.yaml"
    shell:"""
        mkdir -p {output.dir};

      	[ ! -d \"{output.dir}\" ] && mkdir {output.dir}

        STAR --genomeDir {input.index} \
            --readFilesIn {input.fastq} \
            --readFilesCommand zcat \
            --outSAMunmapped Within \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMstrandField intronMotif \
            --outSAMattrIHstart 0 \
            --outFilterIntronMotifs RemoveNoncanonical \
            --runThreadN {threads} \
            --outFileNamePrefix {wildcards.samples}_ 2> {log};

        mv {wildcards.samples}_Aligned.sortedByCoord.out.bam {output.file}
        mv {wildcards.samples}_Log.final.out {wildcards.samples}_Log.progress.out {wildcards.samples}_Log.out {wildcards.samples}_SJ.out.tab {output.dir}
    """
