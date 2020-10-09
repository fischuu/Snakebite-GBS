# vim: set filetype=sh :

rule picard_markDuplicates:
    """
    Mark duplicate reads (PICARD).
    """
    input:
        "%s/BAM/STAR/{samples}.bam" % (config["project-folder"])
    output:
        bam="%s/BAM/PICARD/{samples}_duplicatesMarked.bam" % (config["project-folder"]),
        txt="%s/BAM/PICARD/{samples}_duplicateMetric.txt" % (config["project-folder"])
    log:
        "%s/logs/PICARD/picard_markDuplicates_{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/PICARD/picard_markDuplicates_{samples}.benchmark.tsv" % (config["project-folder"])
    threads: 20
    params:
      tmpdir=config["local-scratch"]
    conda: "envs/picard.yaml"
    shell:"""
      picard MarkDuplicates I={input} O={output.bam} M={output.txt} TMP_DIR={params.tmpdir} &> {log}
    """
