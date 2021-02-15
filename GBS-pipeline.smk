import pandas as pd
from snakemake.utils import validate, min_version
from multiprocessing import cpu_count
import glob
import re
import os

##### GBS-snakemake pipeline #####
##### Daniel Fischer (daniel.fischer@luke.fi)
##### Natural Resources Institute Finland (Luke)
##### This pipeline is build upon the the GBS-SNP-CROP pipeline
##### Version: 0.3.4

##### set minimum snakemake version #####
#min_version("5.1.2")

##### Sample sheets #####

rawsamples = pd.read_table(config["rawsamples"], header=None, sep='\t')[0].tolist()
samples = pd.read_table(config["barcodes"], header=None, sep='\t')[1].tolist()

workdir: config["project-folder"]

wildcard_constraints:
    rawsamples="|".join(rawsamples),
    samples="|".join(samples)

##### Input definitions #####
    
##### run complete pipeline #####

rule all:
    input:
      # OUTPUT: PREPARATION MODULE
#        expand("%s/FASTQ/CONCATENATED/{samples}_R1_001.merged.fastq.gz" % (config["project-folder"]), samples=samples),
#        expand("%s/FASTQ/CONCATENATED/{samples}_R2_001.merged.fastq.gz" % (config["project-folder"]), samples=samples),
#      # QC OF RAW AND CONCATENATED FILES
        "%s/QC/RAW/multiqc_R1/" % (config["project-folder"]),
#        "%s/QC/CONCATENATED/multiqc_R1/" % (config["project-folder"]),
#        "%s/QC/TRIMMED/multiqc_R1/" % (config["project-folder"]),
#      # OUTPUT STEP 2
#        expand("%s/FASTQ/TRIMMED/{samples}.R1.fq.gz" % (config["project-folder"]), samples=samples),
#        expand("%s/FASTQ/TRIMMED/{samples}.R2.fq.gz" % (config["project-folder"]), samples=samples),
#      # OUTPUT STEP 2b
#        expand("%s/FASTQ/SUBSTITUTED/{samples}.R1.fq.gz" % (config["project-folder"]), samples=samples),
#      # OUTPUT STEP 4
#        "%s/FASTQ/TRIMMED/GSC.MR.Genome.fa" % (config["project-folder"]),
#        "%s/BAM/Mockref/mockToRef.sam.flagstat" % (config["project-folder"]),
#        "%s/MPILEUP/mpileup_mockToRef/mockToRef.mpileup" % (config["project-folder"]),
#      # OUTPUT STEP 5
#        expand("%s/FASTQ/TRIMMED/alignments/{samples}.sam.flagstat" % (config["project-folder"]), samples=samples),
#        expand("%s/FASTQ/TRIMMED/alignments/{samples}.sorted.bam" % (config["project-folder"]), samples=samples),
#        expand("%s/FASTQ/TRIMMED/alignments_clusters/{samples}.sam.flagstat" % (config["project-folder"]), samples=samples),
#        expand("%s/FASTQ/TRIMMED/alignments_clusters/{samples}.sorted.bam" % (config["project-folder"]), samples=samples),
#        expand("%s/FASTQ/TRIMMED/alignments_clusters/{samples}.coverage" % (config["project-folder"]), samples=samples),
#        expand("%s/FASTQ/TRIMMED/alignments_reference/{samples}.sorted.bam" % (config["project-folder"]), samples=samples),
#        expand("%s/FASTQ/TRIMMED/alignments_reference/{samples}.sam.flagstat" % (config["project-folder"]), samples=samples),
#        expand("%s/MPILEUP/mpileup_reference/{samples}.mpileup" % (config["project-folder"]), samples=samples),
#      # OUTPUT STEP 6
#        "%s/FASTQ/TRIMMED/GSC.MasterMatrix.txt" % (config["project-folder"]),
#        "%s/MPILEUP/mpileup_reference/GSC.MasterMatrix.txt" % (config["project-folder"]),
#        expand("%s/FASTQ/TRIMMED/{samples}.count.txt" % (config["project-folder"]), samples=samples),
#      # OUTPUT STEP 7
#        "%s/FASTQ/TRIMMED/variants/GSC.GenoMatrix.txt" % (config["project-folder"]),
#        "%s/MPILEUP/mpileup_reference/variants/GSC.GenoMatrix.txt" % (config["project-folder"]),
#      # OUTPUT STEP 8  
#        "%s/FASTQ/TRIMMED/GSC.vcf" % (config["project-folder"]),
#        "%s/FASTQ/TRIMMED/GSC.vcf.fa" % (config["project-folder"]),
#        "%s/MPILEUP/mpileup_reference/GSC.vcf" % (config["project-folder"]),
#        "%s/MPILEUP/mpileup_reference/GSC.vcf.fa" % (config["project-folder"]),
#      # OUTPUT STEP 9
#        "%s/BAM/mockVariantsToReference/mockVariantsToReference.bam" % (config["project-folder"]),
#      # Quality check
#        "%s/finalReport.html" % (config["project-folder"])


### setup report #####

report: "report/workflow.rst"

##### load rules #####
include: "rules/Module0-PreparationsAndIndexing"
include: "rules/Module1-QC"
#include: "rules/Module2-DataPreprocessing"
#include: "rules/Step3-MockReference"
#include: "rules/Step5-AlignReads"
#include: "rules/Step6-ParseMpileup"
#include: "rules/Step7-FilterVariants"
#include: "rules/Step8-CreateVCF"
#include: "rules/Step9-VariantFlanking"
#include: "rules/Step10-Reporting"
