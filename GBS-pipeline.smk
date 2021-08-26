import pandas as pd
from snakemake.utils import validate, min_version
from multiprocessing import cpu_count
import glob
import re
import os

##### GBS-snakemake pipeline #####
##### Daniel Fischer (daniel.fischer@luke.fi)
##### Natural Resources Institute Finland (Luke)
##### This pipeline is build upon the the GBS-SNP-CROP pipeline:
##### https://github.com/halelab/GBS-SNP-CROP
##### Version: 0.8.7
version = "0.8.7"

##### set minimum snakemake version #####
min_version("6.0")

##### Sample sheets #####

rawsamples = pd.read_table(config["rawsamples"], header=None, sep='\t')[0].tolist()
samples = pd.read_table(config["barcodes"], header=None, sep='\t')[1].tolist()

workdir: config["project-folder"]

wildcard_constraints:
    rawsamples="|".join(rawsamples),
    samples="|".join(samples)

##### Complete the input configuration
config["genome-bwa-index"] = config["genome"]+".bwt"
config["mockref-bwa-index"] = config["mockreference"]+".bwt"
config["genome-star-index"] = config["project-folder"]+"/references/STAR2.7.3a"
config["report-script"] = config["pipeline-folder"]+"/scripts/workflow-report.Rmd"
config["mockeval-script"] = config["pipeline-folder"]+"/scripts/mockeval-report.Rmd"
config["refinement-script"] = config["pipeline-folder"]+"/scripts/refineMockReference.R"
config["adapter"]=config["pipeline-folder"]+"/adapter.fa"

##### Singularity container #####
config["singularity"] = {}
config["singularity"]["star"] = "docker://fischuu/star:2.7.3a-0.1"
config["singularity"]["gbs"] = "docker://fischuu/gbs:0.2"
config["singularity"]["cutadapt"] = "docker://fischuu/cutadapt:2.8-0.3"
config["singularity"]["minimap2"] = "docker://fischuu/minimap2:2.17-0.2"
config["singularity"]["r-gbs"] = "docker://fischuu/r-gbs:3.6.3-0.2"

##### Print the welcome screen #####
print("#################################################################################")
print("##### Welcome to the GBS pipeline")
print("##### version: "+version)
print("#####")
print("##### Pipeline configuration")
print("##### --------------------------------")
print("##### project-folder  : "+config["project-folder"])
print("##### pipeline-folder : "+config["pipeline-folder"])
print("##### report-script   : "+config["report-script"])
print("##### pipeline-config (NOT NECESSARILY THE USED ONE!!!): "+config["pipeline-config"])
print("#####")
print("##### Singularity configuration")
print("##### --------------------------------")
print("##### star     : "+config["singularity"]["star"])
print("##### gbs      : "+config["singularity"]["gbs"])
print("##### cutadapt : "+config["singularity"]["cutadapt"])
print("##### minimap2 : "+config["singularity"]["minimap2"])
print("##### r-gbs    : "+config["singularity"]["r-gbs"])
print("#####")
print("##### Runtime-configurations")
print("##### --------------------------------")
print("##### genome         : "+ config["genome"])
print("##### existing mock  : "+ config["mockreference"])
print("##### barcodes-file  : "+ config["barcodes"])
print("##### rawsample file : "+ config["rawsamples"])
print("#####")
print("##### Derived runtime parameters")
print("##### --------------------------------")
print("##### BWA-Genome index    : "+config["genome-bwa-index"])
print("##### STAR-Genome index   : "+config["genome-star-index"])
print("##### Existing Mock index : "+config["mockref-bwa-index"])
print("##### Adapter file        : "+ config["adapter"])
print("#################################################################################")

##### Define conditional input/outputs #####
conditionalOut = list()
if config["mockreference"] != "":
        conditionalOut.append("%s/VCF/FinalSetVariants_existingMock.vcf" % (config["project-folder"]))

    
##### run complete pipeline #####

rule all:
    input:
        conditionalOut,
      # QC OF RAW AND CONCATENATED FILES
        "%s/QC/RAW/multiqc_R1/" % (config["project-folder"]),
        "%s/QC/CONCATENATED/multiqc_R1/" % (config["project-folder"]),
        "%s/QC/TRIMMED/multiqc_R1/" % (config["project-folder"]),
      # OUTPUT STEP 4
        "%s/FASTQ/TRIMMED/GSC.MR.Genome.fa" % (config["project-folder"]),
        "%s/BAM/Mockref/mockToRef.sam.flagstat" % (config["project-folder"]),
        "%s/MPILEUP/mpileup_mockToRef/mockToRef.mpileup" % (config["project-folder"]),
      # OUTPUT STEP 5
        expand("%s/FASTQ/TRIMMED/alignments/{samples}.sam.flagstat" % (config["project-folder"]), samples=samples),
        expand("%s/FASTQ/TRIMMED/alignments/{samples}.sorted.bam" % (config["project-folder"]), samples=samples),
        expand("%s/FASTQ/TRIMMED/alignments_clusters/{samples}.sam.flagstat" % (config["project-folder"]), samples=samples),
        expand("%s/FASTQ/TRIMMED/alignments_clusters/{samples}.sorted.bam" % (config["project-folder"]), samples=samples),
        expand("%s/FASTQ/TRIMMED/alignments_clusters/{samples}.coverage" % (config["project-folder"]), samples=samples),
        expand("%s/FASTQ/TRIMMED/alignments_reference/{samples}.sorted.bam" % (config["project-folder"]), samples=samples),
        expand("%s/FASTQ/TRIMMED/alignments_reference/{samples}.sam.flagstat" % (config["project-folder"]), samples=samples),
        expand("%s/MPILEUP/mpileup_reference/{samples}.mpileup" % (config["project-folder"]), samples=samples),
      # OUTPUT STEP 6
        "%s/FASTQ/TRIMMED/GSC.MasterMatrix.txt" % (config["project-folder"]),
        "%s/MPILEUP/mpileup_reference/GSC.MasterMatrix.txt" % (config["project-folder"]),
        expand("%s/FASTQ/TRIMMED/{samples}.count.txt" % (config["project-folder"]), samples=samples),
      # OUTPUT STEP 7
        "%s/FASTQ/TRIMMED/variants/GSC.GenoMatrix.txt" % (config["project-folder"]),
        "%s/MPILEUP/mpileup_reference/variants/GSC.GenoMatrix.txt" % (config["project-folder"]),
      # OUTPUT STEP 8  
        "%s/FASTQ/TRIMMED/GSC.vcf" % (config["project-folder"]),
        "%s/FASTQ/TRIMMED/GSC.vcf.fa" % (config["project-folder"]),
        "%s/VCF/FinalSetVariants_referenceGenome.vcf" % (config["project-folder"]),
        "%s/MPILEUP/mpileup_reference/GSC.vcf.fa" % (config["project-folder"]),
      # OUTPUT STEP 9
        "%s/BAM/mockVariantsToReference/mockVariantsToReference.bam" % (config["project-folder"]),
      # Quality check
        expand("%s/BAM/alignments_finalMock/{samples}.sam.flagstat" % (config["project-folder"]), samples=samples),
        "%s/MockReference/MockReference.fa" % (config["project-folder"]),
        "%s/VCF/FinalSetVariants_finalMock.vcf" % (config["project-folder"]),
        "%s/finalReport.html" % (config["project-folder"])

rule QC:
    input:
        "%s/QC/RAW/multiqc_R1/" % (config["project-folder"]),
        "%s/QC/CONCATENATED/multiqc_R1/" % (config["project-folder"]),
        "%s/QC/TRIMMED/multiqc_R1/" % (config["project-folder"]),

rule MockEval:
    input:
        "%s/mockEvalReport.html" % (config["project-folder"])

rule FinalReport:
    input:
        "%s/finalReport.html" % (config["project-folder"])

rule MockRefVCF:
    input:
        "%s/VCF/FinalSetVariants_finalMock.vcf" % (config["project-folder"])

### setup report #####

report: "report/workflow.rst"

##### load rules #####
include: "rules/Module0-PreparationsAndIndexing"
include: "rules/Module1-QC"
include: "rules/Module2-DataPreprocessing"
include: "rules/Module3-MockReference"
include: "rules/Module4-ReadAlignment"
include: "rules/Module5-CallVariants"
include: "rules/Module6-PostProcessing"
include: "rules/Module7-Reporting"
include: "rules/Module8-CallNewData"
