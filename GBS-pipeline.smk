import pandas as pd
from snakemake.utils import validate, min_version
from multiprocessing import cpu_count
import glob
import re
import os
import sys

##### GBS-snakemake pipeline #####
##### Daniel Fischer (daniel.fischer@luke.fi)
##### Natural Resources Institute Finland (Luke)
##### This pipeline is build upon the the GBS-SNP-CROP pipeline:
##### https://github.com/halelab/GBS-SNP-CROP
##### Version: 0.16.1
version = "0.16.1"

##### set minimum snakemake version #####
min_version("6.0")

##### Sample sheets #####

##### load config and sample sheets #####

samplesheet = pd.read_table(config["samplesheet-file"]).set_index("rawsample", drop=False)
rawsamples=list(samplesheet.rawsample)
samples=list(set(list(samplesheet.sample_name)))
lane=list(samplesheet.lane)

#### CONTINUE FROM HERE TO ADD PIPE CONFIG ONTO THE FILE
#if '--configfile' in sys.argv:
#    i = sys.argv.index('--configfile')
#    config["pipeline-config"] = sys.argv[i + 1]
#else:
#    config["pipeline-config"] = ""
#
#if '--cluster-config' in sys.argv:
#    i = sys.argv.index('--cluster-config')
#    config["server-config"] = sys.argv[i + 1]
#else:
#    config["server-config"] = ""

workdir: config["project-folder"]

wildcard_constraints:
    rawsamples="|".join(rawsamples),
    samples="|".join(samples)

##### input function definitions ######

def get_lane(wildcards):
    output = samplesheet.loc[wildcards.rawsamples][["lane"]]
    return output.tolist()

def get_sample(wildcards):
    output = samplesheet.loc[wildcards.rawsamples][["sample_name"]]
    return output.tolist()

def get_raw_input_fastqs(wildcards):
    reads = samplesheet.loc[wildcards.rawsamples][["read1", "read2"]]
    path = config["rawdata-folder"]
    output = [path + x for x in reads]
    return output

def get_raw_input_read1(wildcards):
    reads = samplesheet.loc[wildcards.rawsamples][["read1"]]
    path = config["rawdata-folder"]
    output = [path + "/" + x for x in reads]
    return output

def get_raw_input_read2(wildcards):
    reads = samplesheet.loc[wildcards.rawsamples][["read2"]]
    path = config["rawdata-folder"]
    output = [path + "/" + x for x in reads]
    return output
    
def get_fastq_for_concatenating_read1(wildcards):
    r1 = samplesheet.loc[samplesheet["sample_name"] == wildcards.samples]["read1"]
    path = config["rawdata-folder"] + "/"
    output = [path + x for x in r1]
    return output   

def get_fastq_for_concatenating_read2(wildcards):
    r1 = samplesheet.loc[samplesheet["sample_name"] == wildcards.samples]["read2"]
    path = config["rawdata-folder"] + "/"
    output = [path + x for x in r1]
    return output   

##### Complete the input configuration
config["genome-bwa-index"] = config["genome"]+".bwt"
config["mockref-bwa-index"] = config["mockreference"]+".bwt"
config["genome-star-index"] = config["project-folder"]+"/References/STAR2.7.5a"    # Change here to the path from the reference genome!!!!!!!!!!!!!!!!!!
config["barcodes-script"] = config["pipeline-folder"]+"/scripts/prepareBarcodes.R"
config["report-script"] = config["pipeline-folder"]+"/scripts/Final-report.R"
config["qc-script"] = config["pipeline-folder"]+"/scripts/QC-report.R"
config["mockeval-script"] = config["pipeline-folder"]+"/scripts/mockeval-report.Rmd"
config["refinement-script"] = config["pipeline-folder"]+"/scripts/refineMockReference.R"
config["adapter"]=config["pipeline-folder"]+"/adapter.fa"
config["barcodes-file"] = config["project-folder"]+"/barcodesID.txt"

##### Singularity container #####
config["singularity"] = {}
config["singularity"]["star"] = "docker://fischuu/star:2.7.5a"
config["singularity"]["gbs"] = "docker://fischuu/gbs:0.2"
config["singularity"]["cutadapt"] = "docker://fischuu/cutadapt:2.8-0.3"
config["singularity"]["minimap2"] = "docker://fischuu/minimap2:2.17-0.2"
config["singularity"]["samtools"] = "docker://fischuu/samtools:1.9-0.2"
config["singularity"]["r-gbs"] = "docker://fischuu/r-gbs:4.2.1-0.4"
config["singularity"]["stringtie"] = "docker://fischuu/stringtie:2.2.1-0.1"
config["singularity"]["subread"] = "docker://fischuu/subread:2.0.1-0.1"

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
print("##### pipeline-config : "+config["pipeline-config"])
print("##### server-config   : "+config["server-config"])
print("#####")
print("##### Singularity configuration")
print("##### --------------------------------")
print("##### star      : "+config["singularity"]["star"])
print("##### gbs       : "+config["singularity"]["gbs"])
print("##### cutadapt  : "+config["singularity"]["cutadapt"])
print("##### minimap2  : "+config["singularity"]["minimap2"])
print("##### r-gbs     : "+config["singularity"]["r-gbs"])
print("##### samtools  : "+config["singularity"]["samtools"])
print("##### subread   : "+config["singularity"]["subread"])
print("##### stringtie : "+config["singularity"]["stringtie"])
print("#####")
print("##### Runtime-configurations")
print("##### --------------------------------")
print("##### genome           : "+ config["genome"])
print("##### existing mock    : "+ config["mockreference"])
print("##### Sample sheet     : "+ config["samplesheet-file"])
print("##### Rawdata folder   : "+ config["rawdata-folder"])
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
        "%s/finalReport.html" % (config["project-folder"]),
      # Reference Genome and mock related
        "%s/Stringtie/merged_STRG.gtf" % (config["project-folder"]),
        expand("%s/QUANTIFICATION/Reference_contigs/{samples}_reference_contigs_fc.txt" % (config["project-folder"]), samples=samples)

rule preparations:
    input:
        "%s/chrName.txt" % (config["genome-star-index"]),
        expand("%s/FASTQ/CONCATENATED/{samples}_R1_001.merged.fastq.gz" % (config["project-folder"]), samples=samples),
        expand("%s/FASTQ/CONCATENATED/{samples}_R2_001.merged.fastq.gz" % (config["project-folder"]), samples=samples),
        config["barcodes-file"]

rule QC:
    input:
        "%s/QC/RAW/multiqc_R1/" % (config["project-folder"]),
        "%s/QC/CONCATENATED/multiqc_R1/" % (config["project-folder"]),
        "%s/QC/TRIMMED/multiqc_R1/" % (config["project-folder"]),
        expand("%s/QC/RAW/{rawsamples}_R1_qualdist.txt" % (config["project-folder"]), rawsamples=rawsamples),
        expand("%s/QC/RAW/{rawsamples}_R2_qualdist.txt" % (config["project-folder"]), rawsamples=rawsamples),
        expand("%s/QC/CONCATENATED/{samples}_R1_qualdist.txt" % (config["project-folder"]), samples=samples),
        expand("%s/QC/CONCATENATED/{samples}_R2_qualdist.txt" % (config["project-folder"]), samples=samples),
        expand("%s/QC/TRIMMED/{samples}_R1_qualdist.txt" % (config["project-folder"]), samples=samples),
        expand("%s/QC/TRIMMED/{samples}_R2_qualdist.txt" % (config["project-folder"]), samples=samples),
        "%s/QC-Report.html" % (config["project-folder"])

rule preprocessing:
    input:
        expand("%s/FASTQ/TRIMMED/{samples}.R1.fq.gz" % (config["project-folder"]), samples=samples),
        expand("%s/FASTQ/TRIMMED/{samples}.R2.fq.gz" % (config["project-folder"]), samples=samples),
        expand("%s/FASTQ/SUBSTITUTED/{samples}.R1.fq.gz" % (config["project-folder"]), samples=samples),
        expand("%s/FASTQ/SUBSTITUTED/{samples}.R2.fq.gz" % (config["project-folder"]), samples=samples),

rule mockreference:
    input:
        "%s/FASTQ/TRIMMED/GSC.MR.Genome.fa" % (config["project-folder"]),
        "%s/FASTQ/TRIMMED/GSC.MR.Clusters.fa" % (config["project-folder"])

rule readalignment:
    input:
        expand("%s/FASTQ/TRIMMED/alignments_reference/{samples}.sorted.bam" % (config["project-folder"]), samples=samples),
        expand("%s/MPILEUP/mpileup_reference/{samples}.mpileup" % (config["project-folder"]), samples=samples),
        expand("%s/FASTQ/TRIMMED/alignments_reference/{samples}.sam.flagstat" % (config["project-folder"]), samples=samples),
        expand("%s/FASTQ/TRIMMED/alignments/{samples}.sorted.bam" % (config["project-folder"]), samples=samples),
        expand("%s/FASTQ/TRIMMED/alignments_clusters/{samples}.sorted.bam" % (config["project-folder"]), samples=samples),
        expand("%s/FASTQ/TRIMMED/{samples}.mpileup" % (config["project-folder"]), samples=samples),
        expand("%s/FASTQ/TRIMMED/alignments/{samples}.sam.flagstat" % (config["project-folder"]), samples=samples),
        expand("%s/FASTQ/TRIMMED/alignments_clusters/{samples}.sam.flagstat" % (config["project-folder"]), samples=samples),
        expand("%s/BAM/alignments_finalMock/{samples}.sorted.bam" % (config["project-folder"]), samples=samples),
        expand("%s/MPILEUP/mpileup_finalMock/{samples}.mpileup" % (config["project-folder"]), samples=samples),
        expand("%s/BAM/alignments_finalMock/{samples}.sam.flagstat" % (config["project-folder"]), samples=samples)

rule callvariants:
    input:
        expand("%s/FASTQ/TRIMMED/{samples}.ref.txt" % (config["project-folder"]), samples=samples),
        expand("%s/MPILEUP/mpileup_reference/{samples}.count.txt" % (config["project-folder"]), samples=samples),
        "%s/FASTQ/TRIMMED/VerticalRefPos.txt" % (config["project-folder"]),
        "%s/MPILEUP/mpileup_reference/VerticalRefPos.txt" % (config["project-folder"]),
        "%s/FASTQ/TRIMMED/GSC.MasterMatrix.txt" % (config["project-folder"]),
        "%s/MPILEUP/mpileup_reference/GSC.MasterMatrix.txt" % (config["project-folder"]),
        "%s/FASTQ/TRIMMED/variants/GSC.GenoMatrix.txt" % (config["project-folder"]),
        "%s/MPILEUP/mpileup_reference/variants/GSC.GenoMatrix.txt" % (config["project-folder"]),
        "%s/FASTQ/TRIMMED/GSC.vcf" % (config["project-folder"]),
        "%s/VCF/FinalSetVariants_referenceGenome.vcf" % (config["project-folder"]),
        expand("%s/MPILEUP/mpileup_finalMock/{samples}.count.txt" % (config["project-folder"]), samples=samples),
        "%s/MPILEUP/mpileup_finalMock/VerticalRefPos.txt" % (config["project-folder"]),
        "%s/MPILEUP/mpileup_finalMock/GSC.MasterMatrix.txt" % (config["project-folder"]),
        "%s/MPILEUP/mpileup_finalMock/variants/GSC.GenoMatrix.txt" % (config["project-folder"]),
        "%s/VCF/FinalSetVariants_finalMock.vcf" % (config["project-folder"])

rule postprocessing:
    input:
        "%s/FASTQ/TRIMMED/GSC.vcf.fa" % (config["project-folder"]),
        "%s/MPILEUP/mpileup_reference/GSC.vcf.fa" % (config["project-folder"]), 
        "%s/SAM/mockVariantsToReference/mockVariantsToReference.sam" % (config["project-folder"]),
        "%s/BAM/mockVariantsToReference/mockVariantsToReference.sorted.bam" % (config["project-folder"])

rule finalreport:
    input:
        "%s/finalReport.html" % (config["project-folder"])

rule MockEval:
    input:
        "%s/mockEvalReport.html" % (config["project-folder"])

rule FinalReport:
    input:
        "%s/finalReport.html" % (config["project-folder"])

rule MockRefVCF:
    input:
        expand("%s/MPILEUP/mpileup_existingMock/{samples}.vcf.gz" % (config["project-folder"]), samples=samples),
        expand("%s/BAM/alignments_existingMock/{samples}.sam.flagstat" % (config["project-folder"]), samples=samples),
    #    "%s/VCF/FinalSetVariants_finalMock.vcf" % (config["project-folder"]),
                    "%s/VCF/variants_existingMock.vcf" % (config["project-folder"])

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
include: "rules/Module9-ReferenceGenome"
