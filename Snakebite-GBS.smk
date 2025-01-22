import pandas as pd
from snakemake.utils import validate, min_version
from multiprocessing import cpu_count
import glob
import re
import os
import sys
import yaml

##### GBS-snakemake pipeline #####
##### Daniel Fischer (daniel.fischer@luke.fi)
##### Natural Resources Institute Finland (Luke)
##### This pipeline is build upon the the GBS-SNP-CROP pipeline:
##### https://github.com/halelab/GBS-SNP-CROP
##### Version: 0.21.2
version = "0.21.2"

##### set minimum snakemake version #####
min_version("6.0")

##### Fill the configuration lines for relative paths
# project-folder should not end with "/", so remove it

if config["project-folder"][-1] == '/':
   config["project-folder"]=config["project-folder"][:-1]
   
if(config["pipeline-config"][0]!='/'):
    config["pipeline-config"] = config["project-folder"] + '/' + config["pipeline-config"]

if(config["server-config"][0]!='/'):
    config["server-config"] = config["project-folder"] + '/' + config["server-config"]

if(config["rawdata-folder"][0]!='/'):
    config["rawdata-folder"] = config["project-folder"] + '/' + config["rawdata-folder"]

if(config["samplesheet-file"][0]!='/'):
    config["samplesheet-file"] = config["project-folder"] + '/' + config["samplesheet-file"]

if config["sampleinfo-file"] == "":
    pass
else:
    if(config["sampleinfo-file"][0]!='/'):
        config["sampleinfo-file"] = config["project-folder"] + '/' + config["sampleinfo-file"]
    
if config["genome"] == "":
    pass
else:
    if(config["genome"][0]!='/'):
        config["genome"] = config["project-folder"] + '/' + config["genome"]

if(config["tmpdir"][0]!='/'):
    config["tmpdir"] = config["project-folder"] + '/' + config["tmpdir"]

##### load config and sample sheets #####

samplesheet = pd.read_table(config["samplesheet-file"]).set_index("rawsample", drop=False)
rawsamples=list(samplesheet.rawsample)
samples=list(set(list(samplesheet.sample_name)))
lane=list(samplesheet.lane)

# Get the basename fastq inputs
possible_ext = [".fastq", ".fq.gz", ".fastq.gz", ".fasta", ".fa", ".fa.gz", ".fasta.gz"]
ext = ".null"

reads1_tmp = list(samplesheet.read1)
reads1_trim = []
for r in reads1_tmp:
    for e in possible_ext:
        if r.endswith(e):
            addThis = r[:-len(e)]
            reads1_trim += [addThis]
            ext=e

reads2_tmp = list(samplesheet.read2)
reads2_trim = []
if(config["libtype"][0]=='PE'):
    for r in reads2_tmp:
        for e in possible_ext:
            if r.endswith(e):
                addThis = r[:-len(e)]
                reads2_trim += [addThis] 
                ext=e

mockSamples = list(samplesheet.useForMock)
samplesUsedForMock = str(mockSamples.count('YES'))


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

##### Extract the cluster resource requests from the server config #####
cluster=dict()
if os.path.exists(config["server-config"]):
    with open(config["server-config"]) as yml:
        cluster = yaml.load(yml, Loader=yaml.FullLoader)

workdir: config["project-folder"]

wildcard_constraints:
    rawsamples="|".join(rawsamples),
    samples="|".join(samples),
    reads1_trim="|".join(reads1_trim),
    reads2_trim="|".join(reads2_trim)

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

def get_raw_input_read_bs(wildcards):
    reads = wildcards.reads + ext
    output = config["rawdata-folder"] + "/" + reads
    return output

def get_raw_input_read_bs1(wildcards):
    reads = wildcards.reads1 + ext
    output = config["rawdata-folder"] + "/" + reads
    return output

def get_raw_input_read_bs2(wildcards):
    reads = wildcards.reads2 + ext
    output = config["rawdata-folder"] + "/" + reads
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

def get_preparations_files(wildcards):
    if config["genome"] == "":
        return []
    else:
        file1 = config["genome-star-index"] + "/chrName.txt"
        file2 = config["genome-bwa-index"]
        
        output = [file1, file2]
        return output   

##### Complete the input configuration
config["genome-bwa-index"] = config["genome"]+".bwt"
config["mockref-bwa-index"] = config["mockreference"]+".bwt"
config["full-insilico-genome"] = config["project-folder"]+"/References/full_inSilico_reference.fa"
config["selected-insilico-genome"] = config["project-folder"]+"/References/sizeSelected_inSilico_reference.fa"
config["genome-bwa-full_insilico-index"]=config["full-insilico-genome"]+".bwt"
config["genome-bwa-selected_insilico-index"]=config["selected-insilico-genome"]+".bwt"
config["genome-star-index"] = config["project-folder"]+"/References/STAR2.7.5a"    # Change here to the path from the reference genome!!!!!!!!!!!!!!!!!!
config["barcodes-script"] = config["pipeline-folder"]+"/scripts/prepareBarcodes.R"
config["report-script"] = config["pipeline-folder"]+"/scripts/Final-report.R"
config["qc-script"] = config["pipeline-folder"]+"/scripts/QC-report.R"
config["variant-script"] = config["pipeline-folder"]+"/scripts/VariantCalling-report.R"
config["mockeval-script"] = config["pipeline-folder"]+"/scripts/mockeval-report.Rmd"
config["refinement-script"] = config["pipeline-folder"]+"/scripts/refineMockReference.R"
config["insilico-script"] = config["pipeline-folder"]+"/scripts/inSilicoFasta.R"
config["similarity-script"] = config["pipeline-folder"]+"/scripts/mockRef_similarity.R"
config["insilico-report-script"] = config["pipeline-folder"]+"/scripts/Insilico-report.R"
config["adapter"]=config["pipeline-folder"]+"/adapter.fa"
config["barcodes-file"] = config["project-folder"]+"/barcodesID.txt"
config["enz1_clean"] = re.sub(r"[^A-Za-z]", "", config["enz1"])
config["enz2_clean"] = re.sub(r"[^A-Za-z]", "", config["enz2"])

##### Singularity container #####
config["singularity"] = {}
config["singularity"]["bedtools"] = "docker://fischuu/bedtools:2.30-0.1"
config["singularity"]["star"] = "docker://fischuu/star:2.7.5a"
config["singularity"]["gbs"] = "docker://fischuu/gbs:0.2"
config["singularity"]["cutadapt"] = "docker://fischuu/cutadapt:2.8-0.3"
config["singularity"]["minimap2"] = "docker://fischuu/minimap2:2.26-0.1"
config["singularity"]["samtools"] = "docker://fischuu/samtools:1.9-0.2"
config["singularity"]["r-gbs"] = "docker://fischuu/r-gbs:4.2.1-0.7"
config["singularity"]["stringtie"] = "docker://fischuu/stringtie:2.2.1-0.1"
config["singularity"]["subread"] = "docker://fischuu/subread:2.0.1-0.1"
config["singularity"]["transanno"] = "docker://informationsea/transanno:0.2.4"

##### Print the welcome screen #####
print("#################################################################################")
print("##### Welcome to the GBS pipeline")
print("##### version: "+version)
print("#####")
print("##### Pipeline configuration")
print("##### --------------------------------")
print("##### project-folder       : "+config["project-folder"])
print("##### pipeline-folder      : "+config["pipeline-folder"])
print("##### report-script        : "+config["report-script"])
print("##### pipeline-config      : "+config["pipeline-config"])
print("##### server-config        : "+config["server-config"])
print("##### Size selection (min) : "+str(config["minLength"]))
print("##### Size selection (max) : "+str(config["maxLength"]))
print("##### Enzyme 1 recog. site : "+config["enz1"])
print("##### Enzyme 2 recog. site : "+config["enz2"])
print("##### Library type         : "+config["libtype"])
print("#####")
print("##### Singularity configuration")
print("##### --------------------------------")
print("##### bdtools   : "+config["singularity"]["bedtools"])
print("##### star      : "+config["singularity"]["star"])
print("##### gbs       : "+config["singularity"]["gbs"])
print("##### cutadapt  : "+config["singularity"]["cutadapt"])
print("##### minimap2  : "+config["singularity"]["minimap2"])
print("##### r-gbs     : "+config["singularity"]["r-gbs"])
print("##### samtools  : "+config["singularity"]["samtools"])
print("##### subread   : "+config["singularity"]["subread"])
print("##### stringtie : "+config["singularity"]["stringtie"])
print("##### transanno : "+config["singularity"]["transanno"])
print("#####")
print("##### Runtime-configurations")
print("##### --------------------------------")
print("##### genome                : "+ config["genome"])
print("##### existing mock         : "+ config["mockreference"])
print("##### Sample sheet          : "+ config["samplesheet-file"])
print("##### Rawdata folder        : "+ config["rawdata-folder"])
print("##### Samples to build mock : " + samplesUsedForMock)
print("#####")
print("##### Derived runtime parameters")
print("##### --------------------------------")
print("##### BWA-Genome index    : "+config["genome-bwa-index"])
print("##### STAR-Genome index   : "+config["genome-star-index"])
print("##### Existing Mock index : "+config["mockref-bwa-index"])
print("##### Adapter file        : "+config["adapter"])
print("#####")
print("##### Output files")
print("##### --------------------------------")
print("##### Size selected in-silico predictions : "+ config["selected-insilico-genome"])
print("##### All in-silico predictions           : "+ config["full-insilico-genome"])
print("#################################################################################")

##### Define conditional input/outputs #####
conditionalOut = list()
if config["mockreference"] != "":
        conditionalOut.append("%s/VCF/FinalSetVariants_existingMock.vcf" % (config["project-folder"]))

    
##### run-time rules for complete pipeline and submodules #####

rule all:
    input:
        conditionalOut,
      # QC OF RAW AND CONCATENATED FILES
        "%s/QC/RAW/multiqc_R1/" % (config["project-folder"]),
        expand("%s/QC/RAW/{rawsamples}_R1_qualdist.txt" % (config["project-folder"]), rawsamples=rawsamples),
        "%s/QC/CONCATENATED/multiqc_R1/" % (config["project-folder"]),
        "%s/QC/TRIMMED/multiqc_R1/" % (config["project-folder"]),
      # OUTPUT STEP 3  
        "%s/REPORTS/DATA/MockReference_Reference_similarity.report.txt" % (config["project-folder"]),
      # OUTPUT STEP 4
        "%s/FASTQ/TRIMMED/GSC.MR.Genome.fa" % (config["project-folder"]),
      # OUTPUT STEP 5
        expand("%s/FASTQ/TRIMMED/alignments/{samples}.sam.flagstat" % (config["project-folder"]), samples=samples),
        expand("%s/FASTQ/TRIMMED/alignments/{samples}.sorted.bam" % (config["project-folder"]), samples=samples),
        expand("%s/FASTQ/TRIMMED/alignments_clusters/{samples}.sam.flagstat" % (config["project-folder"]), samples=samples),
        expand("%s/FASTQ/TRIMMED/alignments_clusters/{samples}.sorted.bam" % (config["project-folder"]), samples=samples),
        expand("%s/FASTQ/TRIMMED/alignments_clusters/{samples}.coverage" % (config["project-folder"]), samples=samples),
      # OUTPUT STEP 6
        "%s/FASTQ/TRIMMED/GSC.MasterMatrix.txt" % (config["project-folder"]),
        expand("%s/FASTQ/TRIMMED/{samples}.count.txt" % (config["project-folder"]), samples=samples),
      # OUTPUT STEP 7
        "%s/FASTQ/TRIMMED/variants/GSC.GenoMatrix.txt" % (config["project-folder"]),
      # OUTPUT STEP 8  
        "%s/FASTQ/TRIMMED/GSC.vcf" % (config["project-folder"]),
        "%s/FASTQ/TRIMMED/GSC.vcf.fa" % (config["project-folder"]),
      # OUTPUT STEP 9
      #  "%s/VCF/FinalSetVariants_finalMock_liftOver-to-Reference_succeeded.vcf" % (config["project-folder"]),
      # Quality check
        expand("%s/BAM/alignments_finalMock/{samples}.sam.flagstat" % (config["project-folder"]), samples=samples),
        "%s/MockReference/MockReference.fa" % (config["project-folder"]),
        "%s/VCF/FinalSetVariants_finalMock.vcf" % (config["project-folder"]),
        "%s/finalReport.html" % (config["project-folder"]),

rule variantsReference:
    input:
        "%s/VCF/FinalSetVariants_referenceGenome.vcf" % (config["project-folder"])

rule liftOver:
    input:
        "%s/MockReference/MockReference.paf" % (config["project-folder"]),
        "%s/VCF/FinalSetVariants_finalMock_liftOver-to-Reference_succeeded.vcf" % (config["project-folder"])

rule insilico:
    input:
        "%s/References/full_inSilico_reference.fa" % (config["project-folder"]),
        "%s/References/sizeSelected_inSilico_reference.fa" % (config["project-folder"]),
        expand("%s/BAM/Insilico/full/{samples}.bam" % (config["project-folder"]), samples=samples),
        expand("%s/BAM/Insilico/selected/{samples}.bam" % (config["project-folder"]), samples=samples),
        expand("%s/BAM/Insilico/full/{samples}.sam.flagstat" % (config["project-folder"]), samples=samples),
        expand("%s/BAM/Insilico/full/{samples}.coverage" % (config["project-folder"]), samples=samples),
        expand("%s/BAM/Insilico/selected/{samples}.sam.flagstat" % (config["project-folder"]), samples=samples),
        expand("%s/BAM/Insilico/selected/{samples}.coverage" % (config["project-folder"]), samples=samples),
        "%s/Insilico-Report.html" % (config["project-folder"])

rule preparations:
    input:
        get_preparations_files,
        expand("%s/FASTQ/CONCATENATED/{samples}_R1_001.merged.fastq.gz" % (config["project-folder"]), samples=samples),
        config["barcodes-file"]

rule datapublication:
    input:
        expand("%s/FASTQ/CONCATENATED/{samples}_R1_001.merged.fastq.gz" % (config["project-folder"]), samples=samples),
        config["barcodes-file"]


rule QC:
    input:
        "%s/QC/RAW/multiqc_R1/" % (config["project-folder"]),
        "%s/QC/CONCATENATED/multiqc_R1/" % (config["project-folder"]),
        "%s/QC/TRIMMED/multiqc_R1/" % (config["project-folder"]),
        expand("%s/QC/RAW/{rawsamples}_R1_qualdist.txt" % (config["project-folder"]), rawsamples=rawsamples),
        expand("%s/QC/CONCATENATED/{samples}_R1_qualdist.txt" % (config["project-folder"]), samples=samples),
        expand("%s/QC/TRIMMED/{samples}_R1_qualdist.txt" % (config["project-folder"]), samples=samples),
        "%s/QC-Report.html" % (config["project-folder"])

rule preprocessing:
    input:
        expand("%s/FASTQ/TRIMMED/{samples}.R1.fq.gz" % (config["project-folder"]), samples=samples),
        expand("%s/FASTQ/SUBSTITUTED/{samples}.R1.fq.gz" % (config["project-folder"]), samples=samples),

rule mockreference:
    input:
        "%s/FASTQ/TRIMMED/GSC.MR.Genome.fa" % (config["project-folder"]),
        "%s/FASTQ/TRIMMED/GSC.MR.Clusters.fa" % (config["project-folder"]),
        "%s/REPORTS/DATA/MockReference_Reference_similarity.report.txt" % (config["project-folder"])
        
def get_readalignment_expand_files(wildcards):
    if config["genome"] == "":
        return []
    else:    
        samples = set(samplesheet["sample_name"])
        path = config["project-folder"]
        output1 = [path + "/FASTQ/TRIMMED/alignments_reference/" + x + ".sorted.bam" for x in samples]
        output2 = [path + "/MPILEUP/mpileup_reference/" + x + ".mpileup" for x in samples]
        output3 = [path + "/FASTQ/TRIMMED/alignments_reference/" + x + ".sam.flagstat" for x in samples]

        output = output1 + output2 + output3
        return output

rule readalignment:
    input:
        get_readalignment_expand_files,  
        expand("%s/FASTQ/TRIMMED/alignments/{samples}.sorted.bam" % (config["project-folder"]), samples=samples),
        expand("%s/FASTQ/TRIMMED/alignments_clusters/{samples}.sorted.bam" % (config["project-folder"]), samples=samples),
        expand("%s/FASTQ/TRIMMED/{samples}.mpileup" % (config["project-folder"]), samples=samples),
        expand("%s/FASTQ/TRIMMED/alignments/{samples}.sam.flagstat" % (config["project-folder"]), samples=samples),
        expand("%s/FASTQ/TRIMMED/alignments_clusters/{samples}.sam.flagstat" % (config["project-folder"]), samples=samples),
        expand("%s/BAM/alignments_finalMock/{samples}.sorted.bam" % (config["project-folder"]), samples=samples),
        expand("%s/MPILEUP/mpileup_finalMock/{samples}.mpileup" % (config["project-folder"]), samples=samples),
        expand("%s/BAM/alignments_finalMock/{samples}.sam.flagstat" % (config["project-folder"]), samples=samples)

def get_callvariants_files(wildcards):
    if config["genome"] == "":
        return []
    else:
        path = config["project-folder"]
        file1 = path + "/MPILEUP/mpileup_reference/VerticalRefPos.txt"
        file2 = path + "/MPILEUP/mpileup_reference/GSC.MasterMatrix.txt"
        file3 = path + "/VCF/FinalSetVariants_referenceGenome.vcf"
        file4 = path + "/MPILEUP/mpileup_reference/variants/GSC.GenoMatrix.txt"
        output = [file1, file2, file3, file4]
        return output   


def get_callvariants_expand_files(wildcards):
    if config["genome"] == "":
        return []
    else:    
        samples = set(samplesheet["sample_name"])
        path = config["project-folder"]
        output1 = [path + "/MPILEUP/mpileup_reference/" + x + ".count.txt" for x in samples]

        output = output1
        return output

rule callvariants:
    input:
        get_callvariants_files,
        get_callvariants_expand_files,
        expand("%s/FASTQ/TRIMMED/{samples}.ref.txt" % (config["project-folder"]), samples=samples),
        "%s/FASTQ/TRIMMED/VerticalRefPos.txt" % (config["project-folder"]),
        "%s/FASTQ/TRIMMED/GSC.MasterMatrix.txt" % (config["project-folder"]),
        "%s/FASTQ/TRIMMED/variants/GSC.GenoMatrix.txt" % (config["project-folder"]),
        "%s/FASTQ/TRIMMED/GSC.vcf" % (config["project-folder"]),
        expand("%s/MPILEUP/mpileup_finalMock/{samples}.count.txt" % (config["project-folder"]), samples=samples),
        "%s/MPILEUP/mpileup_finalMock/VerticalRefPos.txt" % (config["project-folder"]),
        "%s/MPILEUP/mpileup_finalMock/GSC.MasterMatrix.txt" % (config["project-folder"]),
        "%s/MPILEUP/mpileup_finalMock/variants/GSC.GenoMatrix.txt" % (config["project-folder"]),
        "%s/VCF/FinalSetVariants_finalMock.vcf" % (config["project-folder"]),
        "%s/VariantCalling-Report.html" % (config["project-folder"])

rule postprocessing:
    input:
        "%s/FASTQ/TRIMMED/GSC.vcf.fa" % (config["project-folder"]),
        "%s/MPILEUP/mpileup_reference/GSC.vcf.fa" % (config["project-folder"]), 
        "%s/SAM/mockVariantsToReference/mockVariantsToReference.sam" % (config["project-folder"]),
        "%s/BAM/mockVariantsToReference/mockVariantsToReference.sorted.bam" % (config["project-folder"])

rule MockEval:
    input:
        "%s/mockEvalReport.html" % (config["project-folder"])

rule finalReport:
    input:
        "%s/finalReport.html" % (config["project-folder"])

rule MockRefVCF:
    input:
        expand("%s/MPILEUP/mpileup_existingMock/{samples}.vcf.gz" % (config["project-folder"]), samples=samples),
        expand("%s/BAM/alignments_existingMock/{samples}.sam.flagstat" % (config["project-folder"]), samples=samples),
    #   "%s/VCF/FinalSetVariants_finalMock.vcf" % (config["project-folder"]),
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
include: "rules/Module10-InSilico"
