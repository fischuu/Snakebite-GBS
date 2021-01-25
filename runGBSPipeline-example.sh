# run the Genotype snakemake pipeline on the taito cluster
#
# Before you start, take care that the following two files are available:
# rawsamples, barcodesID.txt, sampleInfo.txt
#
# Example shell scripts on how to create them, you find below

# rawsamples
###########################
# find /scratch/project_2001746/Example/FASTQ/RAW/ -name '*R1_001.fastq.gz' | xargs -n1 basename | sed 's/_R1_001.fastq.gz//g' >  /scratch/project_2001746/Example/rawsamples
# FORMAT: One raw sample file per row

# barcodesID.txt
###########################
# grep -i "Identifier" /scratch/project_2001881/191009_NB551722_0004_AHT2JYAFXY/191009_NB551722_0004_AHT2JYAFXY/SampleSheet.csv | \
# awk 'BEGIN {FS=","; OFS="\t"} {print $8"+"$6, $1 ,"NO"}' > /scratch/project_2001746/Example/barcodesID.txt
# FORMAT:   INDEX \tab SampleName \tab BOOL, include in mock reference?
# $1 is target file

# sampleInfo.txt
##############################
# cut -f2 barcodesID.txt > samples
# cut -f1 samples | cut -d'_' -f1 | sed 's/kaapio/case/g' | sed 's/sisar/control/g' | sed 's/isa/control/g' | sed 's/ema/control/g' > caseControl
# cut -f1 samples | cut -d'_' -f2 | sed 's/[^0-9]*//g' > famID
# cut -f1 samples | cut -d'_' -f1 > role
# paste -d '\t' samples caseControl famID role > sampleInfo.txt
# sed -i '1i sample\tcasecontrol\tfamily\trole' sampleInfo.txt
# rm caseControl famID role samples
#
# 
#

module load bioconda/3
# Uncomment this, once the environemnt is created
source activate /projappl/project_2001746/conda_envs/Genotype
#source activate /projappl/project_2001289/FAANGlncRNA

# Create the rulegraph
#snakemake -s GBS-pipeline.smk \
#          --configfile /scratch/project_2001746/Pipeline-GBS/GBS-pipeline_config-example.yaml \
#          --rulegraph | dot -T png > ./workflow.png

snakemake -s /scratch/project_2001746/Pipeline-GBS/GBS-pipeline.smk \
          -j 150 \
          --use-conda \
          --use-singularity \
          --configfile /scratch/project_2001746/Pipeline-GBS/GBS-pipeline_config-example.yaml \
          --latency-wait 60 \
          --cluster-config /scratch/project_2001746/Pipeline-GBS/GBS-pipeline_puhti-config.yaml \
          --cluster "sbatch -t {cluster.time} --account={cluster.account} --gres=nvme:{cluster.nvme} --job-name={cluster.job-name} --tasks-per-node={cluster.ntasks} --cpus-per-task={cluster.cpus-per-task} --mem-per-cpu={cluster.mem-per-cpu} -p {cluster.partition} -D {cluster.working-directory}" \
          $1 
