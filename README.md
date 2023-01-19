# Snakebite-GBS
The snakemake pipeline based on GBS-SNP-Crop and further extended.

# DAG
The visual overview of the various rules of the pipeline
![alt text](https://github.com/fischuu/Pipeline-GBS/blob/main/workflow.png?raw=true)

# Configuration
The pipeline parameter configuration takes place in the file

`Snakebite-GBS_config.yaml`

Further, the settings important for the computing cluster are defined here:

`Snakebite-GBS_server-config.yaml`

# Starting the pipeline

Currently, the default configuration is such that the pipelines makes use of docker/singularity. If you plan to use local software installation and or conda environments, you would need to install the required tools manually.

Once the environment is created, you are ready to go. 

To start the pipeline, just type

`bash run_Snakebite-GBS.sh`

You can initiate a dry-run by running

`bash run_Snakebite-GBS.sh -np`

# More information
For further information, please visit the wiki page here:
https://github.com/fischuu/Snakebite-GBS/wiki

# Current state
The pipeline is currently under active development, unstable development versions with the latest features
can be found in the `dev` branch (odd version numbers), the `main` branch should be the latest stable version (even version numbers).

# Possible use-cases

## Complete de-novo approach
No previous information available, just a set of FASTQ-files as input

## Reference genome available
Here, we have the FASTQ files and an existing reference genome in FASTA format

## Mock reference anmd variant set available
Here, we ran a previous de-novo approach and obtained from there a set of variants in VCF format and also a mock reference in fasta format

## Reference genome and variant set available
A reference genome as well as a set of variants is available.

# FUNDING
The development of this pipeline received funding from the ARCTAQUA project (https://site.nord.no/arctaqua/) and GENOTYPE, Natural Resources Institute Finland (LUKE) (https://www.luke.fi) 

