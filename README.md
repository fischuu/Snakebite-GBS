# Pipeline-GBS
The snakemake pipeline based on GBS-SNP-Crop and further extended.

# Configuration
The pipeline parameter configuration takes place in the file

`GBS-pipeline_config.yaml` (there is an example here: `GBS-pipeline_config-example.yaml`)

Further, the settings important for the computing cluster are defined here:

`GBS-pipeline_server-config.yaml` (example: `GBS-pipeline_puhti-config.yaml`)

# Starting the pipeline
Before you start the pipeline, make sure you create the corresponding conda environment,
as defined in this file

`GBS.yaml`

Once the environment is created, you are ready to go. 

To start the pipeline, just type

`bash runGBSPipeline-example.sh`

You can initiate a dry-run by running

`bash runGBSPipeline-example.sh -np`

# More information
For further information, please visit the wiki page here:
https://github.com/fischuu/Pipeline-GBS/wiki

# Current state
The pipeline is currently under active development, unstable development versions with the latest features
can be found in the `dev` branch (odd version numbers), the `main` branch should be stable (even version numbers).

