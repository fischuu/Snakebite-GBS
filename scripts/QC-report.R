if(!is.element("snakemake",ls())){
  projFolder <- "/scratch/project_2001746/TestProject"
  pipelineFolder <- "~/git/Pipeline-GBS/"
  pipelineConfig.file <- "/scratch/project_2001746/TestProject/GBS-pipeline_config.yaml"
  refGenome.file <- "hermetiaRef_112020.fasta"
}

createRMD.command <- paste0("cat ",pipelineFolder,"/scripts/QC-header.Rmd ",
                                   pipelineFolder,"/scripts/helpFunctions.Rmd ",
                                   pipelineFolder,"/scripts/generalWorkflow.Rmd ",
                                   pipelineFolder,"/scripts/basicStats.Rmd ",
                                   pipelineFolder,"/scripts/QC.Rmd ",
                                   "> ",projFolder,"/QC-Report.Rmd")

system(createRMD.command)

rmarkdown::render(file.path(projFolder,"QC-Report.Rmd"), output_file=file.path(projFolder,"QC-Report.html"))
