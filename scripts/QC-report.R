if(!is.element("snakemake",ls())){
  projFolder <- "/scratch/project_2001746/TestProject"
  pipelineFolder <- "~/git/Pipeline-GBS/"
  pipelineConfig.file <- "/scratch/project_2001746/TestProject/GBS-pipeline_config.yaml"
  refGenome.file <- "hermetiaRef_112020.fasta"
}

createRMD.command <- paste0("cat ",pipelineFolder,"/scripts/Rmodules/QC-header.Rmd ",
                                   pipelineFolder,"/scripts/Rmodules/helpFunctions.Rmd ",
                                   pipelineFolder,"/scripts/Rmodules/generalWorkflow.Rmd ",
                                   pipelineFolder,"/scripts/Rmodules/basicStats.Rmd ",
                                   pipelineFolder,"/scripts/Rmodules/QC.Rmd ",
                                   "> ",projFolder,"/QC-Report.Rmd")

system(createRMD.command)

rmarkdown::render(file.path(projFolder,"QC-Report.Rmd"), output_file=file.path(projFolder,"QC-Report.html"))
