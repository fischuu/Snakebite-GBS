if(!is.element("snakemake",ls())){
  projFolder <- "/scratch/project_2001746/TestProject"
  pipelineFolder <- "~/git/Pipeline-GBS/"
  pipelineConfig.file <- "/scratch/project_2001746/TestProject/GBS-pipeline_config.yaml"
  refGenome.file <- "hermetiaRef_112020.fasta"
}

createRMD.command <- paste0("cat ",pipelineFolder,"scripts/Rmodules/final-header.Rmd ",
                                   pipelineFolder,"scripts/Rmodules/helpFunctions.Rmd ",
                                   pipelineFolder,"scripts/Rmodules/generalWorkflow.Rmd ",
                                   pipelineFolder,"scripts/Rmodules/basicStats.Rmd ",
                                   pipelineFolder,"scripts/Rmodules/QC.Rmd ",
                                   pipelineFolder,"scripts/Rmodules/clusterAndReference.Rmd ",
                                   pipelineFolder,"scripts/Rmodules/alignments.Rmd ",
                                   pipelineFolder,"scripts/Rmodules/variants.Rmd ",
                                   pipelineFolder,"scripts/Rmodules/homology.Rmd ",
                                   pipelineFolder,"scripts/Rmodules/caseControl.Rmd ",
                                   pipelineFolder,"scripts/Rmodules/warnings.Rmd ",
                                   pipelineFolder,"scripts/Rmodules/benchmarks.Rmd ",
                                   pipelineFolder,"scripts/Rmodules/appendix.Rmd ",
                                   pipelineFolder,"scripts/Rmodules/references.Rmd ",
                                   "> ",pipelineFolder,"scripts/final-report.Rmd")

system(createRMD.command)

rmarkdown::render(file.path(pipelineFolder,"scripts","final-report.Rmd"), output_file=file.path(pipelineFolder,"finalReport.html"))