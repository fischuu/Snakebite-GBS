if(!is.element("snakemake",ls())){
  projFolder <- "/scratch/project_2001746/TestProject"
  pipelineFolder <- "~/git/Pipeline-GBS/"
  pipelineConfig.file <- "/scratch/project_2001746/TestProject/GBS-pipeline_config.yaml"
  refGenome.file <- "hermetiaRef_112020.fasta"
}

createRMD.command <- paste0("cat ",pipelineFolder,"scripts/final-header.Rmd ",
                                   pipelineFolder,"scripts/helpFunctions.Rmd ",
                                   pipelineFolder,"scripts/generalWorkflow.Rmd ",
                                   pipelineFolder,"scripts/basicStats.Rmd ",
                                   pipelineFolder,"scripts/QC.Rmd ",
                                   pipelineFolder,"scripts/clusterAndReference.Rmd ",
                                   pipelineFolder,"scripts/alignments.Rmd ",
                                   pipelineFolder,"scripts/variants.Rmd ",
                                   pipelineFolder,"scripts/homology.Rmd ",
                                   pipelineFolder,"scripts/caseControl.Rmd ",
                                   pipelineFolder,"scripts/warnings.Rmd ",
                                   pipelineFolder,"scripts/benchmarks.Rmd ",
                                   pipelineFolder,"scripts/appendix.Rmd ",
                                   pipelineFolder,"scripts/references.Rmd ",
                                   "> ",pipelineFolder,"scripts/final-report.Rmd",)

system(createRMD.command)

rmarkdown::render(file.path(pipelineFolder,"scripts","final-report.Rmd"), output_file=file.path(pipelineFolder,"finalReport.html"))"