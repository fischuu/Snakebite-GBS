if(!is.element("snakemake",ls())){
  projFolder <- "/scratch/project_2001746/TestProject"
  pipelineFolder <- "~/git/Pipeline-GBS/"
  pipelineConfig.file <- "/scratch/project_2001746/TestProject/GBS-pipeline_config.yaml"
  refGenome.file <- "hermetiaRef_112020.fasta"
}

last_n_char <- function(x, n=1){
  substr(x, nchar(x) - n + 1, nchar(x))
}

if(last_n_char(pipelineFolder)!="/") pipelineFolder <- paste0(pipelineFolder, "/")


if(refGenome.file == ""){
  createRMD.command <- paste0("cat ",pipelineFolder,"scripts/Rmodules/final-header.Rmd ",
                                     pipelineFolder,"scripts/Rmodules/helpFunctions.Rmd ",
                                     pipelineFolder,"scripts/Rmodules/generalWorkflow.Rmd ",
                                     pipelineFolder,"scripts/Rmodules/basicStats.Rmd ",
                                     pipelineFolder,"scripts/Rmodules/QC.Rmd ",
                                     pipelineFolder,"scripts/Rmodules/clusterAndReference_mockPart.Rmd ",
                                     pipelineFolder,"scripts/Rmodules/alignments_mockPart.Rmd ",
                                     pipelineFolder,"scripts/Rmodules/variants_mockPart.Rmd ",
                                     pipelineFolder,"scripts/Rmodules/homology.Rmd ",
                                     pipelineFolder,"scripts/Rmodules/caseControl.Rmd ",
                                     pipelineFolder,"scripts/Rmodules/warnings.Rmd ",
                                     pipelineFolder,"scripts/Rmodules/benchmarks.Rmd ",
                                     pipelineFolder,"scripts/Rmodules/appendix.Rmd ",
                                     pipelineFolder,"scripts/Rmodules/references.Rmd ",
                                     "> ",projFolder,"/finalReport.Rmd")
} else {
  createRMD.command <- paste0("cat ",pipelineFolder,"scripts/Rmodules/final-header.Rmd ",
                                     pipelineFolder,"scripts/Rmodules/helpFunctions.Rmd ",
                                     pipelineFolder,"scripts/Rmodules/generalWorkflow.Rmd ",
                                     pipelineFolder,"scripts/Rmodules/basicStats.Rmd ",
                                     pipelineFolder,"scripts/Rmodules/QC.Rmd ",
                                     pipelineFolder,"scripts/Rmodules/clusterAndReference_mockPart.Rmd ",
                                     pipelineFolder,"scripts/Rmodules/clusterAndReference_refPart.Rmd ",
                                     pipelineFolder,"scripts/Rmodules/alignments_mockPart.Rmd ",
                                     pipelineFolder,"scripts/Rmodules/alignments_refPart.Rmd ",
                                     pipelineFolder,"scripts/Rmodules/Insilico.Rmd ",
                                     pipelineFolder,"scripts/Rmodules/variants_mockPart.Rmd ",
                                     pipelineFolder,"scripts/Rmodules/variants_refPart.Rmd ",
                                     pipelineFolder,"scripts/Rmodules/homology.Rmd ",
                                     pipelineFolder,"scripts/Rmodules/caseControl.Rmd ",
                                     pipelineFolder,"scripts/Rmodules/warnings.Rmd ",
                                     pipelineFolder,"scripts/Rmodules/benchmarks.Rmd ",
                                     pipelineFolder,"scripts/Rmodules/appendix.Rmd ",
                                     pipelineFolder,"scripts/Rmodules/references.Rmd ",
                                     "> ",projFolder,"/finalReport.Rmd")
}

cat(createRMD.command, "\n")

system(createRMD.command)

rmarkdown::render(file.path(projFolder,"finalReport.Rmd"), output_file=file.path(projFolder,"finalReport.html"))