if(!is.element("snakemake",ls())){
  projFolder <- ""
  pipelineFolder <- ""
  pipelineConfig.file <- ""
  refGenome.file <- ""
}

if(refGenome.file == ""){
  createRMD.command <- paste0("cat ",pipelineFolder,"/scripts/Rmodules/VariantCalling-header.Rmd ",
                              pipelineFolder,"/scripts/Rmodules/helpFunctions.Rmd ",
                              pipelineFolder,"/scripts/Rmodules/generalWorkflow.Rmd ",
                              pipelineFolder,"/scripts/Rmodules/basicStats.Rmd ",
                              pipelineFolder,"/scripts/Rmodules/variants_mockPart.Rmd ",
                              "> ",projFolder,"/VariantCalling-Report.Rmd")
} else {
  createRMD.command <- paste0("cat ",pipelineFolder,"/scripts/Rmodules/VariantCalling-header.Rmd ",
                              pipelineFolder,"/scripts/Rmodules/helpFunctions.Rmd ",
                              pipelineFolder,"/scripts/Rmodules/generalWorkflow.Rmd ",
                              pipelineFolder,"/scripts/Rmodules/basicStats.Rmd ",
                              pipelineFolder,"/scripts/Rmodules/variants_mockPart.Rmd ",
                              pipelineFolder,"/scripts/Rmodules/variants_refPart.Rmd ",
                              "> ",projFolder,"/VariantCalling-Report.Rmd")
}


system(createRMD.command)

rmarkdown::render(file.path(projFolder,"VariantCalling-Report.Rmd"), output_file=file.path(projFolder,"VariantCalling-Report.html"))
