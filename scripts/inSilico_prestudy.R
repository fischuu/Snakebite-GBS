# Developing version options, delete later
  refGenome.file <- "GCF_901001165.1_fSalTru1.1_genomic.fna"
  refGenome.file <- "Genomes/GCA_901001165.2_fSalTru1.2_genomic.fna"
  projFolder <- "/scratch/project_2005222/BrownInSilico"
  enz1.in <- "CTGCAG" 
  enz2.in <- "GGATCC"
  pipeFolder <- "/users/fischerd/git/Snakebite-GBS"

# Prepare the output matrix
  cut.sites <- c(seq(150,700,50))
  resultOverview <- as.data.frame(t(combn(cut.sites,2)))
  colnames(resultOverview) <- c("minSize", "maxSize")
  resultOverview$totalFragments <- 0
  resultOverview$recommendedReads <- 0


# Updated importFA file, added here for sake of lazyness, it is updated in the
# GenomicTools.fileHandler package also and should be used from there.
  source(file.path(pipeFolder, "scripts", "modified_importFA.R"))

# Import required libraries
  library("SimRAD")
  library("GenomicTools")

# Set additional parameters and import updated functionalities
  options(scipen = 999)
  pipeFolder <- "/users/fischerd/git/Snakebite-GBS"
  source(file.path(pipeFolder, "scripts", "modified_digestion_methods.R"))

# Adjust the paths
  if(substr(refGenome.file,1,1)!="/") refGenome.file <- file.path(projFolder, refGenome.file)

# The SimRAD default import function for fasta sequences fails for larger genomes,
# so I use an own one here
# WARNING!!!! IN PRODUCTION, SWITCH THAT FUNCTION BACK TO THE REGULAR, UPDATED
# ONE AS IT IS IN THE GENOMICTOOLS.FILEHANDLER PACKAGE AVAILABLE
  fasta <- importFA2(refGenome.file, toupper=TRUE, verbose=TRUE)
  names(fasta) <- NULL

# Adjust the enzyme sequence based on separator symbol
  if(grepl("'", enz1.in)){
    enz1 <- strsplit(enz1.in, "'")
  } else {
    enz1 <- strsplit(enz1.in, "")  
  }

  if(grepl("'", enz2.in)){
    enz2 <- strsplit(enz2.in, "'")
  } else {
    enz2 <- strsplit(enz2.in, "")  
  }

  enz1.1 <- enz1[[1]][1]
  enz1.2 <- paste(enz1[[1]][-1], collapse="")
  enz2.1 <- enz2[[1]][1]
  enz2.2 <- paste(enz2[[1]][-1], collapse="")

# Perform the insilico digestion
# Remember, this function is based on an updated version of the original SimRAD
# package, and is imported with the source command earlier in this script
  ddout <- insilico.dddigest(fasta, cut_site_5prime1 = enz1.1, 
                                    cut_site_3prime1 = enz1.2,
                                    cut_site_5prime2 = enz2.1,
                                    cut_site_3prime2 = enz2.2)

  names(ddout) <- paste0("> Location", 1:length(ddout))
  refseq.selA <- myadapt.ABBAselect(ddout, enz1.1, enz1.2, enz2.1, enz2.2) 

# Wrapper function for the prestudy estimation
  get_prestudy_tests <- function(lowerEnd, upperEnd, 
                                 flanking_uncertainty=30,
                                 fragment_factor_uncertainty=2,
                                 verbose=TRUE){
    
      
      
      ddout.selected <- size.select(refseq.selA, min.size = lowerEnd - flanking_uncertainty, max.size = upperEnd + flanking_uncertainty, graph = FALSE)
      
      number_fragments <- length(ddout.selected)
      
      # Parameters
      N <- number_fragments  # number of clusters
      individual_callrate <- 0.75
      target_global_callrate <- 0.75
      samples <- 150
      r <- 10     # minimum number of realizations per cluster
      # Fi
      p <- 0.8   # desired proportion of clusters with at least r realizations
      reps <- 100
      
      for(p_test in c(0.5, 0.6, 0.7, 0.8, 0.9)){
        call_rate <- c()
        for(i in 1:reps){
          variant_calls <- apply(matrix(rbinom(N*samples,1,p_test), ncol=N),2,sum)/samples > individual_callrate
          call_rate[i] <- sum(variant_calls) / N
        }
        
        global_callrate <- mean(call_rate)
        if(global_callrate >= target_global_callrate){
          p <- p_test
          message("Set desired cluster coverage proportion with at least r realisations to ", p)
          break
        }
      }
      
      # Now determine the suggested number of reads per sample
      
      # Initialize lambda
      lambda <- r
      
      # Function to simulate one round of sampling
      simulate <- function(N, lambda) {
        # Simulate the number of draws per cluster
        draws <- rpois(N, lambda)
        
        # Return the proportion of clusters with at least r realizations
        mean(draws >= r)
      }
      
      
      # Loop until the desired coverage is achieved
      while (TRUE) {
        # Simulate a round of sampling
        prop <- mean(replicate(100, simulate(N, lambda)))
        
        # Check if the desired coverage is achieved
        if (prop >= p) {
          break
        }
        
        # Increase lambda
        lambda <- lambda + 1
      }
      
      # Print the estimated number of total draws
      total_draws <- fragment_factor_uncertainty * N * lambda
      if(verbose) cat("Based on the input, the suggested minimum number of reads per sample should be", prettyNum(total_draws, big.mark = ","), "\n")
      
      list(totalFragments = number_fragments,
           recommendedReads = total_draws)
      
  }

for(i in 1:nrow(resultOverview)){
  tmp <- get_prestudy_tests(lowerEnd = resultOverview$minSize[i], resultOverview$maxSize[i], verbose = FALSE)
  resultOverview$totalFragments[i] <- tmp$totalFragments
  resultOverview$recommendedReads[i] <- tmp$recommendedReads
}
