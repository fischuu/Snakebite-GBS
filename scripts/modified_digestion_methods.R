# Modified from of the SimRAD by Olivier Lepais and Jason Weir(2016),
# https://CRAN.R-project.org/package=SimRAD

# original: insilico.digest
insilico.dddigest <- function (DNAseq, cut_site_5prime1, cut_site_3prime1, cut_site_5prime2, 
                               cut_site_3prime2, verbose = TRUE) 
{
  recognition_code1 <- paste(cut_site_5prime1, cut_site_3prime1, 
                             sep = "")
  digest1 <- insilico.dddigest.internal(DNAseq, recognition_code1)
  if (cut_site_5prime2 == "NULL") {
    RESULT <- digest1
  }
  if (cut_site_5prime2 != "NULL") {
    recognition_code2 <- paste(cut_site_5prime2, cut_site_3prime2, 
                               sep = "")
    digest2 <- insilico.dddigest.internal(digest1, recognition_code2)
    RESULT <- digest2
  }
  
  RESULT <- unlist(RESULT)
  if (verbose == TRUE) {
    cat("Number of restriction sites for the first enzyme: ", 
        length(insilico.dddigest.internal(DNAseq, recognition_code1)) - 1, "\n", 
        sep = "")
    cat("Number of restriction sites for the second enzyme: ", 
        length(insilico.dddigest.internal(DNAseq, recognition_code2)) - 1, "\n", 
        sep = "")
    dig1 <- RESULT[isMatchingStartingAt(recognition_code1, 
                                        RESULT)]
    dg1 <- reverseComplement(DNAStringSet(dig1))
    re2match <- reverseComplement(DNAStringSet(recognition_code2))
    dg2 <- dg1[isMatchingStartingAt(re2match[[1]], dg1)]
    dig2 <- reverseComplement(dg2)
    dig1bis <- RESULT[isMatchingStartingAt(recognition_code2, 
                                           RESULT)]
    dg1bis <- reverseComplement(DNAStringSet(dig1bis))
    re2matchbis <- reverseComplement(DNAStringSet(recognition_code1))
    dg2bis <- dg1bis[isMatchingStartingAt(re2matchbis[[1]], 
                                          dg1bis)]
    dig2bis <- reverseComplement(dg2bis)
    cat("Number of type AB and BA fragments:", 
        length(dig2) + length(dig2bis), "\n", sep = "")
    RE1RE1.re2match <- reverseComplement(DNAStringSet(recognition_code1))
    RE1RE1.dg2 <- dg1[isMatchingStartingAt(RE1RE1.re2match[[1]], 
                                           dg1)]
    RE1RE1.dig2 <- reverseComplement(RE1RE1.dg2)
    cat("Number of type AA fragments:", length(RE1RE1.dig2), 
        "\n", sep = "")
    dig3 <- RESULT[isMatchingStartingAt(recognition_code2, 
                                        RESULT)]
    dg3 <- reverseComplement(DNAStringSet(dig3))
    RE2RE2.re2match <- reverseComplement(DNAStringSet(recognition_code2))
    RE2RE2.dg3 <- dg3[isMatchingStartingAt(RE2RE2.re2match[[1]], 
                                           dg3)]
    RE2RE2.dig3 <- reverseComplement(RE2RE2.dg3)
    cat("Number of type BB fragments:", length(RE2RE2.dig3), 
        "\n", sep = "")
    
    
  }
  return(RESULT)
}


# original: insilico.digest.internal
insilico.dddigest.internal<- function (DNAseq, recognition_code) 
{
  frag1 <- strsplit(DNAseq, split = recognition_code, fixed = FALSE, 
                    perl = FALSE)
  n <- length(frag1)
  for (i in 1:n) {
    ni <- length(frag1[[i]])
    if (ni == 1) {
      frag1[[i]] <- frag1[[i]]
    }
    if (ni > 1) {
      for (y in 1:ni) {
        if (y == 1) {
          frag1[[i]][1] <- paste(frag1[[i]][1], recognition_code, 
                                 sep = "")
        }
        if (y > 1 & y < ni) {
          frag1[[i]][y] <- paste(recognition_code, frag1[[i]][y], 
                                 recognition_code, sep = "")
        }
        if (y == ni) {
          frag1[[i]][y] <- paste(recognition_code, frag1[[i]][y], 
                                 sep = "")
        }
      }
    }
  }
  if (n == 1 & ni == 1) {
    frag2 <- frag1
  }
  frag2 <- unlist(frag1)
  return(frag2)
}

# original: adapt.select
myadapt.ABBAselect <- function (sequences, cut_site_5prime1, 
                              cut_site_3prime1, cut_site_5prime2, cut_site_3prime2) 
{
  recognition_code1 <- paste(cut_site_5prime1, cut_site_3prime1, 
                             sep = "")
  recognition_code2 <- paste(cut_site_5prime2, cut_site_3prime2, 
                             sep = "")
  dig1 <- sequences[isMatchingStartingAt(recognition_code1, 
                                         sequences)]
  dg1 <- reverseComplement(DNAStringSet(dig1))
  re2match <- reverseComplement(DNAStringSet(recognition_code2))
  dg2 <- dg1[isMatchingStartingAt(re2match[[1]], dg1)]
  dig2 <- reverseComplement(dg2)
  dig1bis <- sequences[isMatchingStartingAt(recognition_code2, 
                                            sequences)]
  dg1bis <- reverseComplement(DNAStringSet(dig1bis))
  re2matchbis <- reverseComplement(DNAStringSet(recognition_code1))
  dg2bis <- dg1bis[isMatchingStartingAt(re2matchbis[[1]], 
                                        dg1bis)]
  dig2bis <- reverseComplement(dg2bis)
  cat("Number of type AB and BA fragments:", 
      length(dig2) + length(dig2bis), "\n", sep = "")
  
  dig3 <- c(dig2, dig2bis)
  return(dig3)
}
