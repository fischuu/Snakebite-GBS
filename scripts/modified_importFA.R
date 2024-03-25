importFA2 <- function(file, toupper=TRUE, verbose=TRUE){
  # Read in the fasta file line by line
    res <- readLines(file)
  
    if(toupper){
      res <- toupper(res)
      if(verbose) message("Input nucleotides were capitalised")
    } 
    
  # Getting messages on the imported lines
    if(verbose) message("Number of read lines: ", length(res),"\n")
  
  # Check if the Fasta file is alternating, one line label, the next line sequence
    greplRes <- grepl(">",res)
    if(verbose) message("Number of header lines:", sum(greplRes),"\n")
  
  # Test if header and content rows are alternating
    is_odd_sequence <- function(x) {
      diffs <- diff(x)
      all(diffs == 2)
    }
    
    alternatingRows <- FALSE
    if( is_odd_sequence(which(greplRes))) {
      if(verbose) message("It seems your fasta file has alternating header/sequence rows")
      alternatingRows <- TRUE
    } else {
      if(verbose) message("It seems your fasta file does not alternating header/sequence rows")
    }
  
  # Now populate the output vector with the sequences
    if(alternatingRows){
    # Sequences are alternating, hence we can use this to quickly import the files
      seq <- res[seq(2,length(res),2)]
      names(seq) <- res[seq(1,length(res)-1,2)]
    } else {
      idRows <- which(greplRes)
      numberOfSequences <- length(idRows)
      
      # NOTE: Quick and dirty for now with a loop, fix that to be faster later!!!
      seq <- rep("", numberOfSequences)
      for(i in 1:(numberOfSequences-1)){
        seq[i] <- paste(res[(idRows[i]+1):(idRows[i+1]-1)], collapse="")
      }
      
      seq[numberOfSequences] <- paste(res[(idRows[i+1]+1):(length(res))], collapse="")
      names(seq) <- res[idRows]
    }
  class(seq) <- "fa"
  seq
}