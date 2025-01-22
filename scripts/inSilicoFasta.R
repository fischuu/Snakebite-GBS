library("SimRAD")
library("GenomicTools")

source(file.path(pipeFolder, "scripts", "modified_digestion_methods.R"))

minLength <- as.numeric(minLength)-1
maxLength <- as.numeric(maxLength)+1

# Adjust the paths
if(substr(refGenome.file,1,1)!="/") refGenome.file <- file.path(projFolder, refGenome.file)

# The SimRAD default import function for fasta sequences fails for larger genomes, so I use an own one here
#fasta <- ref.DNAseq(refGenome.file, subselect.contigs = FALSE)
fasta <- importFA(refGenome.file)
names(fasta) <- NULL

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

# Remember, this function is based on an updated version of the original SimRAD package, and is imported
# with the source command earlier in this script
ddout <- insilico.dddigest(fasta, cut_site_5prime1 = enz1.1, 
                           cut_site_3prime1 = enz1.2,
                           cut_site_5prime2 = enz2.1,
                           cut_site_3prime2 = enz2.2)

names(ddout) <- paste0("> Location", 1:length(ddout))

refseq.selA <- myadapt.ABBAselect(ddout, enz1.1, enz1.2, enz2.1, enz2.2) 

ddout.selected <- size.select(refseq.selA, min.size = minLength, max.size = maxLength, graph = FALSE)
ddout.selected.fasta <- as.character(ddout.selected)
names(ddout.selected.fasta) <- ddout.selected@ranges@NAMES

write.table(refseq.selA@ranges@NAMES, file=file.path(projFolder, "References", "AB_contigs.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
exportFA(ddout, file=full.file)
exportFA(ddout.selected.fasta, file=selected.file)
