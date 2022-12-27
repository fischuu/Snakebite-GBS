library("SimRAD")
library("GenomicTools")

minLength <- as.numeric(minLength)-1
maxLength <- as.numeric(maxLength)+1

# Adjust the paths
if(substr(refGenome.file,1,1)!="/") refGenome.file <- file.path(projFolder, refGenome.file)

# The SimRAD default import function for fasta sequences fails for larger genomes, so I use an own one here
#fasta <- ref.DNAseq(refGenome.file, subselect.contigs = FALSE)
fasta <- importFA(refGenome.file)
names(fasta) <- NULL

enz1 <- strsplit(enz1.in, "")
enz2 <- strsplit(enz2.in, "")

enz1.1 <- enz1[[1]][1]
enz1.2 <- paste(enz1[[1]][-1], collapse="")
enz2.1 <- enz2[[1]][1]
enz2.2 <- paste(enz2[[1]][-1], collapse="")

ddout <- insilico.digest(fasta, cut_site_5prime1 = enz1.1, 
                         cut_site_3prime1 = enz1.2,
                         cut_site_5prime2 = enz2.1,
                         cut_site_3prime2 = enz2.2)

names(ddout) <- paste0("> Location", 1:length(ddout))
ddout.selected <- size.select(ddout, min.size = minLength, max.size = maxLength, graph = FALSE)

exportFA(ddout, file=full.file)
exportFA(ddout.selected, file=selected.file)
