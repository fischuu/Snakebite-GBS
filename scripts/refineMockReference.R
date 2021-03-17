# Load required libraries
  library("GenomicTools")

# Import the coverage
  clusterFiles <- list.files(file.path(projFolder, "FASTQ", "TRIMMED", "alignments_clusters"), pattern="*.coverage")
  clusterCoverage <- read.table(file.path(projFolder, "FASTQ", "TRIMMED", "alignments_clusters", clusterFiles[1]))
  names(clusterCoverage)[1:2] <- c("cluster", clusterFiles[1])
  
  for(i in 2:length(clusterFiles)){
    tmp <- read.table(file.path(projFolder, "FASTQ", "TRIMMED", "alignments_clusters", clusterFiles[i]))
    names(tmp)[1:2] <- c("cluster", clusterFiles[i])
    clusterCoverage <- merge(clusterCoverage, tmp, by="cluster")
  }

# Get the total coverage
  totalCoverage <- apply(clusterCoverage[,-1],1,sum)

# Filter the coverage
  clusterCoverage.refined <- clusterCoverage[totalCoverage>minTotalReadCoverage,]
  clusterCoverage.refined <- clusterCoverage.refined[apply(clusterCoverage.refined[,-1]>0,1,sum)>=minSampleCoverage,]

# Import the initial mock clusters
  mockClusters <- importFA(file.path(mockClusters.file))
  
# Refine the mock clusters
  mockClusters.refined <- mockClusters[is.element(gsub(">","",names(mockClusters)), clusterCoverage.refined[,1])]
  
# Export the refined mock reference
  exportFA(mockClusters.refined, file=mockClusters.refined.file)
  