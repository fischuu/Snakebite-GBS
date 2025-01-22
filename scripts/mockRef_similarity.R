options(scipen=999)

# First, handing the mock to ref stats file
  mtrs.file <- file.path(projFolder, "BAM/Mockref/mockToRef.sam.stats")
  mtrs <- readLines(mtrs.file)
  
  mtrs.reads_mapped <- as.numeric(strsplit(mtrs[grep("reads mapped", mtrs)[1]], "\t")[[1]][3])
  mtrs.reads_unmapped <- as.numeric(strsplit(mtrs[grep("reads unmapped", mtrs)[1]], "\t")[[1]][3])
  
  mtrs.based_mapped <- as.numeric(strsplit(mtrs[grep("bases mapped", mtrs)[1]], "\t")[[1]][3])
  mtrs.bases_mapped_cigar <- as.numeric(strsplit(mtrs[grep("bases mapped \\(cigar", mtrs)[1]], "\t")[[1]][3])

  mtrs.mismatches <- as.numeric(strsplit(mtrs[grep("mismatches", mtrs)[1]], "\t")[[1]][3])

# Then the secondary alignments  
  mtrsf.file <- file.path(projFolder, "BAM/Mockref/mockToRef.sam.samflags")
  mtrsf <- read.table(mtrsf.file)
  
  mtrs.secondary_alignments <- mtrsf[3,1] + mtrsf[4,1]

# Then the secondary alignments  
  mtrfs.file <- file.path(projFolder, "BAM/Mockref/mockToRef.sam.flagstat")
  mtrfs <- readLines(mtrfs.file)
  
# Keep in mind, that this is primary and secondary alignments, so it is slightlt biased towards good,
# but for this purpose here it is alright.s
  mtrs.mapping_rate <- as.numeric(strsplit(strsplit(mtrfs[grep("mapped \\(", mtrfs)]," \\(")[[1]][2], "%")[[1]][1])
  
# NOW THE FINAL MOCK

# Second, handing the mock to ref stats file
  fmtrs.file <- file.path(projFolder, "BAM/FinalMockref/mockToRef.sam.stats")
  fmtrs <- readLines(fmtrs.file)
  
  fmtrs.reads_mapped <- as.numeric(strsplit(fmtrs[grep("reads mapped", fmtrs)[1]], "\t")[[1]][3])
  fmtrs.reads_unmapped <- as.numeric(strsplit(fmtrs[grep("reads unmapped", fmtrs)[1]], "\t")[[1]][3])
  
  fmtrs.based_mapped <- as.numeric(strsplit(fmtrs[grep("bases mapped", fmtrs)[1]], "\t")[[1]][3])
  fmtrs.bases_mapped_cigar <- as.numeric(strsplit(fmtrs[grep("bases mapped \\(cigar", fmtrs)[1]], "\t")[[1]][3])
  
  fmtrs.mismatches <- as.numeric(strsplit(fmtrs[grep("mismatches", fmtrs)[1]], "\t")[[1]][3])
  
# Then the secondary alignments  
  fmtrsf.file <- file.path(projFolder, "BAM/FinalMockref/mockToRef.sam.samflags")
  fmtrsf <- read.table(fmtrsf.file)
  
  fmtrs.secondary_alignments <- fmtrsf[3,1] + fmtrsf[4,1]
  
# Then the secondary alignments  
  fmtrfs.file <- file.path(projFolder, "BAM/FinalMockref/mockToRef.sam.flagstat")
  fmtrfs <- readLines(fmtrfs.file)
  
  # Keep in mind, that this is primary and secondary alignments, so it is slightlt biased towards good,
  # but for this purpose here it is alright.s
  fmtrs.mapping_rate <- as.numeric(strsplit(strsplit(fmtrfs[grep("mapped \\(", fmtrfs)]," \\(")[[1]][2], "%")[[1]][1])
  
# And now the alignment rates and so forth
  get_alignment_stats <- function(file){
    tmp <- readLines(file)
    getValue <- function(x){
      as.numeric(strsplit(x, " ")[[1]][1])
    }
    tmp.out <- data.frame(total = getValue(tmp[1]),
                          secondary = getValue(tmp[2]),
                          mapped = getValue(tmp[5]),
                          mapping_rate = getValue(tmp[5]) / getValue(tmp[1]),
                          secondary_rate = getValue(tmp[2]) / getValue(tmp[1]))
    tmp.out <- data.frame(variable= colnames(tmp.out), value = as.numeric(t(tmp.out)))
    tmp.out
  }
  
  mock_alignments.files <- list.files(file.path(projFolder, "FASTQ", "TRIMMED", "alignments_clusters"), pattern="*.flagstat", full.names = TRUE)  
  mock_alignments <- merge(get_alignment_stats(mock_alignments.files[1]),  
                           get_alignment_stats(mock_alignments.files[2]),
                           by="variable")
  
  for(i in 3:length(mock_alignments.files)){
    mock_alignments <- merge(mock_alignments,  
                             get_alignment_stats(mock_alignments.files[i]),
                             by="variable")
  }
  
  final_mock_alignments.files <- list.files(file.path(projFolder, "BAM", "alignments_finalMock"), pattern="*.flagstat", full.names = TRUE)  
  final_mock_alignments <- merge(get_alignment_stats(final_mock_alignments.files[1]),  
                                 get_alignment_stats(final_mock_alignments.files[2]),
                                 by="variable")
  
  for(i in 3:length(final_mock_alignments.files)){
    final_mock_alignments <- merge(final_mock_alignments,  
                                   get_alignment_stats(final_mock_alignments.files[i]),
                                   by="variable")
  }
  
# Now check the mocks vs insilicos
  mfifs.file <- file.path(projFolder, "BAM/MockVsInsilico/mockToFullInsilico.sam.flagstat")
  msifs.file <- file.path(projFolder, "BAM/MockVsInsilico/mockToSelectedInsilico.sam.flagstat")
  ffifs.file <- file.path(projFolder, "BAM/MockVsInsilico/finalToFullInsilico.sam.flagstat")
  fsifs.file <- file.path(projFolder, "BAM/MockVsInsilico/finalToSelectedInsilico.sam.flagstat")
  
  mfifs <- readLines(mfifs.file)
  msifs <- readLines(msifs.file)
  ffifs <- readLines(ffifs.file)
  fsifs <- readLines(fsifs.file)
  
  mfi.mapping_rate <- as.numeric(strsplit(strsplit(mfifs[grep("mapped \\(", mfifs)]," \\(")[[1]][2], "%")[[1]][1])
  msi.mapping_rate <- as.numeric(strsplit(strsplit(msifs[grep("mapped \\(", msifs)]," \\(")[[1]][2], "%")[[1]][1])
  ffi.mapping_rate <- as.numeric(strsplit(strsplit(ffifs[grep("mapped \\(", ffifs)]," \\(")[[1]][2], "%")[[1]][1])
  fsi.mapping_rate <- as.numeric(strsplit(strsplit(fsifs[grep("mapped \\(", fsifs)]," \\(")[[1]][2], "%")[[1]][1])
  
    output <- data.frame(mtrs.reads_mapped = mtrs.reads_mapped,
                       mtrs.secondary_alignments = mtrs.secondary_alignments,
                       mtrs.reads_unmapped = mtrs.reads_unmapped,
                       mtrs.based_mapped = mtrs.based_mapped,
                       mtrs.bases_mapped_cigar = mtrs.bases_mapped_cigar,
                       mtrs.mismatches = mtrs.mismatches,
                       mtrs.mapping_rate = mtrs.mapping_rate,
                       fmtrs.reads_mapped = fmtrs.reads_mapped,
                       fmtrs.secondary_alignments = fmtrs.secondary_alignments,
                       fmtrs.reads_unmapped = fmtrs.reads_unmapped,
                       fmtrs.based_mapped = fmtrs.based_mapped,
                       fmtrs.bases_mapped_cigar = fmtrs.bases_mapped_cigar,
                       fmtrs.mismatches = fmtrs.mismatches,
                       fmtrs.mapping_rate = fmtrs.mapping_rate,
                       sampleMock.alignment_rate = mean(unlist(mock_alignments[2,-1])),
                       sampleMock.secondary_rate = mean(unlist(mock_alignments[4,-1])),
                       sampleFinalMock.alignment_rate = mean(unlist(final_mock_alignments[2,-1])),
                       sampleFinalMock.secondary_rate = mean(unlist(final_mock_alignments[4,-1])),
                       mfi.mapping_rate = mfi.mapping_rate,
                       msi.mapping_rate = msi.mapping_rate,
                       ffi.mapping_rate = ffi.mapping_rate,
                       fsi.mapping_rate = fsi.mapping_rate)
  
  output <- data.frame(variable= colnames(output), value = as.numeric(t(output)))
  
  metric <-
# Ratio of secondary alignments in mock to ref alignment
  (1 - output[2,2] / output[1,2]) *
# Ratio of missingness
  (1 - output[3,2] / output[1,2]) *  
# Clipping ratio  
  (1 - output[5,2] / output[4,2])
  
metric
  
  write.table(output, output.report.file, sep = "\t", quote=FALSE, row.names = FALSE)
  