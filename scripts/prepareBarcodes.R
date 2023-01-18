args <- commandArgs(trailingOnly=TRUE)
samplesheet.file <- args[1]
#message("I will use the sample sheet file: ", samplesheet.file)
#samplesheet.file <- "/path/to/the/file"
samplesheet <- read.table(samplesheet.file, header=TRUE, sep="\t")

mockSamples <- NULL
mockSamples <- unique(samplesheet$sample_name[samplesheet$useForMock=="YES"])

if(length(mockSamples)==0) system("echo 'ERROR!!! No samples for mock reference generation found. Make sure to use YES/NO in the sample sheet column 5 (useForMock)'")

samplesheet$useForMock[is.element(samplesheet$sample_name, mockSamples)] <- "YES"

samplesheet <- cbind("NNNNNN+NNNNNN",samplesheet)
colnames(samplesheet)[1] <- "barcode"

barcodes <- unique(samplesheet[,c("barcode", "sample_name", "useForMock")])

barcodes.file <- args[2]

write.table(barcodes, file=barcodes.file, quote=FALSE, row.names = FALSE, col.names = FALSE, sep="\t")