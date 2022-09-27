# SCript to visualise the size distribution
library("data.table")
projFolder <- "/scratch/project_2003491/ArctAquaTest2/"

projFolder <- "/scratch/project_2003491/SLU_Arctic_charr/"

mm.file <- "MPILEUP/mpileup_reference"

mm <- fread(file.path(projFolder, mm.file, "GSC.MasterMatrix.txt"), select=1:2, sep="\t", nrows = 100000)
#mm <- fread(file.path(projFolder, mm.file, "GSC.MasterMatrix.txt"), select=1:2, sep="\t")

steps <- mm$V2[2:length(mm$V2)] - mm$V2[1:(length(mm$V2)-1)]

# Fixing chromosome steps
# Of course, you could do the analysis also for each chromosome, but I think there should not be a chromosome effect
steps[steps<=0] <- 1


hist(log(steps),prob=TRUE)
step.dens <- density(log(steps))

lines(step.dens, col="red", lty=2)
v <- optimize(approxfun(step.dens$x,step.dens$y),interval=c(1,12))$minimum
abline(v=v, col="blue")
exp(v)
estimated
