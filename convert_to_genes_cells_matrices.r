library(Seurat)
args = commandArgs(trailingOnly=TRUE)
sample <- paste0(args[1])
data <- Read10X(data.dir = sample)
write.csv(data, paste0(args[1],"/matrix",".csv"), quote=F)

#dir <- "/home/msnaveed/sra_local_repo/starsolo_results/AllReads/testis/analysis"

#for (exp in c(3,4,"5.1h","5.1m","5.2h","5.2m","5.3h","5.3m")){
#	sample <- paste0(dir,"/exp",exp)
#	data <- Read10X(data.dir = sample)
#	write.csv(data, paste0(dir,"/exp",exp,"/exp",exp,".csv"), quote=F)
#}
