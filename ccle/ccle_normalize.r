# Usage: Rscript ccle_normalize.r <input_file> <output_file>

require(sleuth)

args <- commandArgs(TRUE)
input_file <- args[1]
output_file <- args[2]

data <- read.table(input_file, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
normalization_factors <- norm_factors(data[,3:ncol(data)])
data_norm <- t(t(data[,3:ncol(data)]) / normalization_factors)
data_norm <- as.data.frame(log2(data_norm+1))
data_norm <- cbind(data[,c(1,2)], data_norm)
write.table(data_norm, file=output_file, sep="\t", quote=FALSE, row.names=FALSE)
