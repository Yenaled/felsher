# Usage: Rscript merge_bed_with_genes.r <input_bed_file> <output_file> <genes_list_file>

args <- commandArgs(TRUE)
bed_file <- args[1]
output_file <- args[2]
genes_files <- args[3:length(args)]

bed_data <- read.table(bed_file, sep="\t", header=FALSE, stringsAsFactors=FALSE, fill=TRUE)
genes_data <- c()
for (f in genes_files) {
    genes_data <- c(genes_data, readLines(f))
}
genes_data <- unique(genes_data)
bed_data <- bed_data[bed_data$V4 %in% genes_data,]
write.table(bed_data, file=output_file, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
