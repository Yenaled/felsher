# Usage: Rscript ccle_pairwise_correlation.r <gene_to_correlate> <input_file> <output_file>
# Example: Rscript ccle_pairwise_correlation.r ENSG00000136997 CCLE_RNAseq_rsem_genes_tpm.deseq.log2.tsv.gz CCLE_myc_cor.csv

args <- commandArgs(TRUE)
gene <- args[1]
input_file <- args[2]
output_file <- args[3]

data <- read.table(input_file, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
rownames(data) <- data$gene_id
data$gene_id <- NULL
data$transcript_ids <- NULL
data.myc <- unlist(data["ENSG00000136997",]) # The MYC row
cors <- apply(data, 1, y=data.myc, function(x, y) cor.test(y,x)$estimate ) # This might spit out a bunch of warnings for cases where std dev of x is zero
cors <- data.frame(cors)
colnames(cors) <- "r"
cors <- cors[order(cors, decreasing=TRUE),,drop=FALSE]
write.csv(cors, file=output_file, quote=FALSE)
