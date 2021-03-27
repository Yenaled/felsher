###############################################################
### File: analyze_signatures.r
### Description: Performs analysis of MYC signatures
### Usage: Rscript analyze_signatures2.r
### Written by Delaney Sullivan
###############################################################

require(openxlsx)
require(Rtsne)
require(UpSetR)

signature_dir <- "./output/integrative_signature/"
output_dir <- "./output/signature_analysis/"

# Signatures
signature <- read.csv(paste(signature_dir, "signature.csv", sep=""), header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, row.names=1)
signature_myc_cor <- read.csv(paste(signature_dir, "signature_myc_correlation.csv", sep=""), header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, row.names=1)
signature_tumorigenesis <- read.csv(paste(signature_dir, "signature_tumorigenesis.csv", sep=""), header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, row.names=1)

# Start logging
sink(paste(output_dir, "log.txt", sep=""))

# CCLE data file
ccle_data_file <- "./ccle/CCLE_RNAseq_rsem_genes_tpm.deseq.log2.tsv.gz"
ccle_data <- read.table(ccle_data_file, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
rownames(ccle_data) <- ccle_data$gene_id
ccle_data$gene_id <- NULL
ccle_data$transcript_ids <- NULL
rownames(ccle_data) <- gsub("\\..*", "", rownames(ccle_data))
  
# CCLE correlations
ccle_cor_file <- "./ccle/CCLE_myc_cor.csv"
ccle_cor <- read.csv(ccle_cor_file, header=TRUE, stringsAsFactors=FALSE, row.names=1)
rownames(ccle_cor) <- gsub("\\..*", "", rownames(ccle_cor))

# TCGA correlation data
tcga_cor_data <- read.csv(paste(signature_dir, "correlation_data.csv", sep=""), header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, row.names=1)
# TCGA LFC tumor-normal
tcga_lfc_file <- "./tcga/TCGA.processed.tumor_normal_lfc.tsv.gz"
tcga_lfc <- read.table(tcga_lfc_file, header=TRUE, stringsAsFactors=FALSE, row.names=1, check.names=FALSE)
rownames(tcga_lfc) <- gsub("\\..*", "", rownames(tcga_lfc))

# CCLE MYC clustering: plot PCA/t-SNE
data_signature <- t(ccle_data[rownames(ccle_data) %in% c("ENSG00000136997",rownames(signature)),])
data_signature <- as.data.frame(data_signature)
myc_raw <- data_signature[,"ENSG00000136997"]
myc <- scale(data_signature[,"ENSG00000136997"], center = TRUE, scale = TRUE)
myc_rank <- rank(myc)
data_signature$ENSG00000136997 <- NULL
lower_tertile_size <- ceiling(length(myc_rank) / 3)
mid_tertile_size <- floor(length(myc_rank) / 3)
upper_tertile_size <- length(myc_rank) - mid_tertile_size - lower_tertile_size
colors <- c(rep("blue", lower_tertile_size), rep("gray", mid_tertile_size), rep("red", upper_tertile_size))
colors <- colors[myc_rank]
set.seed(42)
tsne <- Rtsne(scale(data_signature), dims = 2, perplexity=20, verbose=FALSE, max_iter = 5000, initial_dims=50)
pdf(paste(output_dir, "ccle_clustering.pdf", sep=""))
par(cex.axis=1.75, cex.lab=1.75, mar=c(5,5,5,3))
plot(tsne$Y, col=colors, pch=19, xlab="t-SNE1", ylab="t-SNE2")
axis(side = 1, lwd.ticks=2.25)
axis(side = 2, lwd.ticks=2.25)
box(lwd=2.25)
device <- dev.off()
write.csv(t(data_signature), file=paste(output_dir, "ccle_clustering.csv", sep=""))

# TCGA LFC analysis
rowMedians <- function(df) {
  return(apply(df, 1, median, na.rm = T))
}
tcga_lfc_medians <- rowMedians(tcga_lfc)
tcga_lfc_medians <- data.frame(tcga_lfc_medians)
tcga_data <- merge(tcga_lfc_medians, tcga_cor_data, by.x=0, by.y=0, all=F)
scatterPlot <- function(res, pearson_cutoff, lfcthresh, main="Volcano Plot",
                        legendpos="bottomright", textcx=1, nsamples=0, ...) {
  with(res, plot(median_pearson, tcga_lfc_medians, xaxt='n', pch=20, cex.lab=1.3,
                 cex.axis=1.7, xlab="Median Pearson's r", ylab=paste("Log2 Fold Change (n = ", nsamples, " tumors with matched normal samples)", sep=""),
                 col="gray", xlim=c(-0.3, 0.60), ylim=c(-4.75,4), ...))
  with(subset(res, up & Freq == 4),
       points(median_pearson, tcga_lfc_medians, pch=20, col="purple", ...))
  with(subset(res, up & Freq == 5),
       points(median_pearson, tcga_lfc_medians, pch=20, col="red", ...))
  with(subset(res, down & Freq == 4),
       points(median_pearson, tcga_lfc_medians, pch=20, col="deepskyblue", ...))
  with(subset(res, down & Freq == 5),
       points(median_pearson, tcga_lfc_medians, pch=20, col="blue", ...))
  with(res, abline(lty=2, h=lfcthresh))
  with(res, abline(lty=2, v=pearson_cutoff))
  with(res, axis(side=1,cex.lab=1.5, cex.axis=2, at=c(-0.3, -0.15, 0, 0.15, 0.3, 0.45, 0.60)))
}
pdf(paste(output_dir, "tcga_lfc.pdf", sep=""))
scatterPlot(tcga_data, 0.3, 0, nsamples=ncol(tcga_lfc))
device <- dev.off()
rownames(tcga_data) <- tcga_data[,1]
tcga_data[,1] <- NULL
tcga_data <- tcga_data[,c("tcga_lfc_medians","median_pearson","Symbol")]
colnames(tcga_data) <- c("Median Log2 Fold Change Tumor versus Normal", "Median Pearson Correlation Over All Tumors", "Symbol")
write.csv(tcga_data, file=paste(output_dir, "tcga_lfc.csv", sep=""))
