###############################################################
### File: analyze_signatures.r
### Description: Performs comparison of MYC signatures
### Usage: Rscript analyze_signatures.r
### Written by Delaney Sullivan
###############################################################

require(openxlsx)
require(Rtsne)
require(UpSetR)
require(ggbeeswarm)

signature_dir <- "./output/integrative_signature/"
output_dir <- "./output/signature_analysis/"
genesets_myc_file <- "./genesets/myc_signature_genesets_ensembl.gmx"
geneset_names_file <- "./genesets/geneset_names.csv"
geneset_name_this_study <- "This Study"

# Load MYC gene sets
genesets_myc <- read.table(genesets_myc_file, stringsAsFactors=FALSE, check.names=FALSE, header=TRUE, sep="\t")
geneset_names <- read.csv(geneset_names_file, stringsAsFactors=FALSE, check.names=FALSE, header=TRUE)

# Signatures
signature <- read.csv(paste(signature_dir, "signature.csv", sep=""), header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, row.names=1)
signature_myc_cor <- read.csv(paste(signature_dir, "signature_myc_correlation.csv", sep=""), header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, row.names=1)
signature_tumorigenesis <- read.csv(paste(signature_dir, "signature_tumorigenesis.csv", sep=""), header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, row.names=1)

# Start logging
sink(paste(output_dir, "log.txt", sep=""))

# Print out number of genes in each MYC gene set
listInput <- list()
for (i in colnames(genesets_myc)) {
  curr_geneset <- genesets_myc[,i,drop=TRUE]
  curr_geneset <- curr_geneset[!is.na(curr_geneset)]
  curr_geneset <- unique(curr_geneset[curr_geneset != "" & curr_geneset != "NA"])
  print(paste(i, length(curr_geneset)), sep=" ")
  curr_geneset_name <- geneset_names[toupper(geneset_names[,"Name"]) == toupper(i),"Gene Set"]
  listInput[[curr_geneset_name]] <- curr_geneset
}
genes_in_myc_genesets <- unique(unname(unlist(listInput)))
genes_in_signature <- rownames(signature)
listInput[[geneset_name_this_study]] <- genes_in_signature
genes_not_in_myc_genesets <- signature[genes_in_signature[!(genes_in_signature %in% genes_in_myc_genesets)],"Symbol"]
print(paste("Genes in this study not present in other genesets: ", paste(genes_not_in_myc_genesets, collapse=","), sep=""))
# Order genesets by overlap with this study's signature:
geneset_ordering <- data.frame(row.names=names(listInput))
for (i in rownames(geneset_ordering)) {
  geneset_ordering[i,"Set Size"] <- length(listInput[[i]])
  geneset_ordering[i,"Overlap"] <- length(intersect(listInput[[i]], genes_in_signature))
}
geneset_ordering <- geneset_ordering[order(-geneset_ordering$Overlap),,drop=FALSE]
# Make upset plot:
pdf(paste(output_dir, "upset_genesets.pdf", sep=""), height=7, width=32, onefile=FALSE)
upset(fromList(listInput), order.by = "freq", sets=rev(rownames(geneset_ordering)), nintersects=NA, keep.order = TRUE)
device <- dev.off()
pdf(paste(output_dir, "upset_genesets_top_5.pdf", sep=""), height=7, width=8, onefile=FALSE)
upset(fromList(listInput), text.scale=1.7, order.by = "freq", sets=rev(rownames(geneset_ordering)[1:5]), nintersects=NA, keep.order = TRUE)
device <- dev.off()
# Write table describing the gene sets
output_table <- merge(geneset_names, geneset_ordering, by.x="Gene Set", by.y=0)
output_table[output_table$ID == "","Name"] <- ""
output_table <- output_table[order(-output_table$Overlap),]
colnames(output_table) <- c("Gene Set", "MSigDB Name", "MSigDB ID", "Set Size", "Size of Intersection with This Study's Signature")
write.csv(output_table, file=paste(output_dir, "genesets_info.csv", sep=""), row.names=FALSE)

# TCGA meta-z scores
tcga_z_scores_file <- "./tcga/tcga_zscores.csv"
if (!file.exists(tcga_z_scores_file)) {
  tcga_z_scores_file <- "./tcga/tcga_zscores.csv.gz"
}
tcga_z_scores <- read.csv(tcga_z_scores_file, header=TRUE, stringsAsFactors=FALSE, row.names=1, check.names=FALSE)

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

# TCGA meta-z score analysis: plot ECDF
signature.metaz <- tcga_z_scores[rownames(tcga_z_scores) %in% rownames(signature),"MetaZ",drop=FALSE]
pdf(paste(output_dir, "tcga_metaz.pdf", sep=""))
par(cex=1.5, cex.lab=1.2)
plot(ecdf(as.numeric(tcga_z_scores[,"MetaZ"])), pch=".", cex=2.5, lwd=3, col="black", xlim = range(c(-8, 8)), xaxt = "n", yaxt="n", xlab="Meta-z score", ylab="ECDF", main=NULL)
plot(ecdf(as.numeric(signature.metaz[,"MetaZ"])), col="red", pch=".", cex=2.5, lwd=3.5, add=TRUE, verticals=TRUE, xlim = range(c(-8, 8)), xaxt = "n", yaxt="n")
axis(side = 1, at = c(-7.5, -5, -2.5, 0, 2.5, 5, 7.5), lwd.ticks=2.25)
axis(side = 2, at = c(0,0.25,0.5,0.75,1), lwd.ticks=2.25)
box(lwd=2.25)
device <- dev.off()
print("Meta-z score analysis")
print(signature.metaz[1:5,,drop=FALSE])
print(wilcox.test(as.numeric(signature.metaz[,"MetaZ"]), as.numeric(tcga_z_scores[,"MetaZ"])))

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

# CCLE MYC correlation analysis (comparisons to other gene sets)
genesets.cor <- NULL
genesets.med <- list()
for (i in names(listInput)) {
  curr_geneset.cor <- ccle_cor[rownames(ccle_cor) %in% listInput[[i]],"r",drop=FALSE]
  curr_geneset.cor$group <- i
  curr_geneset.cor <- curr_geneset.cor[complete.cases(curr_geneset.cor),]
  genesets.med[[i]] <- median(curr_geneset.cor$r)
  if (is.null(genesets.cor)) {
    genesets.cor <- curr_geneset.cor
  } else {
    genesets.cor <- rbind(genesets.cor, curr_geneset.cor)
  }
}
genesets.order <- genesets.med[order(unlist(genesets.med))]
curr_geneset.cor <- ccle_cor[!(rownames(ccle_cor) %in% "ENSG00000136997"),"r",drop=FALSE] # Reference distribution (all genes except MYC)
curr_geneset.cor$group <- "Reference Distribution"
curr_geneset.cor <- curr_geneset.cor[complete.cases(curr_geneset.cor),]
genesets.cor <- rbind(genesets.cor, curr_geneset.cor)
genesets.cor$group <- factor(genesets.cor$group, levels = c("Reference Distribution", names(genesets.order)))
p <- ggplot(genesets.cor, aes(group, r)) + geom_boxplot(outlier.shape=NA, coef=0) + geom_beeswarm(data=genesets.cor[genesets.cor$group != "Reference Distribution",], show.legend=FALSE,cex=0.8,size=0.8,colour="#666666",alpha=0.8)
p <- p + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))
p <- p + ylab("Pearson's r")
p <- p + theme(axis.title.x=element_blank(), plot.margin=margin(t = 2, r = 6, b = 0, l = 36, unit = "pt"), axis.text = element_text(size = 22, colour="black"), axis.title.y=element_text(size = 22, colour="black"))
pdf(paste(output_dir, "ccle_genesets_cor.pdf", sep=""), height=7, width=23)
print(p)
device <- dev.off()

# Meta-z analysis (comparisons to other gene sets)
genesets.z <- NULL
genesets.med <- list()
for (i in names(listInput)) {
  curr_geneset.z <- tcga_z_scores[rownames(tcga_z_scores) %in% listInput[[i]],"MetaZ",drop=FALSE]
  curr_geneset.z$group <- i
  curr_geneset.z <- curr_geneset.z[complete.cases(curr_geneset.z),]
  genesets.med[[i]] <- median(curr_geneset.z$MetaZ)
  if (is.null(genesets.z)) {
    genesets.z <- curr_geneset.z
  } else {
    genesets.z <- rbind(genesets.z, curr_geneset.z)
  }
}
genesets.order <- genesets.med[order(unlist(genesets.med))]
curr_geneset.z <- tcga_z_scores[!(rownames(tcga_z_scores) %in% "ENSG00000136997"),"MetaZ",drop=FALSE] # Reference distribution (all genes except MYC)
curr_geneset.z$group <- "Reference Distribution"
curr_geneset.z <- curr_geneset.z[complete.cases(curr_geneset.z),]
genesets.z <- rbind(genesets.z, curr_geneset.z)
genesets.z$group <- factor(genesets.z$group, levels = c("Reference Distribution", names(genesets.order)))
p <- ggplot(genesets.z, aes(group, MetaZ)) + geom_boxplot(outlier.shape=NA, coef=0) + geom_beeswarm(data=genesets.z[genesets.z$group != "Reference Distribution",], show.legend=FALSE,cex=0.6,size=0.9,colour="#666666",alpha=0.8)
p <- p + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))
p <- p + ylab("Meta-z score")
p <- p + theme(axis.title.x=element_blank(), plot.margin=margin(t = 2, r = 6, b = 0, l = 64, unit = "pt"), axis.text = element_text(size = 22, colour="black"), axis.title.y=element_text(size = 22, colour="black"))
pdf(paste(output_dir, "tcga_genesets_metaz.pdf", sep=""), height=7, width=23)
print(p)
device <- dev.off()


