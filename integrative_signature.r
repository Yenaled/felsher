###############################################################
### File: integrative_signature.r
### Description: Integrates mouse experimental and human TCGA 
###              MYC signatures
### Usage: Rscript integrative_signature.r
### Written by Delaney Sullivan
###############################################################

source("analysis_setup.r")
output_dir <- "./output/integrative_signature/"

mouse_studies <- names(sample_mapping_myc)
overlap_cutoff <- 4
pearson_cutoff <- 0.3
fdr_cutoff <- 0.05

sink(paste(output_dir, "log.txt", sep=""))

print(paste("Selecting mouse genes differentially expressed in at least", overlap_cutoff, "out of", length(mouse_studies), "studies"))

mouse_de_up_data <- genes_de_up_ids[,mouse_studies]
mouse_de_up <- table(unlist(mouse_de_up_data))
mouse_de_up <- data.frame(mouse_de_up[mouse_de_up >= overlap_cutoff])
colnames(mouse_de_up) <- c("Mouse", "Freq")
mouse_de_up <- merge(orthologs_ids, mouse_de_up, by="Mouse", all.x=FALSE, all.y=TRUE)
mouse_de_up <- merge(ensembl_human_mapping, mouse_de_up, by.x="ID", by.y="Human", all.x=FALSE, all.y=TRUE)
if (nrow(mouse_de_up[is.na(mouse_de_up$Ensembl),]) > 0) {
  print(paste("Following mouse up-gene IDs can't be mapped to human ENSEMBL IDs:", paste(mouse_de_up[is.na(mouse_de_up$Ensembl),"Mouse"], collapse=",")))
}

mouse_de_down_data <- genes_de_down_ids[,mouse_studies]
mouse_de_down <- table(unlist(mouse_de_down_data))
mouse_de_down <- data.frame(mouse_de_down[mouse_de_down >= overlap_cutoff])
colnames(mouse_de_down) <- c("Mouse", "Freq")
mouse_de_down <- merge(orthologs_ids, mouse_de_down, by="Mouse", all.x=FALSE, all.y=TRUE)
mouse_de_down <- merge(ensembl_human_mapping, mouse_de_down, by.x="ID", by.y="Human", all.x=FALSE, all.y=TRUE)
if (nrow(mouse_de_down[is.na(mouse_de_down$Ensembl),]) > 0) {
  print(paste("Following mouse down-gene IDs can't be mapped to human ENSEMBL IDs:", paste(mouse_de_down[is.na(mouse_de_down$Ensembl),"Mouse"], collapse=",")))
}

tcga_rra_data <- read.csv(tcga_rra_file, stringsAsFactors=FALSE, row.names=1)
rownames(tcga_rra_data) <- gsub("\\..*", "", rownames(tcga_rra_data))
tcga_rra_data$up <- FALSE
tcga_rra_data$down <- FALSE
tcga_rra_data[rownames(tcga_rra_data) %in% mouse_de_up$Ensembl,"up"] <- TRUE
tcga_rra_data[rownames(tcga_rra_data) %in% mouse_de_down$Ensembl,"down"] <- TRUE
tcga_rra_data <- merge(tcga_rra_data, rbind(mouse_de_up, mouse_de_down)[,c("Ensembl","Freq")], by.x="row.names", by.y="Ensembl", all.x=TRUE, all.y=FALSE)
rownames(tcga_rra_data) <- tcga_rra_data[,1]
tcga_rra_data[,1] <- NULL

if (nrow(tcga_rra_data[tcga_rra_data$up,]) != nrow(mouse_de_up[!is.na(mouse_de_up$Ensembl),])) {
  aa <- mouse_de_up[!is.na(mouse_de_up$Ensembl),"Ensembl"]
  bb <- rownames(tcga_rra_data[tcga_rra_data$up,])
  print(paste("Up-Gene ENSEMBL IDs", paste(aa[!(aa %in% bb)], collapse=","), "not found in TCGA data file"))
}
if (nrow(tcga_rra_data[tcga_rra_data$down,]) != nrow(mouse_de_down[!is.na(mouse_de_down$Ensembl),])) {
  aa <- mouse_de_down[!is.na(mouse_de_down$Ensembl),"Ensembl"]
  bb <- rownames(tcga_rra_data[tcga_rra_data$down,])
  print(paste("Down-Gene ENSEMBL IDs", paste(aa[!(aa %in% bb)], collapse=","), "not found in TCGA data file"))
}

# Volcano Plot
volcanoplot <- function(res, pearson_cutoff, sigthresh=0.05, main="Volcano Plot",
                         legendpos="bottomright", textcx=1, ...) {
  res$pvalue <- res$padj
  res[!is.na(res$pvalue) & res$pvalue < 10^-40, "pvalue"] <- 10^-40
  with(res, plot(median_pearson, -log10(pvalue), xaxt='n', pch=20, cex.lab=1.3,
                 cex.axis=1.7, xlab="Median Pearson's r", ylab="-Log10(FDR-adjusted p-value)",
                 col="gray", xlim=c(-0.3, 0.60), ylim=c(0,40), ...))
  with(subset(res, up & Freq == 4),
       points(median_pearson, -log10(pvalue), pch=20, col="purple", ...))
  with(subset(res, up & Freq == 5),
       points(median_pearson, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, down & Freq == 4),
       points(median_pearson, -log10(pvalue), pch=20, col="deepskyblue", ...))
  with(subset(res, down & Freq == 5),
       points(median_pearson, -log10(pvalue), pch=20, col="blue", ...))
  with(res, abline(lty=2, h=-log10(sigthresh)))
  with(res, abline(lty=2, v=pearson_cutoff))
  with(res, axis(side=1,cex.lab=1.5, cex.axis=2, at=c(-0.3, -0.15, 0, 0.15, 0.3, 0.45, 0.60)))
}
pdf(paste(output_dir, "volcano.pdf", sep=""))
volcanoplot(tcga_rra_data, pearson_cutoff, fdr_cutoff)
device <- dev.off()

tcga_rra_data$MYC_Correlated <- FALSE
tcga_rra_data[tcga_rra_data$padj < fdr_cutoff & tcga_rra_data$median_pearson > pearson_cutoff,"MYC_Correlated"] <- TRUE
tcga_rra_data <- merge(tcga_rra_data, ensembl_human_symbols_mapping, by.x="row.names", by.y="Ensembl", all.x=TRUE, all.y=FALSE)
rownames(tcga_rra_data) <- tcga_rra_data[,1]
tcga_rra_data[,1] <- NULL
tcga_rra_data <- tcga_rra_data[order(-tcga_rra_data$median_pearson),]

print("Genes with highest median correlation:")
print(tcga_rra_data[1:25,c("Symbol","padj","median_pearson","up")])
print(tcga_rra_data[tcga_rra_data$down & tcga_rra_data$Freq == 5,c("Symbol","padj","median_pearson","up")])

# Write out TCGA correlation data
write.csv(tcga_rra_data, file=paste(output_dir, "correlation_data.csv", sep=""))

# Write out formatted TCGA correlation data
write.csv(tcga_rra_data[,c("Symbol","pval","padj","median_pearson")], file=paste(output_dir, "correlation_data_formatted.csv", sep=""))

# Write out final signature
final_signature <- tcga_rra_data[tcga_rra_data$MYC_Correlated & tcga_rra_data$up,]
final_signature <- final_signature[,c("Symbol","pval","padj","median_pearson","Freq")]
colnames(final_signature) <- c("Symbol","P-Value","Adjusted P-Value","Median Pearson","Number of MYC mouse experiments in which the gene is differentially expressed")
write.csv(final_signature, file=paste(output_dir, "signature.csv", sep=""))

# Write out other signatures
write.csv(tcga_rra_data[!(rownames(tcga_rra_data) %in% rownames(final_signature)) & tcga_rra_data$up,"Symbol",drop=FALSE], file=paste(output_dir, "signature_tumorigenesis.csv", sep=""))
write.csv(tcga_rra_data[!(rownames(tcga_rra_data) %in% rownames(final_signature)) & tcga_rra_data$MYC_Correlated,"Symbol",drop=FALSE], file=paste(output_dir, "signature_myc_correlation.csv", sep=""))

# Print out info about number of genes in signatures
print("Number of genes in signatures")
num_myc_cor_genes <- nrow(tcga_rra_data[tcga_rra_data$MYC_Correlated,])
num_signature_genes <- nrow(final_signature)
num_mouse_de_genes <- nrow(mouse_de_up)
num_mouse_de_genes_mapped_to_human <- nrow(tcga_rra_data[tcga_rra_data$up,])
print(paste("Number of MYC-correlated genes:", num_myc_cor_genes))
print(paste("Number of MYC-correlated (but not mouse-DE) genes:", num_myc_cor_genes - num_signature_genes))
print(paste("Number of MYC-correlated AND mouse-DE genes:", num_signature_genes))
print(paste("Number of total mouse-DE (but not MYC-correlated) genes:", num_mouse_de_genes - num_signature_genes))
print(paste("Number of human-ortholog-mapped mouse-DE (but not MYC-correlated) genes:", num_mouse_de_genes_mapped_to_human - num_signature_genes))
print(paste("Number of total mouse-DE genes:", num_mouse_de_genes))
print(paste("Number of mouse-DE genes that are human-ortholog-mapped:", num_mouse_de_genes_mapped_to_human))
print(paste("Number of mouse-DE genes without human ortholog:", num_mouse_de_genes - num_mouse_de_genes_mapped_to_human))

sink()


# EnrichR (GO Biological Process + Mouse Gene Atlas) for final signature
signature_go_bp <- do_enrichment("GO_Biological_Process_2018", NULL, list(Signature=final_signature$Symbol), output_file=paste(output_dir, "signature_go_bp.xlsx", sep=""))
signature_go_mouse_cells <- do_enrichment("Mouse_Gene_Atlas", NULL, list(Signature=final_signature$Symbol), output_file=paste(output_dir, "signature_go_mouse_cells.xlsx", sep=""))

# EnrichR (GO Biological Process) for "tumorigenesis" signature (i.e. upregulated in mouse models but doesn't meet correlation threshold in TCGA)
signature_tumorigenesis_go_bp <- do_enrichment("GO_Biological_Process_2018", NULL, list(Signature=tcga_rra_data[!(rownames(tcga_rra_data) %in% rownames(final_signature)) & tcga_rra_data$up,"Symbol"]), output_file=paste(output_dir, "signature_tumorigenesis_go_bp.xlsx", sep=""))

# GO Biological Process pathway analysis for the top 5 pathways in the final signature and the tumorigenesis signature

num_pathways <- 5
mouse_genes_up <- tcga_rra_data[tcga_rra_data$up,c("Symbol","median_pearson")]
mouse_genes_up <- mouse_genes_up[order(mouse_genes_up$median_pearson),]
top_pathways_signature <- rownames(signature_go_bp[1:num_pathways,])
top_pathways_signature_tumorigenesis <- rownames(signature_tumorigenesis_go_bp[1:num_pathways,])

# Load Gene Ontology genesets
raw_genesets_go <- readLines("data/GO_Biological_Process_2018.txt")
genesets_go <- list()
for (line in raw_genesets_go) {
    for (geneset_name in c(top_pathways_signature, top_pathways_signature_tumorigenesis)) {
        if (startsWith(line, paste(geneset_name, "\t", sep=""))) {
            geneset_data <- unlist(strsplit(line, "\t"))
            geneset_data <- c(geneset_name, gsub("(.*),.*", "\\1", geneset_data[2:length(geneset_data)]))
            genesets_go[[geneset_name]] <- geneset_data
        }
    }
}
genelist_top_pathways_signature <- c()
for (pathway in top_pathways_signature) {
    genelist_top_pathways_signature <- c(genelist_top_pathways_signature, genesets_go[[pathway]])
}
genelist_top_pathways_signature <- unique(genelist_top_pathways_signature)
genelist_top_pathways_signature_tumorigenesis <- c()
for (pathway in top_pathways_signature_tumorigenesis) {
    genelist_top_pathways_signature_tumorigenesis <- c(genelist_top_pathways_signature_tumorigenesis, genesets_go[[pathway]])
}
genelist_top_pathways_signature_tumorigenesis <- unique(genelist_top_pathways_signature_tumorigenesis)

# Output barcode plots of top pathways
pdf(paste(output_dir, "barcode_top_", num_pathways, "_pathways.pdf", sep=""), width=4.5, height=2)
pathways <- list(genelist_top_pathways_signature, genelist_top_pathways_signature_tumorigenesis)
par(mfrow=c(length(pathways),1), mar=c(0.5,0.5,0.5,0.5))
for (pathway in pathways) {
    inPathway <- toupper(mouse_genes_up$Symbol) %in% toupper(pathway)
    barplot(height=inPathway, col="black", border=NA, space=0, axes=FALSE)
    segments(0,0, 0,1)
    segments(0,1, length(inPathway),1)
    segments(length(inPathway),1, length(inPathway),0)
    segments(0,0, length(inPathway),0)
    abline(v = num_below_cor, col="red")
}
device <- dev.off()
