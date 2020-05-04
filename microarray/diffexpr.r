###############################################################
### File: diffexpr.r
### Description: Gets differential gene expression
###              from microarray data. Outputs into 
###              supplied <input_directory>/tissue/.
### Usage: Rscript diffexpr.r <input_directory> <annotation_file>
### Written by Delaney Sullivan
###############################################################

# Get command-line arguments
args <- commandArgs(TRUE)
input_directory <- paste(args[1], "tissue/", sep="/")
annotated_probes <- read.csv(args[2], stringsAsFactors=FALSE, header=TRUE) # Read annotation file
annotated_probes <- annotated_probes[!is.na(annotated_probes$ID),] # Only consider probes that don't have NA IDs

# Load R packages
pkgs <- c("limma", "aggregation")
invisible(lapply(pkgs, function(x) suppressWarnings(suppressMessages(library(x, character.only=TRUE)))))

### Function: getMicroarrayData
### Retrieves normalized microarray data for a given sample
getMicroarrayData <- function(input_directory, tissue_type) {
    read.table(paste(input_directory, tissue_type, ".txt", sep=""), sep="\t", stringsAsFactors=FALSE, header=TRUE)
}

### Function: diffExpr
### Performs differential gene expression analysis
diffExpr <- function(data.control, data.myc, batch = NULL) {
    n.control <- ncol(data.control)
    n.myc <- ncol(data.myc)
    data.combined <- cbind(data.control, data.myc)
    timecourse <- factor(c(rep(0, n.control), rep(1, n.myc)))
    columns <- unique(c(paste("t", timecourse, sep="")))
    if (!is.null(batch)) {
        design <- model.matrix(~0+timecourse+batch)
        columns <- c(columns, paste("batch", 1:(length(unique(batch))-1), sep=""))
    } else {
        design <- model.matrix(~0+timecourse)
    }
    colnames(design) <- columns
    contrast.list <- paste(columns[2], columns[1], sep="-")
    contrast <- do.call(makeContrasts, c(contrast.list, list(levels=design)))
    fit <- lmFit(as.matrix(data.combined), design)
    fit2 <- contrasts.fit(fit, contrast)
    fit2 <- eBayes(fit2)
    results <- topTable(fit2,adjust.method="BH",n=100000,p.value=1, coef=1)
    return(results)
}

### Function: geneDE
### Performs gene-level differential gene expression analysis
geneDE <- function(diffexpr, annotated_probes) {
  diffexpr <- merge(annotated_probes,diffexpr, by.x="ProbeID", by.y=0)
  diffexpr_geneExpr <- aggregate(diffexpr[,"logFC",drop=FALSE], by=list(diffexpr$ID), FUN=function(x) { x[which.max(abs(x))] })
  diffexpr_genePVal <- aggregate(diffexpr[,"P.Value",drop=FALSE], by=list(diffexpr$ID), FUN=function(x) { fisher(x) })
  colnames(diffexpr_geneExpr) <- c("ID", "logFC")
  colnames(diffexpr_genePVal) <- c("ID", "P.Value")
  diffexpr <- merge(diffexpr_geneExpr, diffexpr_genePVal, by="ID")
  diffexpr <- merge(annotated_probes[!duplicated(annotated_probes$ID),c("ID", "Symbol")], diffexpr, by="ID")
  diffexpr[,"adj.P.Val"] <- p.adjust(diffexpr[,"P.Value"], method="BH")
  diffexpr <- diffexpr[order(diffexpr[,"adj.P.Val"]),]
  return(diffexpr)
}
### Function: geneVals
### Gets aggregated probe-intensity values at gene-level
geneVals <- function(df, annotated_probes) {
  df_cols <- colnames(df)
  df <- merge(annotated_probes,df, by.x="ProbeID", by.y=0)
  df <- aggregate(df[,df_cols,drop=FALSE], by=list(df$ID), FUN=mean)
  colnames(df)[1] <- "ID"
  df <- merge(annotated_probes[!duplicated(annotated_probes$ID),c("ID", "Symbol")], df, by="ID")
  return(df)
}


#############################
#### GET MICROARRAY DATA ####
#############################

# Retrieve lung tissue microarray data
lung.series <- "lung"
lung_control <- c("lung_normal_rep1", "lung_normal_rep2")
lung_myc <- c("lung_tumor_myc_rep1", "lung_tumor_myc_rep2", "lung_tumor_myc_rep3", "lung_tumor_myc_rep4", "lung_tumor_myc_rep5")
lung_mycras <- c("lung_tumor_myckras_rep1", "lung_tumor_myckras_rep2", "lung_tumor_myckras_rep3", "lung_tumor_myckras_rep4")
lung_ras <- c("lung_tumor_kras_rep1", "lung_tumor_kras_rep2", "lung_tumor_kras_rep3", "lung_tumor_kras_rep4", "lung_tumor_kras_rep5")
data <- getMicroarrayData(input_directory, lung.series)
data.lung_control <- data[,lung_control]
data.lung_myc <- data[,lung_myc]
data.lung_mycras <- data[,lung_mycras]
data.lung_ras <- data[,lung_ras]
# Lung MYC diffexpr
lung.diffexpr <- diffExpr(data.lung_control, data.lung_myc)
write.csv(lung.diffexpr, file=paste(input_directory, "lung_diffexpr.csv", sep=""))
write.csv(data.lung_control, file=paste(input_directory, "lung_control_values.csv", sep=""))
write.csv(data.lung_myc, file=paste(input_directory, "lung_myc_values.csv", sep=""))
write.csv(geneDE(lung.diffexpr, annotated_probes), file=paste(input_directory, "lung_diffexpr_gene.csv", sep=""), row.names=FALSE)
write.csv(geneVals(data.lung_control, annotated_probes), file=paste(input_directory, "lung_control_values_gene.csv", sep=""), row.names=FALSE)
write.csv(geneVals(data.lung_myc, annotated_probes), file=paste(input_directory, "lung_myc_values_gene.csv", sep=""), row.names=FALSE)

# Lung MYC+KRAS diffexpr
lungmycras.diffexpr <- diffExpr(data.lung_control, data.lung_mycras)
write.csv(lungmycras.diffexpr, file=paste(input_directory, "lung_mycras_diffexpr.csv", sep=""))
write.csv(data.lung_mycras, file=paste(input_directory, "lung_mycras_values.csv", sep=""))
write.csv(geneDE(lungmycras.diffexpr, annotated_probes), file=paste(input_directory, "lung_mycras_diffexpr_gene.csv", sep=""), row.names=FALSE)
write.csv(geneVals(data.lung_mycras, annotated_probes), file=paste(input_directory, "lung_mycras_values_gene.csv", sep=""), row.names=FALSE)

# Lung KRAS diffexpr
lungras.diffexpr <- diffExpr(data.lung_control, data.lung_ras)
write.csv(lungras.diffexpr, file=paste(input_directory, "lung_ras_diffexpr.csv", sep=""))
write.csv(data.lung_ras, file=paste(input_directory, "lung_ras_values.csv", sep=""))
write.csv(geneDE(lungras.diffexpr, annotated_probes), file=paste(input_directory, "lung_ras_diffexpr_gene.csv", sep=""), row.names=FALSE)
write.csv(geneVals(data.lung_ras, annotated_probes), file=paste(input_directory, "lung_ras_values_gene.csv", sep=""), row.names=FALSE)

# Retrieve kidney tissue microarray data
kidney.series <- "kidney"
kidney_control.batch1 <- c("kidney_normal_batch1_rep1", "kidney_normal_batch1_rep2")
kidney_myc.batch1 <- c("kidney_preneo_myc_1wk_batch1", "kidney_preneo_myc_2wk_batch1", "kidney_tumor_myc_4wk_batch1_rep1", "kidney_tumor_myc_4wk_batch1_rep2")
kidney_myc.batch1.relevant <- c("kidney_tumor_myc_4wk_batch1_rep1", "kidney_tumor_myc_4wk_batch1_rep2")
kidney_control.batch2 <- c("kidney_normal_batch2")
kidney_myc.batch2 <- c("kidney_preneo_myc_1wk_batch2_rep1", "kidney_preneo_myc_1wk_batch2_rep2", "kidney_preneo_myc_2wk_batch2_rep1", "kidney_preneo_myc_2wk_batch2_rep2", "kidney_tumor_myc_4wk_batch2")
kidney_myc.batch2.relevant <- c("kidney_tumor_myc_4wk_batch2")
kidney_myc_factors <- c("1wk", "2wk", "4wk", "4wk", "1wk", "1wk", "2wk", "2wk", "4wk")
data <- getMicroarrayData(input_directory, kidney.series)
data.kidney_control.batch1 <- data[,kidney_control.batch1]
data.kidney_control.batch2 <- data[,kidney_control.batch2]
data.kidney_myc.batch1 <- data[,kidney_myc.batch1]
data.kidney_myc.batch2 <- data[,kidney_myc.batch2]
data.kidney_control <- cbind(data.kidney_control.batch1, data.kidney_control.batch2)
data.kidney_myc <- cbind(data.kidney_myc.batch1, data.kidney_myc.batch2)
colnames(data.kidney_control) <- c(kidney_control.batch1, kidney_control.batch2)
colnames(data.kidney_myc) <- c(kidney_myc.batch1, kidney_myc.batch2)
kidney.batches <- factor(c(rep("Batch1", length(kidney_control.batch1)), rep("Batch2", length(kidney_control.batch2)), rep("Batch1", length(kidney_myc.batch1)), rep("Batch2", length(kidney_myc.batch2))))
kidney.batches.relevant <- factor(c(rep("Batch1", length(kidney_control.batch1)), rep("Batch2", length(kidney_control.batch2)), rep("Batch1", length(kidney_myc.batch1.relevant)), rep("Batch2", length(kidney_myc.batch2.relevant))))
kidney.myc.output.order <- c("kidney_preneo_myc_1wk_batch1", "kidney_preneo_myc_1wk_batch2_rep1", "kidney_preneo_myc_1wk_batch2_rep2", "kidney_preneo_myc_2wk_batch1", "kidney_preneo_myc_2wk_batch2_rep1", "kidney_preneo_myc_2wk_batch2_rep2", "kidney_tumor_myc_4wk_batch1_rep1", "kidney_tumor_myc_4wk_batch1_rep2", "kidney_tumor_myc_4wk_batch2")
# Kidney diffexpr
kidney.diffexpr <- diffExpr(data.kidney_control, data.kidney_myc[,c(kidney_myc.batch1.relevant,kidney_myc.batch2.relevant)], kidney.batches.relevant)
write.csv(kidney.diffexpr, file=paste(input_directory, "kidney_diffexpr.csv", sep=""))
write.csv(geneDE(kidney.diffexpr, annotated_probes), file=paste(input_directory, "kidney_diffexpr_gene.csv", sep=""), row.names=FALSE)
# batch-correct the data before outputting
condition <- factor(c(rep("Normal", ncol(data.kidney_control)), kidney_myc_factors))
design <- model.matrix(~0+condition)
kidney.data.batchcorrect <- cbind(data.kidney_control, data.kidney_myc)
kidney.data.batchcorrect <- removeBatchEffect(kidney.data.batchcorrect, design=design, batch=kidney.batches)
write.csv(kidney.data.batchcorrect[,1:ncol(data.kidney_control)], file=paste(input_directory, "kidney_control_values.csv", sep=""))
write.csv(kidney.data.batchcorrect[,kidney.myc.output.order], file=paste(input_directory, "kidney_myc_values.csv", sep=""))
write.csv(geneVals(kidney.data.batchcorrect[,1:ncol(data.kidney_control)], annotated_probes), file=paste(input_directory, "kidney_control_values_gene.csv", sep=""), row.names=FALSE)
write.csv(geneVals(kidney.data.batchcorrect[,kidney.myc.output.order], annotated_probes), file=paste(input_directory, "kidney_myc_values_gene.csv", sep=""), row.names=FALSE)


# Retrieve liver tissue microarray data
liver.series <- "liver"
liver_control.batch2 <- c("liver_normal_batch2_rep1", "liver_normal_batch2_rep2")
liver_control.batch3 <- c("liver_normal_batch3_rep1", "liver_normal_batch3_rep2")
liver_myc.batch1 <- c("liver_tumor_myc_batch1_rep1", "liver_tumor_myc_batch1_rep2")
liver_myc.batch2 <- c("liver_tumor_myc_batch2_rep1", "liver_tumor_myc_batch2_rep2")
liver_myc.batch3 <- c("liver_tumor_myc_batch3_rep1", "liver_tumor_myc_batch3_rep2")
data <- getMicroarrayData(input_directory, liver.series)
data.liver_control.batch2 <- data[,liver_control.batch2]
data.liver_control.batch3 <- data[,liver_control.batch3]
data.liver_myc.batch1 <- data[,liver_myc.batch1]
data.liver_myc.batch2 <- data[,liver_myc.batch2]
data.liver_myc.batch3 <- data[,liver_myc.batch3]
data.liver_control <- cbind(data.liver_control.batch2, data.liver_control.batch3)
data.liver_myc <- cbind(data.liver_myc.batch1, data.liver_myc.batch2, data.liver_myc.batch3)
colnames(data.liver_control) <- c(liver_control.batch2, liver_control.batch3)
colnames(data.liver_myc) <- c(liver_myc.batch1, liver_myc.batch2, liver_myc.batch3)
liver.batches <- factor(c(rep("Batch2", length(liver_control.batch2)), rep("Batch3", length(liver_control.batch3)), rep("Batch1", length(liver_myc.batch1)), rep("Batch2", length(liver_myc.batch2)), rep("Batch3", length(liver_myc.batch3))))
# Liver diffexpr
liver.diffexpr <- diffExpr(data.liver_control, data.liver_myc, liver.batches)
write.csv(liver.diffexpr, file=paste(input_directory, "liver_diffexpr.csv", sep=""))
write.csv(geneDE(liver.diffexpr, annotated_probes), file=paste(input_directory, "liver_diffexpr_gene.csv", sep=""), row.names=FALSE)
# batch-correct the data before outputting
condition <- factor(c(rep("Normal", ncol(data.liver_control)), rep("MYC", ncol(data.liver_myc))))
design <- model.matrix(~0+condition)
liver.data.batchcorrect <- cbind(data.liver_control, data.liver_myc)
liver.data.batchcorrect <- removeBatchEffect(liver.data.batchcorrect, design=design, batch=liver.batches)
write.csv(liver.data.batchcorrect[,1:ncol(data.liver_control)], file=paste(input_directory, "liver_control_values.csv", sep=""))
write.csv(liver.data.batchcorrect[,c(liver_myc.batch1, liver_myc.batch2, liver_myc.batch3)], file=paste(input_directory, "liver_myc_values.csv", sep=""))
write.csv(geneVals(liver.data.batchcorrect[,1:ncol(data.liver_control)], annotated_probes), file=paste(input_directory, "liver_control_values_gene.csv", sep=""), row.names=FALSE)
write.csv(geneVals(liver.data.batchcorrect[,c(liver_myc.batch1, liver_myc.batch2, liver_myc.batch3)], annotated_probes), file=paste(input_directory, "liver_myc_values_gene.csv", sep=""), row.names=FALSE)
