###############################################################
### File: processTCGAbatches.r
### Description: Selects relevant TCGA samples and 
###              batch-corrects them.
### Written by Delaney Sullivan
###############################################################

# Load R packages
pkgs <- c("limma")
invisible(lapply(pkgs, function(x) suppressWarnings(suppressMessages(library(x, character.only=TRUE)))))

# Set up some variables
args <- commandArgs(TRUE) # Extract command-line arguments
input_file <- args[1] # "TCGA.rsem.deseq2.log2.tsv.gz"
bcr_file <- args[2] # Maps barcode accession to TCGA study name "bcrBatchCode.tsv"
tss_file <-args[3] # "tissueSourceSite.tsv"
output_prefix <- args[4] # "TCGA.processed."
batch_correct <- FALSE
if (length(args) > 4) {
  batch_file <- args[5] # "tcga_batches.csv"
  batch_correct <- TRUE
}

# Which samples to keep:
# 01: Primary Solid Tumor
# 03: Primary Blood Derived Cancer - Peripheral Blood
# 11: Solid Tissue Normal
sample_types_to_keep_tumor <- c("01", "03")
sample_types_to_keep_normal <- c("11")

# Iterate through folders containing the data:

data <- read.table(input_file, stringsAsFactors=FALSE, check.names=FALSE, sep="\t", header=TRUE)
bcr_map <- read.table(bcr_file, stringsAsFactors=FALSE, check.names=FALSE, sep="\t", header=TRUE, na.strings="", quote="")
tss_map <- read.table(tss_file, stringsAsFactors=FALSE, check.names=FALSE, sep="\t", header=TRUE, na.strings="", quote="")
if (batch_correct) {
  batches <- read.csv(batch_file, stringsAsFactors=FALSE, header=TRUE, na.strings="",row.names=1)
}

# Map tissue source site (TSS) IDs to TCGA study abbrebiations
bcr_map[,"Study Name"] <- trimws(bcr_map[,"Study Name"])
tss_map[,"Study Name"] <- trimws(tss_map[,"Study Name"])
tss_map <- merge(tss_map[,c("TSS Code", "Study Name")], bcr_map[,c("Study Name", "Study Abbreviation")], all=FALSE, by="Study Name")
tss_map <- tss_map[!duplicated(tss_map),]

# Tumor samples
tumor_samples <- c()
for (i in paste(sample_types_to_keep_tumor, "$", sep="")) {
  tumor_samples <- c(tumor_samples, colnames(data)[grepl(i, colnames(data))])
}
# Normal samples
normal_samples <- c()
for (i in paste(sample_types_to_keep_normal, "$", sep="")) {
  normal_samples <- c(normal_samples, colnames(data)[grepl(i, colnames(data))])
}

# Create data frame to annotate samples by tumor/normal and by study abbreviation and by batches
annotation <- data.frame(row.names=c(tumor_samples,normal_samples))
annotation$Sample <- rownames(annotation)
annotation[rownames(annotation) %in% normal_samples,"Tumor"] <- 0
annotation[rownames(annotation) %in% tumor_samples,"Tumor"] <- 1
annotation[,"TSS Code"] <- substr(rownames(annotation), 6, 7)
annotation <- merge(annotation, tss_map, by="TSS Code", all.x=TRUE, all.y=FALSE)
annotation[,"SampleID"] <- gsub("(.*)\\-.*", "\\1", annotation$Sample)
if (batch_correct) {
  annotation <- merge(annotation, batches[,"batch.number",drop=FALSE], by.x="SampleID", by.y="row.names", all.x=TRUE, all.y=FALSE)
}

# Filter data to only contain the samples specified in the annotation
data <- data[,annotation$Sample]

# Run tumor minus normal analysis
data_lfc <- NULL
for (study in unique(annotation[,"Study Name"])) {
  study_samples <- annotation[annotation[,"Study Name"] == study & annotation$Tumor == 0,"SampleID"]
  study_samples <- intersect(study_samples, annotation[annotation[,"Study Name"] == study & annotation$Tumor == 1,"SampleID"])
  if (length(study_samples) > 0) {
    print(paste("Computing fold change tumor/normal for:", study))
    normal_sampleNames <- annotation[annotation[,"SampleID"] %in% study_samples & annotation$Tumor == 0,"Sample"]
    tumor_sampleNames <- annotation[annotation[,"SampleID"] %in% study_samples & annotation$Tumor == 1,"Sample"]
    study_lfc <- data[,tumor_sampleNames] - data[,normal_sampleNames]
    colnames(study_lfc) <- study_samples
    if (is.null(data_lfc)) {
      data_lfc <- study_lfc
    } else {
      data_lfc <- cbind(data_lfc, study_lfc)
    }
  }
}
gz_file <- gzfile(paste(output_prefix, "tumor_normal_lfc.tsv.gz", sep=""), "w")
write.table(data_lfc, file=gz_file, sep="\t", quote=FALSE)
close(gz_file)

# Output tumor expression data without batch correction
gz_file <- gzfile(paste(output_prefix, "tumors.tsv.gz", sep=""), "w")
write.table(data[,annotation[annotation$Tumor == 1,"Sample"]], file=gz_file, sep="\t", quote=FALSE)
close(gz_file)

# Output annotation
write.csv(annotation, file=paste(output_prefix, "info.csv", sep=""), quote=FALSE)

if (!batch_correct) {
  quit()
}

# Remove where there's an unidentified batch 
unidentified_batches <- annotation[is.na(annotation[,"batch.number"]),"Sample"]
if (length(unidentified_batches) > 0) {
  print(paste("Removing the following unidentified batches:", paste(unidentified_batches, collapse=",")))
  annotation <- annotation[!is.na(annotation[,"batch.number"]),]
}

# Filter data to only contain the samples specified in the annotation
data <- data[,annotation$Sample]

# Do batch correction
results_tumor <- NULL
pdf(paste(output_prefix, "batch_effects_view.pdf", sep=""))
for (study in unique(annotation[,"Study Name"])) {
  # Batch correct tumor samples for current study
  study_cols <- colnames(data)
  study_cols <- study_cols[study_cols %in% annotation[annotation[,"Study Name"] == study & annotation$Tumor == 1,"Sample"]]
  if (length(study_cols) > 0) {
    study_data <- data[,study_cols]
    study_data_corrected <- study_data
    study_batches <- annotation[annotation$Sample %in% study_cols,c("SampleID","TSS Code","batch.number")]
    if (length(unique(study_batches[,"batch.number"])) > 1 || length(unique(study_batches[,"TSS Code"])) > 1) {
      print(paste("Correcting:", study))
      study_batches1 <- study_batches[,"batch.number"]
      study_batches2 <- study_batches[,"TSS Code"]
      study_data_corrected <- removeBatchEffect(study_data, batch=study_batches1, batch2=study_batches2)
      
      # Colors for plotting PCA
      colors_fixed <- c("black", "purple", "orange", "gray", "yellow", "brown", "green")
      colorPalette <- c("blue", "white", "red")
      # Plot uncorrected data
      study_data_t <- t(study_data[apply(study_data, 1, var) > 0.1, ]) # Transpose and remove near-zero variance features
      pca <- prcomp(study_data_t, center = TRUE, scale. = TRUE)
      labels <- as.numeric(factor(study_batches1))
      colors <- c(colors_fixed, colorRampPalette(colorPalette)(max(0,length(unique(labels)) - length(colors_fixed))))
      plot(pca$x[,1:2], bg=colors[labels], pch=21, cex.lab=1.5, cex.axis=2, cex=2, main=paste(study, "Uncorrected - Batch Number"))
      labels <- as.numeric(factor(study_batches2))
      colors <- c(colors_fixed, colorRampPalette(colorPalette)(max(0,length(unique(labels)) - length(colors_fixed))))
      plot(pca$x[,1:2], bg=colors[labels], pch=21, cex.lab=1.5, cex.axis=2, cex=2, main=paste(study, "Uncorrected - Tissue Source Site"))
      # Plot corrected data
      study_data_corrected_t <- t(study_data_corrected[apply(study_data_corrected, 1, var) > 0.1, ]) # Transpose and remove near-zero variance features
      pca <- prcomp(study_data_corrected_t, center = TRUE, scale. = TRUE)
      labels <- as.numeric(factor(study_batches1))
      colors <- c(colors_fixed, colorRampPalette(colorPalette)(max(0,length(unique(labels)) - length(colors_fixed))))
      plot(pca$x[,1:2], bg=colors[labels], pch=21, cex.lab=1.5, cex.axis=2, cex=2, main=paste(study, "Corrected - Batch Number"))
      labels <- as.numeric(factor(study_batches2))
      colors <- c(colors_fixed, colorRampPalette(colorPalette)(max(0,length(unique(labels)) - length(colors_fixed))))
      plot(pca$x[,1:2], bg=colors[labels], pch=21, cex.lab=1.5, cex.axis=2, cex=2, main=paste(study, "Corrected - Tissue Source Site"))
    }
  }
  if (is.null(results_tumor)) {
    results_tumor <- study_data_corrected
  } else {
    results_tumor <- cbind(results_tumor, study_data_corrected)
  }
  }
device <- dev.off()

results_tumor <- round(results_tumor, 4)
gz_file <- gzfile(paste(output_prefix, "tumors_corrected.tsv.gz", sep=""), "w")
write.table(results_tumor, file=gz_file, sep="\t", quote=FALSE)
close(gz_file)
write.csv(annotation, file=paste(output_prefix, "info_corrected.csv", sep=""), quote=FALSE)
