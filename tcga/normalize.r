###############################################################
### File: normalize.r
### Description: Gets normalized gene expression values
###              from processed GDC gene data.
### Written by Delaney Sullivan
###############################################################

# Load R packages
pkgs <- c("DESeq2", "BiocParallel", "sva")
invisible(lapply(pkgs, function(x) suppressWarnings(suppressMessages(library(x, character.only=TRUE)))))

# Set up some variables
args <- commandArgs(TRUE) # Extract command-line arguments
project.name <- args[1]
current_directory <- args[2]
project.path <- paste(current_directory, "/processed/",project.name,"/",sep="")
output.dir <- paste(current_directory, "/normalized/",project.name,"/",sep="")
num.cores <- 4
doParallel <- TRUE
counts_prefix <- "counts"

# Set up
if (doParallel) { # Register multiprocessing
    register(MulticoreParam(num.cores))
}
dir.create(output.dir, showWarnings=FALSE, recursive=TRUE) # Create output directory if it doesn't already exist

# Read in counts data (rows=genes; columns=subjects) and prepare counts matrix
cts <- read.csv(paste(project.path, paste(counts_prefix, "_experimental.csv",sep=""),sep=""), check.names=FALSE)
rownames(cts) <- gsub("\\..*", "", cts[,1]) # Format ENSEMBL IDs by removing periods and everything after periods
cts[,1] <- NULL
cts <- cts[!startsWith(rownames(cts), "__"), ] # Get rid of some metadata present in the counts datafile

# Read in batch info
batch.data <- read.table(paste(project.path, "batch.txt",sep=""), header=TRUE)
rownames(batch.data) <- batch.data[,1]
batch.data[,1] <- NULL

# Prepare coldata (rows=subjects; columns=batch and amplification) for DESeq2
coldata <- data.frame(row.names=colnames(cts))
coldata <- transform(merge(coldata, batch.data, by=0, all.x=TRUE), row.names=Row.names, Row.names=NULL)
coldata$Batch <- factor(coldata$Batch)

# DESeq2 analysis to retrieve normalized, batch-corrected counts for correlation analysis
# Note that zero variance genes are removed for the purposes of combat
var.data <- apply(cts, 1, var)
dds <- DESeqDataSetFromMatrix(countData = cts[-which(var.data == 0),], colData = coldata, design= ~ Batch)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
assay(vsd) <- ComBat(assay(vsd), vsd$Batch)
write.csv(as.data.frame(assay(vsd)), file=paste(output.dir, "normalized_", counts_prefix, ".csv", sep=""))

# DESeq2 analysis to retrieve normalized, batch-corrected counts for correlation analysis
# Note that zero variance genes are removed for the purposes of combat
# This also looks at matched controls
ctrl_file <- paste(project.path, paste(counts_prefix, "_control.csv",sep=""), sep="")
if (!file.exists(ctrl_file)) {
  quit()
}
cts1 <- read.csv(ctrl_file, check.names=FALSE)
rownames(cts1) <- gsub("\\..*", "", cts1[,1]) # Format ENSEMBL IDs by removing periods and everything after periods
cts1[,1] <- NULL
cts1 <- cts1[!startsWith(rownames(cts1), "__"), ] # Get rid of some metadata present in the counts datafile
patient_ids <- intersect(colnames(cts), colnames(cts1))
cts <- cts[,patient_ids]
cts1 <- cts1[,patient_ids]
coldata <- data.frame(row.names=colnames(cts))
coldata <- transform(merge(coldata, batch.data, by=0, all.x=TRUE), row.names=Row.names, Row.names=NULL)
coldata <- coldata[colnames(cts),,drop=FALSE]
coldata$ID <- rownames(coldata)
coldata$Experimental <- "experimental"
coldata_control <- coldata
coldata_control$Experimental <- "control"
coldata <- rbind(coldata, coldata_control)
coldata$ID <- factor(coldata$ID)
coldata$Experimental <- factor(coldata$Experimental)
coldata$Batch <- factor(coldata$Batch)
cts <- transform(merge(cts, cts1, by=0, all=FALSE), row.names=Row.names, Row.names=NULL)
colnames(cts) <- rownames(coldata)
var.data <- apply(cts, 1, var)
cts <- cts[-which(var.data == 0),]
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design= ~ Experimental + Batch)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
assay(vsd) <- ComBat(assay(vsd), vsd$Batch)
df <- as.data.frame(assay(vsd))
num_patients <- ncol(df) / 2
df1 <- df[,1:num_patients]
df2 <- df[,(num_patients+1):(num_patients*2)]
colnames(df2) <- colnames(df1)
write.csv(df1, file=paste(output.dir, "experimental_normalized_", counts_prefix, ".csv", sep=""))
write.csv(df2, file=paste(output.dir, "control_normalized_", counts_prefix, ".csv", sep=""))
df_lfc <- df1 - df2
write.csv(df_lfc, file=paste(output.dir, "lfc_normalized_", counts_prefix, ".csv", sep=""))
