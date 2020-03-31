###############################################################
### File: download.r
### Description: Downloads microarray normalized probe
###              intensity values from GEO.
### Usage: Rscript download.r <output_directory>
### Written by Delaney Sullivan
###############################################################

# Load R packages
pkgs <- c("GEOquery")
invisible(lapply(pkgs, function(x) suppressWarnings(suppressMessages(library(x, character.only=TRUE)))))

# Get command line arguments to get output directory that files will be created in:
args <- commandArgs(TRUE)
output_dir <- args[1]
output_dir_tissue <- paste(output_dir, "/", "tissue/", sep="")
output_dir_cells <- paste(output_dir, "/", "cells/", sep="")
dir.create(output_dir_tissue, recursive=TRUE)
dir.create(output_dir_cells, recursive=TRUE)

# Obtain GEO records
cell_accession <- "GSE143250"
tissue_accession <- "GSE143253"

# Tissue data:
data <- getGEO(tissue_accession)[[1]]
data_mapping <- pData(data)[,c("title","geo_accession")]
data <- exprs(data)
sample_type <- "lung"
data_to_write <- data[,data_mapping[grepl(sample_type, data_mapping$title),"geo_accession"]]
colnames(data_to_write) <- data_mapping$title[match(colnames(data_to_write), data_mapping$geo_accession)]
write.table(data_to_write, file=paste(output_dir_tissue, sample_type, ".txt", sep=""), sep="\t", quote=FALSE, row.names=TRUE)
sample_type <- "kidney"
data_to_write <- data[,data_mapping[grepl(sample_type, data_mapping$title),"geo_accession"]]
colnames(data_to_write) <- data_mapping$title[match(colnames(data_to_write), data_mapping$geo_accession)]
write.table(data_to_write, file=paste(output_dir_tissue, sample_type, ".txt", sep=""), sep="\t", quote=FALSE, row.names=TRUE)
sample_type <- "liver"
data_to_write <- data[,data_mapping[grepl(sample_type, data_mapping$title),"geo_accession"]]
colnames(data_to_write) <- data_mapping$title[match(colnames(data_to_write), data_mapping$geo_accession)]
write.table(data_to_write, file=paste(output_dir_tissue, sample_type, ".txt", sep=""), sep="\t", quote=FALSE, row.names=TRUE)

# Cell-line data (myc-only):
data <- getGEO(cell_accession)[[1]]
data_mapping <- pData(data)[,c("title","geo_accession")]
data <- exprs(data)
sample_type <- "osteosarcoma"
data_to_write <- data[,data_mapping[grepl(sample_type, data_mapping$title) & grepl("myc", data_mapping$title),"geo_accession"]]
colnames(data_to_write) <- data_mapping$title[match(colnames(data_to_write), data_mapping$geo_accession)]
write.table(data_to_write, file=paste(output_dir_cells, sample_type, ".txt", sep=""), sep="\t", quote=FALSE, row.names=TRUE)
sample_type <- "kidney"
data_to_write <- data[,data_mapping[grepl(sample_type, data_mapping$title) & grepl("myc", data_mapping$title),"geo_accession"]]
colnames(data_to_write) <- data_mapping$title[match(colnames(data_to_write), data_mapping$geo_accession)]
write.table(data_to_write, file=paste(output_dir_cells, sample_type, ".txt", sep=""), sep="\t", quote=FALSE, row.names=TRUE)
sample_type <- "liver"
data_to_write <- data[,data_mapping[grepl(sample_type, data_mapping$title) & grepl("myc", data_mapping$title),"geo_accession"]]
colnames(data_to_write) <- data_mapping$title[match(colnames(data_to_write), data_mapping$geo_accession)]
write.table(data_to_write, file=paste(output_dir_cells, sample_type, ".txt", sep=""), sep="\t", quote=FALSE, row.names=TRUE)
sample_type <- "lymphoma"
data_to_write <- data[,data_mapping[grepl(sample_type, data_mapping$title) & grepl("myc", data_mapping$title),"geo_accession"]]
colnames(data_to_write) <- data_mapping$title[match(colnames(data_to_write), data_mapping$geo_accession)]
write.table(data_to_write, file=paste(output_dir_cells, sample_type, ".txt", sep=""), sep="\t", quote=FALSE, row.names=TRUE)
