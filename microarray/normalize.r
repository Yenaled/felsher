###############################################################
### File: normalize.r
### Description: Produces quantile-normalized log2-transformed
###              microarray data from raw average probe
###              intensity values.
### Usage: Rscript normalize.r <output_directory> <annotation_file>
### Written by Delaney Sullivan
###############################################################

# Load R packages
pkgs <- c("lumi", "GEOquery", "aggregation")
invisible(lapply(pkgs, function(x) suppressWarnings(suppressMessages(library(x, character.only=TRUE)))))

# Get command line arguments to get output directory that files will be created in as well as the annotated probes file:
args <- commandArgs(TRUE)
output_dir <- args[1]
annotated_probes <- read.csv(args[2], stringsAsFactors=FALSE, header=TRUE) # Read annotation file
annotated_probes <- annotated_probes[!is.na(annotated_probes$ID), "ProbeID"] # Only consider probes that don't have NA IDs
output_dir_raw <- paste(output_dir, "/", "raw/", sep="")
output_dir_tissue_geo <- paste(output_dir, "/", "tissue_geo/", sep="")
output_dir_cells_geo <- paste(output_dir, "/", "cells_geo/", sep="")
output_dir_tissue <- paste(output_dir, "/", "tissue/", sep="")
output_dir_cells <- paste(output_dir, "/", "cells/", sep="")
output_dir_filtered <- paste(output_dir, "/", "raw_filtered/", sep="")
dir.create(output_dir_raw, recursive=TRUE)
dir.create(output_dir_tissue, recursive=TRUE)
dir.create(output_dir_cells, recursive=TRUE)
dir.create(output_dir_tissue_geo, recursive=TRUE)
dir.create(output_dir_cells_geo, recursive=TRUE)
dir.create(output_dir_filtered, recursive=TRUE)
pval_threshold <- 0.01

# Obtain GEO records
cell_accession <- "GSE143250"
tissue_accession <- "GSE143253"
geo_accessions <- c(cell_accession, tissue_accession)
filenames <- sapply(geo_accessions, function(x, output_dir_raw) {
  rownames(getGEOSuppFiles(x, fetch_files = TRUE, baseDir=output_dir_raw, filter_regex="*non-normalized.txt.gz"))
}, output_dir_raw=output_dir_raw)

# Unzip and process files
filenames <- sapply(filenames, function(x) {
  unzipped_filename_prefix <- gsub("[.]txt.gz$", "", x)
  gzip_file <- gzfile(x)
  data <- read.table(gzip_file, sep="\t", stringsAsFactors=FALSE, check.names=TRUE, header=TRUE)
  sample_names <- colnames(data)[c(FALSE, TRUE)]
  colnames(data)[c(FALSE, TRUE)] <- paste(sample_names, "AVG_Signal", sep=".")
  colnames(data)[c(TRUE, FALSE)] <- c("PROBE_ID", paste(sample_names, "Detection Pval", sep="."))
  tissue_types <- table(gsub("\\_.*", "", colnames(data[,2:ncol(data)])))
  col_index <- 2
  unzipped_filenames <- c()
  for (tissue_type in names(tissue_types)) {
    unzipped_filename <- paste(unzipped_filename_prefix, "_", tissue_type, ".txt", sep="")
    unzipped_filename <- paste(output_dir_raw, basename(unzipped_filename), sep="")
    unzipped_filenames <- c(unzipped_filenames, unzipped_filename)
    data_to_write <- data[,c(1,col_index:(col_index+tissue_types[tissue_type]-1))]
    if (grepl(cell_accession, basename(unzipped_filename_prefix))) {
      # For cell line data, only consider myc-driven cell lines (for the purposes of the paper)
      data_to_write <- data_to_write[,c(colnames(data_to_write)[1], colnames(data_to_write)[grepl("myc", colnames(data_to_write))])]
    }
    write.table(data_to_write, file=unzipped_filename, sep="\t", quote=FALSE, row.names=FALSE)
    col_index <- col_index + tissue_types[tissue_type]
  }
  return(unzipped_filenames)
}
)
filenames <- unlist(filenames)

# Iterate through filenames and apply lumi to normalize them (for GEO)
# Also get probes that exceed detection p-value threshold for filtering later on
probes_highP_tissue <- NULL
for (f in filenames) {
  x.lumi <- lumiR(f)
  lumi.Tlog2 <- lumiT(x.lumi, method="log2")
  batch_num <- 0
  data_final <- NULL # Normalized probe intensity values
  data_final_p <- NULL # Detection P-values
  while (TRUE) { # Iterate through batches present in data
    batch_num <- batch_num + 1
    batch_name <- paste("batch", batch_num, sep="")
    samples <- sampleNames(lumi.Tlog2)
    if (TRUE %in% grepl(batch_name, samples)) {
      samples <- samples[grepl(batch_name, samples)]
      lumi.Nquantile <- lumiN(lumi.Tlog2[,samples], method = "quantile")
      data <- exprs(lumi.Nquantile)
      data_p <- detection(lumi.Nquantile)
      if (is.null(data_final)) {
        data_final <- data
        data_final_p <- data_p
      } else {
        ordering <- rownames(data)
        data_final <- transform(merge(data_final,data,by=0), row.names=Row.names, Row.names=NULL)
        data_final_p <- transform(merge(data_final_p,data_p,by=0), row.names=Row.names, Row.names=NULL)
        data_final <- data_final[ordering,]
        data_final_p <- data_final_p[ordering,]
      }
    } else {
      break
    }
  }
  if (is.null(data_final)) {
    lumi.Nquantile <- lumiN(lumi.Tlog2, method = "quantile")
    data_final <- exprs(lumi.Nquantile)
    data_final_p <- detection(lumi.Nquantile)
  }
  output_file <- gsub("(.*)non-normalized_","",basename(f))
  if (grepl(cell_accession, basename(f))) { # Write out processed cell line expression data
    write.table(data_final, file=paste(output_dir_cells_geo, output_file, sep=""), sep="\t", quote=FALSE, row.names=TRUE)
    write.table(data_final_p, file=paste(output_dir_cells_geo, "Detection_Pval_", output_file, sep=""), sep="\t", quote=FALSE, row.names=TRUE)
  } else { # Write out processed tissue expression data
    badprobes <- rownames(data_final_p[rowSums(data_final_p > pval_threshold) == ncol(data_final_p), ])
    if (is.null(probes_highP_tissue)) {
      probes_highP_tissue <- badprobes
    } else {
      probes_highP_tissue <- c(probes_highP_tissue[!(probes_highP_tissue %in% rownames(data_final_p))], intersect(probes_highP_tissue, badprobes))
    }
    write.table(data_final, file=paste(output_dir_tissue_geo, output_file, sep=""), sep="\t", quote=FALSE, row.names=TRUE)
    write.table(data_final_p, file=paste(output_dir_tissue_geo, "Detection_Pval_", output_file, sep=""), sep="\t", quote=FALSE, row.names=TRUE)
  }
}

probes_highP_tissue <- unique(unlist(probes_highP_tissue))
print(paste("Annotated probe set size:", length(annotated_probes)))
print(paste("Probes with high p-values:", length(probes_highP_tissue)))

# Iterate through filenames and apply lumi to normalize them (for paper)
# (Unlike the GEO submission, probes are filtered beforehand)
for (f in filenames) {
  filename <- paste(output_dir_filtered, basename(f), sep="")
  data <- read.table(f, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, sep="\t")
  data <- data[data$PROBE_ID %in% annotated_probes,]
  data <- data[!(data$PROBE_ID %in% probes_highP_tissue),]
  write.table(data, file=filename, sep="\t", quote=FALSE, row.names=FALSE)
  
  x.lumi <- lumiR(filename)
  batch_num <- 0
  data_final <- NULL # Normalized probe intensity values
  while (TRUE) { # Iterate through batches present in data
    batch_num <- batch_num + 1
    batch_name <- paste("batch", batch_num, sep="")
    samples <- sampleNames(x.lumi)
    if (TRUE %in% grepl(batch_name, samples)) {
      samples <- samples[grepl(batch_name, samples)]
      lumi.Nquantile <- lumiT(lumiN(x.lumi[,samples], method = "quantile"), method="log2")
      data <- exprs(lumi.Nquantile)
      if (is.null(data_final)) {
        data_final <- data
      } else {
        ordering <- rownames(data)
        data_final <- transform(merge(data_final,data,by=0), row.names=Row.names, Row.names=NULL)
        data_final <- data_final[ordering,]
      }
    } else {
      break
    }
  }
  if (is.null(data_final)) {
    lumi.Nquantile <- lumiT(lumiN(x.lumi[,samples], method = "quantile"), method="log2")
    data_final <- exprs(lumi.Nquantile)
  }
  output_file <- gsub("(.*)non-normalized_","",basename(f))
  if (grepl(cell_accession, basename(f))) { # Write out processed cell line expression data
    write.table(data_final, file=paste(output_dir_cells, output_file, sep=""), sep="\t", quote=FALSE, row.names=TRUE)
  } else { # Write out processed tissue expression data
    write.table(data_final, file=paste(output_dir_tissue, output_file, sep=""), sep="\t", quote=FALSE, row.names=TRUE)
  }
}
