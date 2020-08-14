###############################################################
### File: rankaggregation.r
### Description: Pair-wise pearson correlations with MYC and
###              rank aggregation of correlation coefficients.
### Written by Delaney Sullivan
###############################################################

# Load R packages
pkgs <- c("RobustRankAggreg")
invisible(lapply(pkgs, function(x) suppressWarnings(suppressMessages(library(x, character.only=TRUE)))))

# Set up some variables
args <- commandArgs(TRUE) # Extract command-line arguments
tcga_file <- args[1] # TCGA data file (.tsv.gz)
info_file <- args[2] # TCGA annotation file (.csv)
output_prefix <- args[3]

# Iterate through folders containing the data:

data <- read.table(tcga_file, stringsAsFactors=FALSE, check.names=FALSE, sep="\t", header=TRUE)
annotation <- read.csv(info_file, stringsAsFactors=FALSE, check.names=FALSE, header=TRUE, na.strings="",row.names=1)

# Run correlation analysis
data_cor <- NULL
for (study in unique(annotation[,"Study Abbreviation"])) {
  print(paste("Computing pair-wise correlations for:", study))
  study_cols <- annotation[annotation[,"Study Abbreviation"] == study & annotation$Tumor == 1,"Sample"]
  study_data <- data[,study_cols]
  study_data.myc <- unlist(study_data["ENSG00000136997",]) # The MYC row
  cors <- apply(study_data, 1, y=study_data.myc, function(x, y) cor.test(y,x)$estimate ) # This might spit out a bunch of warnings for cases where std dev of x is zero
  cors <- data.frame(cors)
  colnames(cors) <- study
  if (is.null(data_cor)) {
    data_cor <- cors
  } else {
    data_cor <- cbind(data_cor, cors)
  }
}

ranked_lists <- list()
data_rank <- data_cor[rowSums(is.na(data_cor)) <= 0.90*ncol(data_cor),] # Remove rows where 90% are NA
data_rank <- data_rank[!(rownames(data_rank) %in% rownames(data_rank["ENSG00000136997",])),] # Remove MYC
median_cor <- apply(data_rank, 1, function(x) { median(x, na.rm=TRUE) })
for (study in colnames(data_rank)) {
  curr_list <- data_rank[!is.na(data_rank[,study]),study,drop=FALSE]
  ranked_lists[[study]] <- rownames(curr_list[order(-curr_list[,study]),study,drop=FALSE])
}
ranks_aggregated <- aggregateRanks(ranked_lists)
ranks_aggregated$padj <- p.adjust(ranks_aggregated$Score, method="BH")
ranks_aggregated <- merge(ranks_aggregated, data.frame(median_cor), by="row.names", all.x=TRUE, all.y=FALSE)
rownames(ranks_aggregated) <- ranks_aggregated[,1]
ranks_aggregated[,1] <- NULL
ranks_aggregated <- ranks_aggregated[,c("Score", "padj", "median_cor")]
colnames(ranks_aggregated) <- c("pval", "padj", "median_pearson")

# Output rankings and correlations

write.csv(data_cor, file=paste(output_prefix, "myc_correlations.csv", sep=""), quote=FALSE)
write.csv(ranks_aggregated, file=paste(output_prefix, "myc_rra.csv", sep=""), quote=FALSE)

