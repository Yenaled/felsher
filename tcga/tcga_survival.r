# TCGA survival analysis
# Usage: Rscript tcga_survival.r <tcga_survival_xlsx_file> <tcga_gene_expression_tsv_file> <tcga_info_csv_file> <num_threads> <output_file>
# Example: Rscript tcga_survival.r TCGA-CDR-SupplementalTableS1.xlsx TCGA.processed.tumors_corrected.tsv.gz TCGA.processed.info_corrected.csv 4 tcga_zscores.csv

# Load R packages
require(openxlsx)
require(survival)
require(parallel)

# Read in command-line arguments
args <- commandArgs(TRUE)
tcga_survival_file <- args[1]
gene_expression_file <- args[2]
tcga_info_file <- args[3]
num_cores <- as.numeric(args[4])
output_file <- args[5]

# Read in files
survival_data <- read.xlsx(tcga_survival_file)
gene_expression_data <- read.table(gene_expression_file, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
rownames(gene_expression_data) <- gsub('\\..*', '', rownames(gene_expression_data))
colnames(gene_expression_data) <- substr(colnames(gene_expression_data), 0, 12)
tcga_info <- read.csv(tcga_info_file, stringsAsFactors=FALSE, row.names=1, check.names=FALSE)


# Set up some stuff
tcga_studies <- unique(tcga_info[,"Study Name"])
all_genes <- rownames(gene_expression_data)
is.error <- function(x) inherits(x, "try-error")

# Print out data set dimensions
print(paste("Gene expression data contains ", length(all_genes), " genes and ", ncol(gene_expression_data), " patients", sep=""))

zscores_data <- mclapply(tcga_studies, function(study) {
   sampleIDs <- tcga_info[tcga_info[,"Study Name"] == study,"SampleID"]
   sampleIDs <- intersect(sampleIDs, colnames(gene_expression_data))
   sampleIDs <- intersect(sampleIDs, survival_data$bcr_patient_barcode)
   if (length(sampleIDs) < 1) {
       return(NULL)
   }
   gene_expression_data_study <- t(gene_expression_data[,sampleIDs])
   survival_data_study <- survival_data[survival_data$bcr_patient_barcode %in% sampleIDs,]
   print(paste(study, ": ", length(sampleIDs), " patients", sep=""))
   data_study <- merge(gene_expression_data_study, survival_data_study, by.x="row.names", by.y="bcr_patient_barcode", all=FALSE)
   zscores <- lapply(all_genes, function(x) {
      current_gene <- x
      design_string <- paste("Surv(OS.time, OS) ~ ", current_gene, sep="")
      success <- try({ fit <- coxph(formula = eval(parse(text=design_string)), data = data_study) })
      if (!is.error(success)) {
         cox <- summary(fit)
         return(as.numeric((cox$coefficients)[,"z"]))
      } else {
         return(NA)
      }
   })
   print(paste("Successfully completed:", study))
   return(zscores)
}, mc.cores=num_cores)

processZScoreMatrix <- function(zscores_df, output_csv=NULL) {
    zscores_df <- t(zscores_df)
    is.na(zscores_df) <- do.call(cbind,lapply(zscores_df, is.infinite))
    zscores_df <- zscores_df[rowSums(is.na(zscores_df)) != ncol(zscores_df), ]
    metaz <- data.frame(rowSums(zscores_df, na.rm=TRUE) / sqrt(rowSums(!is.na(zscores_df)))) # Stouffer's method (unweighted)
    colnames(metaz) <- "MetaZ"
    zscores_df <- cbind(metaz, zscores_df)
    zscores_df <- zscores_df[order(-zscores_df$MetaZ),]
    if (!is.null(output_csv)) {
        write.csv(zscores_df, file=output_csv, quote=FALSE)
    }
    return(zscores_df)
}

zscores_1 <- NULL
study_names <- c()
for (i in 1:length(zscores_data)) {
       curr_zscores_1 <- data.frame(zscores_data[[i]])[1,]
       if (ncol(curr_zscores_1) != length(all_genes) || rowSums(is.na(curr_zscores_1)) == length(all_genes)) {
          next
       }
       study_names <- c(study_names, tcga_studies[i])
       colnames(curr_zscores_1) <- all_genes
       if (is.null(zscores_1)) {
          zscores_1 <- curr_zscores_1
       } else {
          zscores_1 <- rbind(zscores_1, curr_zscores_1)
      }
}
rownames(zscores_1) <- study_names
zscores <- processZScoreMatrix(zscores_1, output_file)
