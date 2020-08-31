# Organize downloaded SE BED files from dbSUPER
# Usage: Rscript process_bed_zip.r data/all_mm9_bed.zip data/dbSUPER_mm9.bed data/superenhancer_gene_annotations.csv

args <- commandArgs(TRUE)
input_bed_zip <- args[1]
output_bed <- args[2]
gene_annotations <- read.csv(args[3], stringsAsFactors=FALSE) # Mapping between SE ID and gene symbol
gene_annotations <- gene_annotations[grep(" ", gene_annotations$Gene, invert=TRUE),] # Remove all records containing a space (indicative of multimapping or unmapped)
unzipped_files <- unzip(input_bed_zip)

data <- NULL
for (f in unzipped_files) {
    tissue_name <- sub('\\.bed$', '', basename(f))
    tissue_name <- gsub(" ", "_", tissue_name) # Replace spaces with underscores
    bed_data <- read.table(f, sep="\t", stringsAsFactors=FALSE, header=FALSE)
    bed_data$V5 <- tissue_name
    bed_data <- merge(bed_data, gene_annotations, by.x="V4", by.y="Name", all.x=TRUE, all.y=FALSE)
    bed_data <- bed_data[,c("V1","V2","V3","Gene","V5")]
    bed_data <- bed_data[complete.cases(bed_data),]
    if (is.null(data)) {
        data <- bed_data
    } else {
        data <- rbind(data, bed_data)
    }
}

write.table(data, file=output_bed, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
