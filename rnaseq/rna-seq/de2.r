# Differential gene expression analysis using sleuth

# Load R packages
pkgs <- c("GenomicFeatures", "sleuth")
invisible(lapply(pkgs, function(x) suppressWarnings(suppressMessages(library(x, character.only=TRUE)))))

# Process command-line arguments
args <- commandArgs(TRUE)
num_cores <- as.numeric(args[1])
input.path <- args[2] # Path to input files
prefix <- args[3] # output file prefix
samples <- unlist(strsplit(args[4],",")) # Sample names separated by commas
treatment <- unlist(strsplit(args[5],",")) # Which samples are treatment or control samples (e.g. c,c,c,t,t,t)
gtf_file <- args[6]

# Prepare metadata data frame
metadata <- data.frame(sample=samples, treatment=treatment, path=paste(input.path, "/", samples, "/", "abundance.h5", sep=""), stringsAsFactors=FALSE)

# Ensure that all abundance.h5 files have the same transcripts and obtain the transcript identifiers
identifiers <- NULL
for (file in metadata$path) {
    abundances <- read_kallisto_h5(file, read_bootstrap = FALSE)$abundance
    if (is.null(identifiers)) {
        identifiers <- abundances$target_id
    } else {
        if (!(setequal(identifiers, abundances$target_id) && length(identifiers) == length(abundances$target_id))) {
            print("Error: Transcript identifiers in supplied files are not identical")
            quit()
        }
    }
}

length_identifiers <- length(identifiers)
print(paste("Found", length(identifiers), "transcripts"))
identifiers <- data.frame(target_id=identifiers, transcript_id=gsub("\\|.*","",identifiers), stringsAsFactors=FALSE)

# Map transcript IDs to genes for gene-level aggregation, based on supplied GTF file
txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")
mapping <- read.table(gtf_file, header=FALSE, stringsAsFactors=FALSE, sep="\t")
mapping <- mapping[,9]
mapping <- lapply(mapping, function(x) {
    a <- trimws(unlist(strsplit(x, ";")))
    transcript_id = substring(a[grepl("transcript_id", a)], 15)
    if (length(transcript_id) == 0) {
        return(NULL)
    }
    return(c("gene_id"=substring(a[grepl("gene_id", a)], 9), "transcript_id"=transcript_id, "gene_name"=gsub("\"", "", substring(a[grepl("gene_name", a)], 11))))
})
mapping <- as.data.frame(do.call(rbind, mapping))
mapping <- mapping[!duplicated(mapping),]
mapping <- merge(mapping, identifiers, by="transcript_id", all.x=FALSE, all.y=FALSE)
if (nrow(mapping) != length_identifiers) {
    print("Error: GTF file transcript identifier does not coincide with identifiers present in kallisto abundances files")
    quit()
}

# Perform differential gene expression analysis

## Transcript-level analysis ##
so <- sleuth_prep(metadata, target_mapping = mapping, aggregation_column = 'gene_name', num_cores=num_cores, transform_fun_counts = function(x) log2(x + 0.5))
so <- sleuth_fit(so, ~treatment, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_wt(so, 'treatmentt', 'full')
# P-value aggregation
sleuth_table <- sleuth_results(so, 'treatmentt', test_type='wt', pval_aggregate = TRUE, show_all = TRUE)
write.csv(sleuth_table, paste(input.path, "/", prefix, "_sleuth_results_pval_aggregated.csv", sep=""), quote=FALSE)
# Transcript differential expression results
sleuth_table <- sleuth_results(so, 'treatmentt', test_type='wt', pval_aggregate = FALSE, show_all = TRUE)
write.csv(sleuth_table, paste(input.path, "/", prefix, "_sleuth_results_transcripts.csv", sep=""), quote=FALSE)
id_mapping <- sleuth_table[,c("transcript_id","gene_id","gene_name","target_id")]
# TPM and counts
sleuth_table <- sleuth_to_matrix(so, "obs_raw", "tpm")
sleuth_table <- merge(id_mapping, sleuth_table, by.x="target_id", by.y=0)
write.csv(sleuth_table, paste(input.path, "/", prefix, "_sleuth_tpm_transcripts.csv", sep=""), quote=FALSE)
sleuth_table <- sleuth_to_matrix(so, "obs_raw", "est_counts")
sleuth_table <- merge(id_mapping, sleuth_table, by.x="target_id", by.y=0)
write.csv(sleuth_table, paste(input.path, "/", prefix, "_sleuth_counts_transcripts.csv", sep=""), quote=FALSE)
sleuth_table <- sleuth_to_matrix(so, "obs_norm", "tpm")
sleuth_table <- merge(id_mapping, sleuth_table, by.x="target_id", by.y=0)
write.csv(sleuth_table, paste(input.path, "/", prefix, "_sleuth_normalized_tpm_transcripts.csv", sep=""), quote=FALSE)
sleuth_table <- sleuth_to_matrix(so, "obs_norm", "est_counts")
sleuth_table <- merge(id_mapping, sleuth_table, by.x="target_id", by.y=0)
write.csv(sleuth_table, paste(input.path, "/", prefix, "_sleuth_normalized_counts_transcripts.csv", sep=""), quote=FALSE)

## Gene-level analysis ##
so_gene <- sleuth_prep(metadata, gene_mode=TRUE, target_mapping = mapping, aggregation_column = 'gene_name', num_cores=num_cores, transform_fun_counts = function(x) log2(x + 0.5))
so_gene <- sleuth_fit(so_gene, ~treatment, 'full')
so_gene <- sleuth_fit(so_gene, ~1, 'reduced')
so_gene <- sleuth_wt(so_gene, 'treatmentt', 'full')
sleuth_table_gene <- sleuth_results(so_gene, 'treatmentt', test_type='wt', show_all = TRUE, gene_mode=TRUE)
names(sleuth_table_gene)[names(sleuth_table_gene) == 'target_id'] <- 'gene_name'
sleuth_table_gene$transcript_id <- NULL
sleuth_table_gene <- sleuth_table_gene[!duplicated(sleuth_table_gene$gene_name),]
write.csv(sleuth_table_gene, paste(input.path, "/", prefix, "_sleuth_results_genes.csv", sep=""), quote=FALSE)
id_mapping <- sleuth_table_gene[,c("gene_id","gene_name")]
# TPM and counts
sleuth_table <- sleuth_to_matrix(so_gene, "obs_norm", "scaled_reads_per_base")
sleuth_table <- merge(id_mapping, sleuth_table, by.x="gene_name", by.y=0)
write.csv(sleuth_table, paste(input.path, "/", prefix, "_sleuth_normalized_counts_genes.csv", sep=""), quote=FALSE)
sleuth_table <- sleuth_to_matrix(so_gene, "obs_norm", "scaled_reads_per_base")
sleuth_table <- log2(sleuth_table+1)
sleuth_table <- merge(id_mapping, sleuth_table, by.x="gene_name", by.y=0)
write.csv(sleuth_table, paste(input.path, "/", prefix, "_sleuth_log_normalized_counts_genes.csv", sep=""), quote=FALSE)
sleuth_table <- sleuth_to_matrix(so_gene, "obs_norm", "tpm")
sleuth_table <- merge(id_mapping, sleuth_table, by.x="gene_name", by.y=0)
write.csv(sleuth_table, paste(input.path, "/", prefix, "_sleuth_normalized_tpm_genes.csv", sep=""), quote=FALSE)
