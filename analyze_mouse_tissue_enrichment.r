###############################################################
### File: analyze_mouse_tissue_enrichment.r
### Description: Performs mouse tissue enrichment analysis
### Usage: Rscript analyze_mouse_tissue_enrichment.r
### Written by Delaney Sullivan
###############################################################

source("analysis_setup.r")
output_dir <- "./output/enrichr_mouse_tissue/"

# Do Mouse Tissue Enrichment analysis

bonferonni_tests <- length(sample_mapping)*2*length(mouse_tissue_geneset_names)
alpha_cutoff <- 0.05
db_name <- "Mouse_Gene_Atlas"
output_excel_file <- paste(output_dir, "tissue_enrichment_table.xlsx", sep="")
enrichment_data <- do_enrichment(db_name, sample_mapping, genes_de_up, genes_de_down, output_excel_file, mouse_tissue_geneset_names)
enrichment_data_up <- enrichment_data[["up"]]
enrichment_data_dn <- enrichment_data[["down"]]

# Output p-values data
printSigPVals <- function(df, threshold){
   ind <- which(!is.na(df) & df<threshold, arr.ind=TRUE)
   print(paste(colnames(df)[ind[,"col"]], rownames(df)[ind[,"row"]], sep=', '))
}
enrichment_data_up_pval <- enrichment_data_up[,grepl("_P.value$", colnames(enrichment_data_up))]
enrichment_data_dn_pval <- enrichment_data_dn[,grepl("_P.value$", colnames(enrichment_data_dn))]
colnames(enrichment_data_up_pval) <- sample_mapping
colnames(enrichment_data_dn_pval) <- sample_mapping
write.csv(enrichment_data_up_pval, file=paste(output_dir, "enrichr_mouse_gene_atlas_pvals_up.csv", sep=""))
write.csv(enrichment_data_dn_pval, file=paste(output_dir, "enrichr_mouse_gene_atlas_pvals_down.csv", sep=""))
sink(file=paste(output_dir, "log_enrichr_mouse_gene_atlas.txt", sep=""), split=TRUE)
print(paste("UP genes: Significant gene sets at p-value < ", alpha_cutoff / bonferonni_tests, " with Bonferroni correction at alpha ", alpha_cutoff, " for ", bonferonni_tests, " tests", sep=""))
printSigPVals(enrichment_data_up_pval, alpha_cutoff / bonferonni_tests)
print(paste("DOWN genes: Significant gene sets at p-value < ", alpha_cutoff / bonferonni_tests, " with Bonferroni correction at alpha ", alpha_cutoff, " for ", bonferonni_tests, " tests", sep=""))
printSigPVals(enrichment_data_dn_pval, alpha_cutoff / bonferonni_tests)
sink()

# Output odds ratio data
enrichment_data_up <- enrichment_data_up[,grepl("_Odds.Ratio$", colnames(enrichment_data_up))]
enrichment_data_dn <- enrichment_data_dn[,grepl("_Odds.Ratio$", colnames(enrichment_data_dn))]
colnames(enrichment_data_up) <- sample_mapping
colnames(enrichment_data_dn) <- sample_mapping
enrichment_data_up[is.na(enrichment_data_up)] <- 0
enrichment_data_dn[is.na(enrichment_data_dn)] <- 0
write.csv(enrichment_data_up, file=paste(output_dir, "enrichr_mouse_gene_atlas_OR_up.csv", sep=""))
write.csv(enrichment_data_dn, file=paste(output_dir, "enrichr_mouse_gene_atlas_OR_down.csv", sep=""))
