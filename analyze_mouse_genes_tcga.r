###############################################################
### File: analyze_mouse_genes_tcga.r
### Description: Analyze mouse DEGs in human TCGA datasets
### Usage: Rscript analyze_mouse_genes_tcga.r
### Written by Delaney Sullivan
###############################################################

source("analysis_setup.r")
output_dir <- "./output/mouse_tcga/"

tcga_directory <- "./tcga/normalized/"
getTcgaFC <- function(tcga_directory, tcga_study, ensembl_mapping, ortholog_mapping) {
    fname <- paste(tcga_directory, tcga_study, "/lfc_normalized_counts.csv.gz", sep="")
    if (!file.exists(fname)) {
      fname <- paste(tcga_directory, tcga_study, "/lfc_normalized_counts.csv", sep="")
    } else {
      fname <- gzfile(fname)
    }
    tcga_data <- read.csv(fname, stringsAsFactors=FALSE, row.names=1)
    tcga_data <- rowMeans(tcga_data)
    tcga_data <- merge(ensembl_mapping, tcga_data, by.x="Ensembl", by.y=0, all.x=FALSE, all.y=TRUE)
    tcga_data <- merge(ortholog_mapping, tcga_data, by.x="Human", by.y="ID", all.x=FALSE, all.y=TRUE)
    tcga_data <- tcga_data[!is.na(tcga_data$Mouse),]
    tcga_data$Ensembl <- NULL
    tcga_data$Human <- NULL
    return(tcga_data)
}

outputTcgaData <- function(mouse_tissue, tcga_studies, up_ids, down_ids, tcga_directory, ensembl_mapping, ortholog_mapping) {
    data_final_up <- NULL
    data_final_down <- NULL
    for (study in tcga_studies) {
        data <- getTcgaFC(tcga_directory, study, ensembl_mapping, ortholog_mapping)
        data_up <- data[data$Mouse %in% up_ids[,mouse_tissue],]
        data_down <- data[data$Mouse %in% down_ids[,mouse_tissue],]
        if (is.null(data_final_up)) {
            data_final_up <- data_up
            data_final_down <- data_down
        } else {
            data_final_up <- merge(data_final_up, data_up, by="Mouse")
            data_final_down <- merge(data_final_down, data_down, by="Mouse")
        }
    }
    colnames(data_final_up) <- c("ID", toupper(paste("tcga", tcga_studies, sep="-")))
    colnames(data_final_down) <- c("ID", toupper(paste("tcga", tcga_studies, sep="-")))
    return(list(up=data_final_up,down=data_final_down))
}

tcga_studies <- c("lihc", "kirc", "luad")
mouse_studies <- c("liver_myc", "kidney_myc", "lung_myc")
wb <- createWorkbook("Mouse TCGA")
for (i in mouse_studies) {
  data <- outputTcgaData(i, tcga_studies, genes_de_up_ids, genes_de_down_ids, tcga_directory, ensembl_human_mapping, orthologs_ids)
  ws_up_name <- paste(sample_mapping_myc[i], " - UP", sep="")
  ws_down_name <- paste(sample_mapping_myc[i], " - DOWN", sep="")
  addWorksheet(wb, ws_up_name)
  writeData(wb, sheet = ws_up_name, data[["up"]], rowNames = FALSE)
  addWorksheet(wb, ws_down_name)
  writeData(wb, sheet = ws_down_name, data[["down"]], rowNames = FALSE)
  write.csv(data[["up"]], paste(output_dir, i, "_up.csv", sep=""), row.names=FALSE)
  write.csv(data[["down"]], paste(output_dir, i, "_down.csv", sep=""), row.names=FALSE)
}
saveWorkbook(wb, paste(output_dir, "mouse_tcga.xlsx", sep=""), overwrite = TRUE)
