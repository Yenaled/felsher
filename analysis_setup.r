###############################################################
### File: analysis_setup.r
### Description: Contains setup information + functions for 
###              downstream analyses of mouse data.
### Written by Delaney Sullivan
###############################################################

# Load R packages
pkgs <- c("fgsea", "openxlsx", "enrichR", "ComplexHeatmap", "survival")
invisible(lapply(pkgs, function(x) suppressWarnings(suppressMessages(library(x, character.only=TRUE)))))

# Set random seed
set.seed(42)

# Read in some files
raw_genesets_mouse_tissue <- readLines("data/Mouse_Gene_Atlas.txt")
genes_de_up <- read.table("./output/mouse_de/de_genes_up_symbols.txt", stringsAsFactors=FALSE, header=TRUE, sep="\t")
genes_de_down <- read.table("./output/mouse_de/de_genes_down_symbols.txt", stringsAsFactors=FALSE, header=TRUE, sep="\t")
genes_de_up_ids <- read.table("./output/mouse_de/de_genes_up_ids.txt", stringsAsFactors=FALSE, header=TRUE, sep="\t")
genes_de_down_ids <- read.table("./output/mouse_de/de_genes_down_ids.txt", stringsAsFactors=FALSE, header=TRUE, sep="\t")
orthologs_ids <- read.csv("./rnaseq/annotation/human_mouse_orthologs_ids.csv", stringsAsFactors=FALSE)
ensembl_human_mapping <- read.csv("./rnaseq/annotation/human_ensembl_ids.csv", stringsAsFactors=FALSE)

# Set some paths
microarray_directory <- "./microarray/processed_microarray/tissue/"
tcga_directory <- "./tcga/normalized/"

# Setup some variables
mouse_tissue_geneset_names <- "liver,kidney,spleen,lung,embryonic stem line Bruce4 p13,embryonic stem line V26 2 p16"
mouse_tissue_geneset_names <- unlist(strsplit(mouse_tissue_geneset_names, ","))
sample_mapping_myc <- c(liver_myc="HCC (MYC)", kidney_myc="RCC (MYC)", tall_myc="T-ALL (MYC)", eumyc_myc="BCL (MYC)", lung_myc="LAC (MYC)")
sample_mapping <- c(sample_mapping_myc, lung_mycras="LAC (MYC+KRAS)", lung_ras="LAC (KRAS)")
sample_mapping_extended <- c(sample_mapping, kidney_myc_preneo_early="RCC Early Pre-Neo (MYC)", kidney_myc_preneo_late="RCC Late Pre-Neo (MYC)")
microarray_sample_expression_control <- c(liver_myc="liver_control", kidney_myc="kidney_control", kidney_myc_preneo_early=c("kidney_myc",1,2,3), kidney_myc_preneo_late=c("kidney_myc",4,5,6), lung_myc="lung_control", lung_ras="lung_control", lung_mycras="lung_control")

# Setup some functions
do_enrichment <- function(db_name, sample_mapping, genes_de_up, genes_de_down=NULL, output_file=NULL, geneset_names=NULL) { # do enrichR analysis
    enrichment <- list()
    enrichment_data_up <- NULL
    enrichment_data_dn <- NULL
    if (is.null(sample_mapping)) {
      sample_mapping <- names(genes_de_up)
      names(sample_mapping) <- names(genes_de_up)
    }
    for (i in names(sample_mapping)) {
        if (typeof(genes_de_up) == "list") {
          genes_up_list <- genes_de_up[[i]]
        } else {
          genes_up_list <- genes_de_up[,i]
        }
        enrich_up <- enrichr(genes_up_list, c(db_name))[[db_name]]
        if (!is.null(genes_de_down)) {
            if (typeof(genes_de_down) == "list") {
              genes_down_list <- genes_de_down[[i]]
            } else {
              genes_down_list <- genes_de_down[,i]
            }
            enrich_dn <- enrichr(genes_down_list, c(db_name))[[db_name]]
            enrichment[[paste(sample_mapping[i], "UP", sep=" - ")]] <- enrich_up
            enrichment[[paste(sample_mapping[i], "DOWN", sep=" - ")]] <- enrich_dn
            if (!is.null(geneset_names)) {
                enrich_dn <- enrich_dn[enrich_dn[,"Term"] %in% geneset_names,]
            }
            enrich_dn <- enrich_dn[,c("Term","Odds.Ratio","P.value", "Adjusted.P.value")]
            colnames(enrich_dn) <- c("Term", paste(i, colnames(enrich_dn)[2:ncol(enrich_dn)], sep="_"))
        } else {
            enrich_dn <- NULL
            enrichment[[sample_mapping[i]]] <- enrich_up
        }
        if (!is.null(geneset_names)) {
            enrich_up <- enrich_up[enrich_up[,"Term"] %in% geneset_names,]
        }
        enrich_up <- enrich_up[,c("Term","Odds.Ratio","P.value", "Adjusted.P.value")]
        colnames(enrich_up) <- c("Term", paste(i, colnames(enrich_up)[2:ncol(enrich_up)], sep="_"))
        if (is.null(enrichment_data_up)) {
            enrichment_data_up <- enrich_up
            enrichment_data_dn <- enrich_dn
        } else {
            enrichment_data_up <- merge(enrichment_data_up, enrich_up, by="Term", all=TRUE)
            if (!is.null(genes_de_down)) {
                enrichment_data_dn <- merge(enrichment_data_dn, enrich_dn, by="Term", all=TRUE)
            }
        }
    }
    rownames(enrichment_data_up) <- enrichment_data_up$Term
    enrichment_data_up$Term <- NULL
    if (!is.null(geneset_names)) {
        enrichment_data_up <- enrichment_data_up[geneset_names,]
    }
    if (!is.null(genes_de_down)) {
        rownames(enrichment_data_dn) <- enrichment_data_dn$Term
        enrichment_data_dn$Term <- NULL
        if (!is.null(geneset_names)) {
            enrichment_data_dn <- enrichment_data_dn[geneset_names,]
        }
    }
    wb <- createWorkbook("Enrichment Results")
    for (i in names(enrichment)) {
        addWorksheet(wb, i)
        writeData(wb, sheet = i, enrichment[[i]], rowNames = FALSE)
    }
    saveWorkbook(wb, output_file, overwrite = TRUE)
    if (is.null(genes_de_down)) {
        return(enrichment_data_up)
    } else {
        return(list("up"=enrichment_data_up, "down"=enrichment_data_dn))
    }
}

# Load mouse tissue genesets
genesets_mouse_tissue <- list()
for (line in raw_genesets_mouse_tissue) {
    for (geneset_name in mouse_tissue_geneset_names) {
        if (startsWith(line, paste(geneset_name, "\t", sep=""))) {
            geneset_data <- unlist(strsplit(line, "\t"))
            geneset_data <- c(geneset_name, gsub("(.*),.*", "\\1", geneset_data[2:length(geneset_data)]))
            genesets_mouse_tissue[[geneset_name]] <- geneset_data
        }
    }
}
