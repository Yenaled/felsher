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
ensembl_human_symbols_mapping <- read.csv("./rnaseq/annotation/hg19_ensembl.csv", stringsAsFactors=FALSE)

# Set some paths
microarray_directory <- "./microarray/processed_microarray/tissue/"
tcga_lfc_file <- "./tcga/TCGA.processed.tumor_normal_lfc.tsv.gz"
tcga_info_file <- "./tcga/TCGA.processed.info_corrected.csv"
tcga_rra_file <- "./output/tcga_correlation/myc_rra.csv"

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

# Manually update a few genes in the human Ensembl-to-symbol mapping, which is slightly out-of-date (note: NOT COMPREHENSIVE!)
updated_ensembl_mapping <- list(ENSG00000275700="AATF",
                                ENSG00000276043="UHRF1",
                                ENSG00000276234="TADA2A",
                                ENSG00000170468="RIOX1",
                                ENSG00000275714="HIST1H3A",
                                ENSG00000249859="PVT1",
                                ENSG00000261236="BOP1",
                                ENSG00000213553="RPLP0P6",
                                ENSG00000213866="YBX1P10",
                                ENSG00000234851="RPL23AP42",
                                ENSG00000233476="EEF1A1P6",
                                ENSG00000196205="EEF1A1P5",
                                ENSG00000281398="SNHG4",
                                ENSG00000214485="RPL7P1",
                                ENSG00000137970="RPL7P9",
                                ENSG00000249353="NPM1P27",
                                ENSG00000224861="YBX1P1",
                                ENSG00000214199="EEF1A1P12",
                                ENSG00000249264="EEF1A1P9",
                                ENSG00000250182="EEF1A1P13",
                                ENSG00000249855="EEF1A1P19",
                                ENSG00000204745="AC083899.1",
                                ENSG00000172974="VDAC2P5",
                                ENSG00000163597="SNHG16",
                                ENSG00000234964="FABP5P7",
                                ENSG00000234743="EIF5AP4",
                                ENSG00000278845="MRPL45",
                                ENSG00000220842="RPL21P16",
                                ENSG00000234741="GAS5",
                                ENSG00000255717="SNHG1",
                                ENSG00000189343="RPS2P46",
                                ENSG00000204253="HNRNPCP2",
                                ENSG00000224578="HNRNPA1P48",
                                ENSG00000243199="AC115223.1",
                                ENSG00000232472="EEF1B2P3",
                                ENSG00000139239="RPL14P1",
                                ENSG00000250899="AC125807.2",
                                ENSG00000183604="SMG1P5",
                                ENSG00000228649="SNHG26",
                                ENSG00000213300="HNRNPA3P6",
                                ENSG00000213131="YWHAZP4",
                                ENSG00000241506="PSMC1P1",
                                ENSG00000232024="LSM12P1")
updated_ensembl_mapping <- data.frame(Ensembl=names(updated_ensembl_mapping), Symbol=unlist(unname(updated_ensembl_mapping)))
ensembl_human_symbols_mapping <- rbind(updated_ensembl_mapping, ensembl_human_symbols_mapping)
ensembl_human_symbols_mapping <- ensembl_human_symbols_mapping[!duplicated(ensembl_human_symbols_mapping$Ensembl),]
