###############################################################
### File: overlap.r
### Description: Performs overlap analysis of results from 
###              microarray and RNA-seq differential gene
###              expression analysis.
### Usage: Rscript overlap.r <log2FC> <fdr> <output_dir> <labels> <display_labels> <files> <microarray_annotation> <rnaseq_annotation>
### Written by Delaney Sullivan
###############################################################

# Rscript overlap.r 2 0.05 output/mouse_de/ liver_myc,kidney_myc,lung_myc,lung_mycras,lung_ras,tall_myc,eumyc_myc "HCC,RCC,LAC (MYC),LAC (MYC+KRAS),LAC (KRAS),T-ALL,BCL" microarray/processed_microarray/tissue/liver_diffexpr_gene.csv,microarray/processed_microarray/tissue/kidney_diffexpr_gene.csv,microarray/processed_microarray/tissue/lung_diffexpr_gene.csv,microarray/processed_microarray/tissue/lung_mycras_diffexpr_gene.csv,microarray/processed_microarray/tissue/lung_ras_diffexpr_gene.csv,rnaseq/tall_rnaseq/kallisto_aligned/tallmycon_sleuth_results_genes.csv,rnaseq/eumyc_rnaseq/kallisto_aligned/eumycon_sleuth_results_genes.csv microarray/annotation/MouseWG-6v2.csv rnaseq/annotation/gencode_GRCm38_vM15.csv

# Load R packages
pkgs <- c("UpSetR", "VennDiagram", "openxlsx")
invisible(lapply(pkgs, function(x) suppressWarnings(suppressMessages(library(x, character.only=TRUE)))))

# Read command-line arguments
args <- commandArgs(TRUE)
lfcThreshold <- log2(as.numeric(args[1]))
fdrThreshold <- as.numeric(args[2])
outputFolder <- args[3]
labels <- unlist(strsplit(args[4], ","))
display_labels <- unlist(strsplit(args[5], ","))
files <- unlist(strsplit(args[6], ","))
microarray_annotation <- read.csv(args[7], stringsAsFactors=FALSE, header=TRUE)
microarray_annotation <- microarray_annotation[!is.na(microarray_annotation$ID),]
microarray_annotation <- microarray_annotation[!duplicated(microarray_annotation$ID),]
rnaseq_annotation <- read.csv(args[8], stringsAsFactors=FALSE, header=TRUE)
rnaseq_annotation <- rnaseq_annotation[!is.na(rnaseq_annotation$ID),]

if (!all(c(length(labels), length(display_labels)) == length(files))) {
    stop("<labels>, <display_labels>, and <files> must all have same number of comma-separated values")
}

data <- NULL
data_padj <- NULL
for (i in 1:length(files)) {
    f <- files[i]
    label <- labels[i]
    df <- read.csv(f, stringsAsFactors=FALSE, header=TRUE)
    if ("b" %in% colnames(df)) { # Sleuth RNA-seq data
        df <- df[,c("gene_name","b","qval")]
        df <- merge(rnaseq_annotation, df, by.x="Symbol", by.y="gene_name")
        df$Symbol <- NULL
        rownames(df) <- df$ID
        df$ID <- NULL
        colnames(df) <- c("logFC","padj")
    } else { # Limma microarray data
        df <- df[,c("ID","logFC","adj.P.Val")]
        rownames(df) <- df$ID
        df$ID <- NULL
        colnames(df) <- c("logFC","padj")
    }

    df_logFC <- df[,"logFC",drop=FALSE]
    colnames(df_logFC) <- label
    df_padj <- df[,"padj",drop=FALSE]
    colnames(df_padj) <- label

    if (is.null(data)) {
        data <- df_logFC
        data_padj <- df_padj
    } else {
        data <- transform(merge(data,df_logFC,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
        data_padj <- transform(merge(data_padj,df_padj,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
    }
}

# Filter differentially expressed genes
data_significant_up <- data
data_significant_dn <- data
data_significant_up[!is.na(data) & data > lfcThreshold & !is.na(data_padj) & data_padj < fdrThreshold] <- TRUE
data_significant_up[is.na(data) | data <= lfcThreshold | is.na(data_padj) | data_padj >= fdrThreshold] <- FALSE
data_significant_dn[!is.na(data) & data < -lfcThreshold & !is.na(data_padj) & data_padj < fdrThreshold] <- TRUE
data_significant_dn[is.na(data) | data >= -lfcThreshold | is.na(data_padj) | data_padj >= fdrThreshold] <- FALSE

# Create excel spreadsheet table
data_order <- data_significant_up*2 + data_significant_dn
data_ordering <- order(rowSums(data_order != 0), apply(data_order, 1, function(x) { sort(table(x),decreasing=TRUE)[1] }), rowSums(data_order == 2), rowSums(data_order == 1), apply(data_order, 1, paste, collapse="-"), rowSums(!is.na(data)), decreasing=TRUE)
data_order <- data_order[data_ordering,]
data_excel <- data[data_ordering,]
colnames(data_excel) <- display_labels
annotation_intersection <- intersect(microarray_annotation$ID, rnaseq_annotation$ID)
microarray_remaining <- microarray_annotation[!(microarray_annotation$ID %in% annotation_intersection),]
microarray_remaining <- microarray_remaining[!duplicated(microarray_remaining$ID),c("Symbol","ID")]
annotation_all <- rbind(rnaseq_annotation, microarray_remaining)
data_excel <- merge(annotation_all, data_excel, by.x="ID", by.y=0, all.x=FALSE, all.y=TRUE)
data_excel <- data_excel[match(rownames(data_order), data_excel$ID),]
wb <- createWorkbook("Differentially Expressed Genes")
addWorksheet(wb, "Differentially Expressed Genes") 
writeData(wb, sheet = 1, data_excel, rowNames = FALSE)
styleUp <- createStyle(fgFill = "red")
styleDn <- createStyle(fgFill = "blue")
for (i in 1:ncol(data_order)) {
    addStyle(wb, sheet = 1, styleUp, rows = which(data_order[,i] == 2) + 1, cols = i+2)
    addStyle(wb, sheet = 1, styleDn, rows = which(data_order[,i] == 1) + 1, cols = i+2)
}
saveWorkbook(wb, paste(outputFolder, "/de_table.xlsx", sep=""), overwrite = TRUE)

# Prepare to make UpSet plots and also output information about DE genes
sink(file=paste(outputFolder, "/log.txt", sep=""), split=TRUE)
listInputUp <- list()
listInputDn <- list()
de_genes_list_up <- list()
de_genes_list_dn <- list()
for (i in 1:ncol(data)) {
    de_genes_list_up[[labels[i]]] <- rownames(data_significant_up[data_significant_up[,i] == TRUE,i,drop=FALSE])
    de_genes_list_dn[[labels[i]]] <- rownames(data_significant_dn[data_significant_dn[,i] == TRUE,i,drop=FALSE])
    display_label <- display_labels[i] # Label to display for UpSet plot
    display_label <- unlist(strsplit(display_label, " \\(.*\\)")) # Remove stuff in parantheses in display label
    if (!(display_label %in% names(listInputUp))) {
        listInputUp[[display_label]] <- de_genes_list_up[[labels[i]]]
        listInputDn[[display_label]] <- de_genes_list_dn[[labels[i]]]
    }
    print(paste(display_labels[i], "-", labels[i], "-", (sum(data_significant_up[,i]) + sum(data_significant_dn[,i])), "DE genes:", sum(data_significant_up[,i]), "up and", sum(data_significant_dn[,i]), "down", sep=" "))
}
print(paste("Number of genes total:", nrow(data)))
sink()
outputList <- function(de_genes_list, output_file, annotation=NULL) {
    max.length <- max(sapply(de_genes_list, length))
    df_genes <- lapply(de_genes_list, function(v) { c(v, rep(NA, max.length-length(v)))})
    df_genes <- t(do.call(rbind, df_genes))
    if (!is.null(annotation)) {
        col_names <- colnames(df_genes)
        df_genes <- as.data.frame(matrix(with(annotation, Symbol[match(df_genes[ ,1:ncol(df_genes)], ID)]), nrow = nrow(df_genes)))
        colnames(df_genes) <- col_names
    }
    write.table(df_genes, file=output_file, sep="\t", quote=FALSE, row.names=FALSE, na="")
}
outputList(de_genes_list_up, paste(outputFolder, "/de_genes_up_ids.txt", sep=""))
outputList(de_genes_list_dn, paste(outputFolder, "/de_genes_down_ids.txt", sep=""))
outputList(de_genes_list_up, paste(outputFolder, "/de_genes_up_symbols.txt", sep=""), annotation_all)
outputList(de_genes_list_dn, paste(outputFolder, "/de_genes_down_symbols.txt", sep=""), annotation_all)

# Make UpSet plots
pdf(paste(outputFolder, "/upset_up.pdf", sep=""), width=10, height=7.5, onefile=FALSE)
upset(fromList(listInputUp), mb.ratio=c(0.65, 0.35), text.scale=c(2, 2, 1.5, 1.5, 2, 1.6), scale.intersections="identity", nsets=length(names(listInputUp)), order.by="degree", point.size=3.5, line.size=1.5)
device <- dev.off()
pdf(paste(outputFolder, "/upset_down.pdf", sep=""), width=10, height=7.5, onefile=FALSE)
upset(fromList(listInputDn), mb.ratio=c(0.65, 0.35), text.scale=c(2, 2, 1.5, 1.5, 2, 1.6), scale.intersections="identity", nsets=length(names(listInputDn)), order.by="degree", point.size=3.5, line.size=1.5)
device <- dev.off()


# Output summary table
if (length(listInputUp) == 5) { # If we have 5 sets, draw a Venn Diagram as well
    pdf(paste(outputFolder, "/venn.pdf", sep=""), width=14, height=7)
    listInputUp <<- listInputUp
    listInputDn <<- listInputDn
    overlap <- function(...) {
        set_indices <- list(...)
        total <- NULL
        for (i in set_indices) {
           if (is.null(total)) {
               total <- listInputUp[[i]]
           } else {
               total <- intersect(total, listInputUp[[i]])
           }
        }
        total2 <- NULL
        for (i in set_indices) {
           if (is.null(total2)) {
               total2 <- listInputDn[[i]]
           } else {
               total2 <- intersect(total2, listInputDn[[i]])
           }
        }
        return(length(total) + length(total2))
    }
    draw.quintuple.venn(overlap(1), overlap(2), overlap(3), overlap(4), overlap(5), overlap(1,2), overlap(1,3), overlap(1,4), overlap(1,5),
        overlap(2,3), overlap(2,4), overlap(2,5), overlap(3,4), overlap(3,5), overlap(4,5), overlap(1,2,3), overlap(1,2,4), overlap(1,2,5), overlap(1,3,4),
        overlap(1,3,5), overlap(1,4,5), overlap(2,3,4), overlap(2,3,5), overlap(2,4,5), overlap(3,4,5), overlap(1,2,3,4), overlap(1,2,3,5),
        overlap(1,2,4,5), overlap(1,3,4,5), overlap(2,3,4,5), overlap(1,2,3,4,5), category = names(listInputUp),
        lwd = rep(2, 5), lty = rep("solid", 5), col =
        rep("black", 5), fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"), alpha = rep(0.5, 5),
        label.col = rep("black", 31), cex = c(rep(2.15, 5), rep(1.75, 20), rep(2.3, 5), 2.6),
        fontface = c(rep("plain", 25), rep("bold", 6)), fontfamily = rep("serif",
        31), cat.pos = c(0, 287.5+50, 215+25, 145, 70-50), cat.dist =
        rep(0.22, 5), cat.col = rep("black", 5), cat.cex =
        rep(2.5, 5), cat.fontface = rep("plain", 5),
        cat.fontfamily = rep("serif", 5), cat.just =
        rep(list(c(0.5, 0.5)), 5), rotation.degree = 0,
        rotation.centre = c(0.5, 0.5), ind = TRUE, cex.prop =
        NULL, print.mode = "raw", sigdigs = 3, direct.area =
        FALSE, area.vector = 0, margin=0.022)
    device <- dev.off()
}
                 
