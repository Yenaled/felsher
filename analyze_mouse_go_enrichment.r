###############################################################
### File: analyze_mouse_go_enrichment.r
### Description: Performs mouse GO enrichment analysis
### Usage #1: Rscript analyze_mouse_go_enrichment.r
### Usage #2: Rscript analyze_mouse_go_enrichment.r output/enrichr_mouse_go/go_gsea.xlsx data/go_term_collapse.csv
### Written by Delaney Sullivan
###############################################################

source("analysis_setup.r")
output_dir <- "./output/enrichr_mouse_go/"
threshold_num_terms <- 100 # The number of GO terms in each GSEA to be considered "top"

# Get command-line arguments
args <- commandArgs(TRUE)
input_gsea <- NULL # GSEA results previously computed
go_term_collapse <- NULL # Collapsed GO term annotation (for heatmap display)
if (length(args) >= 2) {
  input_gsea <- args[1]
  go_term_collapse <- read.csv(args[2], stringsAsFactors=FALSE, check.names=FALSE)
  go_term_collapse$Other <- NULL
}

# Function to get top GO terms
selectTopGO <- function(data, threshold_num_terms) {
  top_terms <- c()
  for (i in colnames(data)) {
    curr_row <- data[!is.na(data[,i]), ]
    curr_row <- curr_row[order(curr_row[,i]),]
    curr_row <- curr_row[1:threshold_num_terms,]
    curr_row <- curr_row[complete.cases(curr_row),]
    top_terms <- c(top_terms, rownames(curr_row))
  }
  top_terms <- unique(top_terms)
  return(top_terms)
}

# Read in differentially expressed genes
de_table <- read.xlsx("output/mouse_de/de_table.xlsx", sep.names=" ", check.names=FALSE)
de_table$ID <- NULL
rownames(de_table) <- de_table$Symbol
de_table$Symbol <- NULL
# Skip the RAS samples
de_table[,sample_mapping["lung_ras"]] <- NULL
de_table[,sample_mapping["lung_mycras"]] <- NULL
colnames(de_table) <- unlist(strsplit(colnames(de_table), " \\(.*\\)")) # Remove stuff in parantheses in display label

# Read in GO gene lists
genesets_mouse_go <- gmtPathways("data/GO_Biological_Process_2018.txt")

if (is.null(input_gsea)) { # If we have yet to perform GSEA (i.e. GSEA results file not supplied)
  # Do GSEA and write output to .xlsx file
  wb <- createWorkbook("GSEA Results")
  gsea_results <- NULL
  for (i in colnames(de_table)) {
    expression_rnk_list <- setNames(de_table[,i], toupper(rownames(de_table)))
    expression_rnk_list <- expression_rnk_list[!is.na(expression_rnk_list)]
    res <- fgsea(pathways=genesets_mouse_go, stats=expression_rnk_list, gseaParam=1, nperm=250000)
    res <- as.data.frame(res)[,c("pathway","ES", "NES", "pval", "padj")]
    res <- res[order(res$pval, abs(res$NES)),]
    addWorksheet(wb, i)
    writeData(wb, sheet = i, res, rowNames = FALSE)
    colnames(res) <- c("pathway", paste(i, colnames(res)[2:ncol(res)], sep="_"))
    if (is.null(gsea_results)) {
      gsea_results <- res
    } else {
      gsea_results <- merge(gsea_results, res, by="pathway", all=TRUE)
    }
  }
  saveWorkbook(wb, paste(output_dir, "go_gsea.xlsx", sep=""), overwrite = TRUE)
  rownames(gsea_results) <- gsea_results$pathway
  
  # Get p-value data
  gsea_results_pval <- gsea_results[,grepl("_pval$", colnames(gsea_results))]
  colnames(gsea_results_pval) <- substr(colnames(gsea_results_pval),1,nchar(colnames(gsea_results_pval))-5)
  gsea_results_pval[is.na(gsea_results_pval)] <- 1
  
  # Get top terms based on p-value
  topTerms <- selectTopGO(gsea_results_pval, threshold_num_terms)
  writeLines(topTerms, paste(output_dir, "topterms.txt", sep=""))
  quit()
}

# Do tissue GO enrichment
tissue_go <- do_enrichment("GO_Biological_Process_2018", NULL, genesets_mouse_tissue, output_file=paste(output_dir, "tissue_go.xlsx", sep=""))
tissue_go_or <- tissue_go[,grepl("_Odds.Ratio$", colnames(tissue_go))]
colnames(tissue_go_or) <- substr(colnames(tissue_go_or),1,nchar(colnames(tissue_go_or))-11)
tissue_go_or <- tissue_go_or[,mouse_tissue_geneset_names]
tissue_go_or[is.na(tissue_go_or)] <- 0

gsea_results <- NULL
for (i in colnames(de_table)) {
  res <- read.xlsx(input_gsea, sheet=i, sep.names=" ", check.names=FALSE)
  colnames(res) <- c("pathway", paste(i, colnames(res)[2:ncol(res)], sep="_"))
  if (is.null(gsea_results)) {
    gsea_results <- res
  } else {
    gsea_results <- merge(gsea_results, res, by="pathway", all=TRUE)
  }
}
rownames(gsea_results) <- gsea_results$pathway

# Get the padj and NES data into separate data frames
gsea_results_padj <- gsea_results[,grepl("_padj$", colnames(gsea_results))]
colnames(gsea_results_padj) <- substr(colnames(gsea_results_padj),1,nchar(colnames(gsea_results_padj))-5)
gsea_results <- gsea_results[,grepl("_NES$", colnames(gsea_results))]
colnames(gsea_results) <- substr(colnames(gsea_results),1,nchar(colnames(gsea_results))-4)

# Fix the NA's in the padj and NES
gsea_results_padj[is.na(gsea_results_padj)] <- 1
gsea_results[is.na(gsea_results)] <- 0

# Select the GO terms we want to analyze
topTerms <- unname(unlist(go_term_collapse))
topTerms <- topTerms[!is.na(topTerms) & topTerms != ""]
print(paste("Processing", length(topTerms), "terms"))
gsea_results_padj <- gsea_results_padj[rownames(gsea_results_padj) %in% topTerms,]
gsea_results <- gsea_results[rownames(gsea_results) %in% topTerms,]

# Convert GSEA padj into signed log10 padj
for (i in rownames(gsea_results_padj)) {
  for (j in colnames(gsea_results_padj)) {
    gsea_results[i,j] <- sign(gsea_results[i,j])*-log10(gsea_results_padj[i,j])
  }
}

# Cluster GO terms by supplied category
go_cluster_assignment <- list()
for (go_term in topTerms) {
  go_term_cluster <- apply(go_term_collapse, 1, function(x) match(go_term, x))
  go_term_cluster <- colnames(go_term_collapse)[go_term_cluster[!is.na(go_term_cluster)][1]]
  if (!is.na(go_term_cluster)) {
    go_cluster_assignment[go_term] <- go_term_cluster
  }
}
gsea_results <- gsea_results[rownames(gsea_results) %in% names(go_cluster_assignment),]
gsea_results <- as.matrix(gsea_results[names(go_cluster_assignment),])

# Function to set color scheme for heatmap
col_fun <- function(x, breaks) {
  a <- rep("white", length(x))
  cutoff <- 0.25
  a[abs(x) < -log10(cutoff)] <- "white"
  a[x >= -log10(cutoff)] <- "#F6CAC8"
  a[x <= log10(cutoff)] <- "#CEC9FF"
  cutoff <- 0.10
  a[x >= -log10(cutoff)] <- "#EF908F"
  a[x <= log10(cutoff)] <- "#958EF7"
  cutoff <- 0.05
  a[x >= -log10(cutoff)] <- "#EB4848"
  a[x <= log10(cutoff)] <- "#3E3EF5"
  cutoff <- 0.01
  a[x >= -log10(cutoff)] <- "red"
  a[x <= log10(cutoff)] <- "blue"
  return(a)
}
attr(col_fun,'breaks') <- c(-2, 0, 2)

# Setup heatmap
ha_plot = HeatmapAnnotation(graph=anno_empty(border=TRUE), height=unit(2, "cm"), gap=unit(c(0.2,0), "cm"))
hm <- Heatmap(as.matrix(t(gsea_results)), bottom_annotation = ha_plot, column_split=factor(paste0(go_cluster_assignment), levels=(colnames(go_term_collapse))), border=TRUE, show_row_dend=FALSE, show_column_names=FALSE, column_gap=unit(0.1, "mm"), column_dend_height=unit(7.5, "mm"), row_names_side="left", row_names_gp = gpar(fontsize = 9), height=unit(20, "mm"), width=unit(200, "mm"), rect_gp = gpar(col = "black", lwd = 0.15), column_title_gp = gpar(fontsize = 3), col=col_fun, cluster_column_slices=FALSE, show_heatmap_legend=FALSE)

# Function to draw bar graph below the heatmap
drawBars <- function(values, num_cells) {
  pos <- floor(num_cells / 2) - 1.65
  box_width <- 1
  colors <- c("black", "red","blue","purple","#A9A87E","#A37C4B")
  for (v in 1:length(values)) {
    grid.rect(x = pos+(v-1)*box_width, y = 0, height=values[v], width=box_width, just="bottom", gp = gpar(fill = colors[v], lwd=0), default.units = "native")
  }
}

# Draw heatmap and the bar graph
pdf(paste(output_dir, "heatmap_go.pdf", sep=""), width=16, height=8)
pushViewport(viewport(layout=grid.layout(nr=1, nc=1)))
pushViewport(viewport(layout.pos.row=1, layout.pos.col=1))
draw(hm, newpage=FALSE)
for(i in 1:ncol(go_term_collapse)) {
  decorate_annotation("graph", slice = i, {
    num_cells <- sum(go_cluster_assignment == colnames(go_term_collapse)[i])
    pushViewport(viewport(xscale=c(0, num_cells), yscale = c(0, 7)))
    terms_in_cluster <- names(go_cluster_assignment[(go_cluster_assignment == colnames(go_term_collapse)[i])])
    drawBars(apply(tissue_go_or[rownames(tissue_go_or) %in% terms_in_cluster,], 2, median), num_cells)
    if (i == 1) {
      grid.yaxis(at = c(0, 3.5, 7), gp=gpar(fontsize = 8))
    }
    popViewport()
  })
}
upViewport()
upViewport()
heatmap_order <- column_order(hm)
device <- dev.off()

# Write out GO terms from heatmap
heatmap_go_terms <- data.frame(Category=character(), Term=character(), stringsAsFactors=FALSE)
for (term_category in names(heatmap_order)) {
  terms_to_add <- rownames(gsea_results)[heatmap_order[[term_category]]]
  heatmap_go_terms <- rbind(heatmap_go_terms, data.frame(Category=rep(term_category, length(terms_to_add)), Term=terms_to_add))
}
write.csv(heatmap_go_terms, file=paste(output_dir, "heatmap_go_terms.csv", sep=""), row.names=FALSE)
