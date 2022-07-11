require(ComplexHeatmap)
# Select genes in the 67-gene signature that are part of both (i.e. within the intersection of) "embryonic stem line V26 2 p16" and "embryonic stem line Bruce4 p13"
genes <- "APEX1,GART,SRM,RUVBL1,NOP58,MYBBP1A,IMPDH2,NCL,NOP2,PWP1,WDR43,RRP9,TAF1D,RCC1,PRMT5,WDR75,RRP1B,GNL3,RSL1D1,CIRH1A,POLR1B,GRWD1,GEMIN4,GEMIN5"
genes <- unlist(strsplit(genes, ","))
negative_genes <- "ADH1"
negative_genes <- unlist(strsplit(toupper(negative_genes),","))
data <- read.table("MYC_off_vs_On_All_p05.txt", stringsAsFactors=F, sep="\t", skip=2, header=T, row.names=1)
data <- data[rownames(data) %in% c(negative_genes,genes),]
data <- data[intersect(c(negative_genes,genes), rownames(data)),1:18]
data <- data[,c(2,1,3,7,9,8,4,6,5, 10,11,12,13,15,14,16,18,17)] # Select first 18 columns (mildly rearranging them for clustering display)
rownames(data) <- tolower(rownames(data))
rownames(data) <- paste(toupper(substring(rownames(data), 1, 1)), substring(rownames(data), 2), sep = "")
hm <- Heatmap(as.matrix(data), border=TRUE, show_row_dend=FALSE, show_column_dend=FALSE, show_column_names=FALSE, column_gap=unit(0.1, "mm"), row_names_side="left", row_names_gp = gpar(fontsize = 9), height=unit(80, "mm"), width=unit(100, "mm"), rect_gp = gpar(col = "black", lwd = 0.15), column_title_gp = gpar(fontsize = 3), cluster_column_slices=FALSE, cluster_columns=FALSE, cluster_rows=FALSE, cluster_row_slices=FALSE, show_heatmap_legend=TRUE)
draw(hm, newpage=FALSE)



