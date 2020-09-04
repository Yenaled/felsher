library(ggplot2)
library(ggbeeswarm)
library(openxlsx)

output_folder <<- "./output/figures/"

plotSE <- function(fname, bed_up_file, bed_down_file, tissue_highlight, output_file, wb, chipseq_study, de_study, output_xlsx=TRUE) {
  bed_cols <- c("chr","start","end","gene","tissue")
  up_label <- "Upregulated genes"
  down_label <- "Downregulated genes"
  colors <- c("black","red","blue","purple")
  bed_up <- read.table(bed_up_file, sep="\t", stringsAsFactors=FALSE,header=FALSE, col.names=bed_cols)
  bed_down <- read.table(bed_down_file, sep="\t", stringsAsFactors=FALSE,header=FALSE, col.names=bed_cols)
  data_mat <- read.table(gzfile(fname), stringsAsFactors=FALSE, header=FALSE, sep="\t", skip=1)
  data_up <- merge(bed_up, data_mat, by.x=c("chr","start","end"), by.y=c("V1","V2","V3"), all=FALSE)
  data_down <- merge(bed_down, data_mat, by.x=c("chr","start","end"), by.y=c("V1","V2","V3"), all=FALSE)
  data_up_mean <- rowMeans(data_up[,9:ncol(data_up)])
  data_down_mean <- rowMeans(data_down[,9:ncol(data_down)])
  data <- data.frame(chr=c(data_up$chr, data_down$chr), start=c(data_up$start, data_down$start), end=c(data_up$end, data_down$end), mean=c(data_up_mean, data_down_mean), gene=c(data_up$gene, data_down$gene), tissue=c(data_up$tissue, data_down$tissue), regulation=c(rep(up_label, length(data_up_mean)), rep(down_label, length(data_down_mean))), stringsAsFactors=FALSE)
  data <- data[!duplicated(data),]
  data$color <- "0"
  data[data$tissue == tissue_highlight,"color"] <- "1"
  print(output_file)
  print(paste("Upregulated genes: ", length(unique(data[data$regulation == up_label,"gene"])), " genes (n = ", nrow(data[data$regulation == up_label,]), " superenhancers)", sep=""))
  print(paste("Downregulated genes: ", length(unique(data[data$regulation == down_label,"gene"])), " genes (n = ", nrow(data[data$regulation == down_label,]), " superenhancers)", sep=""))
  print("Upregulated genes top fold changes:")
  data_up <- data[data$regulation == up_label,]
  data_up <- data_up[order(-data_up$mean),c("chr", "start", "end", "gene","mean","tissue")]
  colnames(data_up) <- c("chr", "start", "end", "Gene", "H3K27ac Log2 Fold Change", "Superenhancer Tissue")
  print(data_up[1:10,c("Gene", "H3K27ac Log2 Fold Change", "Superenhancer Tissue")])
  title <- paste(chipseq_study, " (", de_study, " up genes SEs)", sep="")
  if (output_xlsx) {
    addWorksheet(wb, title)
    writeData(wb, sheet = title, data_up, rowNames = FALSE)
  }
  print("Downregulated genes top fold changes:")
  data_down <- data[data$regulation == down_label,]
  data_down <- data_down[order(data_down$mean),c("chr", "start", "end", "gene","mean","tissue")]
  colnames(data_down) <- c("chr", "start", "end", "Gene", "H3K27ac Log2 Fold Change", "Superenhancer Tissue")
  print(data_down[1:10,c("Gene", "H3K27ac Log2 Fold Change", "Superenhancer Tissue")])
  title <- paste(chipseq_study, " (", de_study, " down genes SEs)", sep="")
  if (output_xlsx) {
    addWorksheet(wb, title)
    writeData(wb, sheet = title, data_down, rowNames = FALSE)
  }
  print(t.test(data[data$regulation == up_label,"mean"],data[data$regulation == down_label,"mean"]))
  print("--------------")
  pdf(paste(output_folder, output_file, sep=""))
  p <- ggplot(data, aes(regulation, mean, color=color)) + geom_beeswarm(show.legend=FALSE,cex=0.8,size=1) + 
    scale_color_manual(values=colors) + theme_classic() + 
    theme(axis.title.x = element_blank(), axis.text=element_text(size=16,color="black"), axis.title=element_text(size=16,color="black"), axis.text.x = element_text(size = 18)) +
    geom_hline(yintercept=0, linetype="dashed") + 
    labs(y = expression(paste("H3K27ac ", Log[2]," Fold Change",sep=""))) + ylim(-2.15, 2.5)
  print(p)
  device <- dev.off()
}

wb <- createWorkbook("Superenhancer H3K27ac profiles")
sink("./output/se_log.txt")

plotSE("./output/mat/hcc_h3k27ac_dbsuper_log2FC.mat.gz", "dbSUPER_hcc_up.bed", "dbSUPER_hcc_down.bed", "Liver", "plot_hcc_h3k27ac_dbsuper.pdf", wb, "HCC", "HCC")
plotSE("./output/mat/eumyc_h3k27ac_dbsuper_log2FC.mat.gz", "dbSUPER_eumyc_up.bed", "dbSUPER_eumyc_down.bed", "Spleen", "plot_eumyc_h3k27ac_dbsuper.pdf", wb, "BCL", "BCL")

# Flip genes (e.g. looking at eumyc genes in HCC chip-seq)
plotSE("./output/mat/hcc_h3k27ac_dbsuper_log2FC.mat.gz", "dbSUPER_eumyc_up.bed", "dbSUPER_eumyc_down.bed", "None", "plot_hcc_h3k27ac_dbsuper_using_eumycgenes.pdf", wb, "HCC", "BCL", FALSE)
plotSE("./output/mat/eumyc_h3k27ac_dbsuper_log2FC.mat.gz", "dbSUPER_hcc_up.bed", "dbSUPER_hcc_down.bed", "None", "plot_eumyc_h3k27ac_dbsuper_using_hccgenes.pdf", wb, "BCL", "HCC", FALSE)

sink()
saveWorkbook(wb, "./output/se_h3k27ac_data.xlsx", overwrite = TRUE)

