load("aggregated_genesets_myc_up.RData")
ensembl_mapping <- read.csv("../../data/rnaseq/annotation/hg19_ensembl.csv", header=TRUE, stringsAsFactors=FALSE)
updated_ensembl <- list(ENSG00000265479="DTX2P1-UPK3BP1-PMS2P11",
                        ENSG00000196668="LINC00173",
                        ENSG00000199004="MIR21",
                        ENSG00000255717="SNHG1",
                        ENSG00000224411="HSP90AA2P",
                        ENSG00000211592="IGKC",
                        ENSG00000226950="DANCR",
                        ENSG00000275993="SIK1B",
                        ENSG00000196756="SNHG17",
                        ENSG00000273808="ABCC6P2",
                        ENSG00000280924="LINC00628",
                        ENSG00000146001="PCDHB18P")
ensembl_mapping <- rbind(ensembl_mapping, data.frame(Ensembl=names(updated_ensembl), Symbol=unlist(updated_ensembl)))
ensembl_mapping[ensembl_mapping$Symbol == "AATF","Ensembl"] <- "ENSG00000275700"
ensembl_mapping[ensembl_mapping$Symbol == "MARCKS","Ensembl"] <- "ENSG00000277443"
ensembl_mapping[ensembl_mapping$Symbol == "RPS17","Ensembl"] <- "ENSG00000182774"
ensembl_mapping[ensembl_mapping$Symbol == "ABCC6P2","Ensembl"] <- "ENSG00000255277"
ensembl_mapping[ensembl_mapping$Symbol == "PIGW","Ensembl"] <- "ENSG00000277161"

all_genes <- ensembl_mapping$Symbol

for (d in names(data)) {
  # Remove genes:
  data[[d]] <- (data[[d]])[data[[d]] != "MYC"]
  # Change genes:
  data[[d]][data[[d]] == "AC005255.1"] <- "OR7C1"
  data[[d]][data[[d]] == "AC114267.1"] <- "OR10H1"
  data[[d]][data[[d]] == "AC010605.1"] <- "FBXO27"
  data[[d]][data[[d]] == "ARL17P1"] <- "ARL17A"
  data[[d]][data[[d]] == "ADGRE2"] <- "EMR2"
  data[[d]][data[[d]] == "ABCC6P2"] <- "ABCC6P2" 
  data[[d]][data[[d]] == "ADGRF1"] <- "GPR110"
  #data[[d]][data[[d]] == "C14orf166B"] <- "LRRC74A"
  data[[d]][data[[d]] == "COP1"] <- "RFWD2"
  data[[d]][data[[d]] == "COQ8A"] <- "ADCK3"
  data[[d]][data[[d]] == "DEPP1"] <- "C10orf10"
  data[[d]][data[[d]] == "EEF1AKMT1"] <- "N6AMT2"
  data[[d]][data[[d]] == "CFAP70"] <- "TTC18"
  #data[[d]][data[[d]] == "TTC18"] <- "CFAP70"
  data[[d]][data[[d]] == "MTERF3"] <- "MTERFD1"
  #data[[d]][data[[d]] == "MTERFD1"] <- "MTERF3"
  data[[d]][data[[d]] == "SLC25A25-AS1"] <- "RP11-395P17.3"
  #data[[d]][data[[d]] == "RP11-395P17.3"] <- "SLC25A25-AS1"
  data[[d]][data[[d]] == "NFIA-AS1"] <- "RP5-833A20.1"
  #data[[d]][data[[d]] == "RP5-833A20.1"] <- "NFIA-AS1"
  data[[d]][data[[d]] == "RIOX2"] <- "MINA"
  #data[[d]][data[[d]] == "MINA"] <- "RIOX2"
  data[[d]][data[[d]] == "SINHCAF"] <- "FAM60A"
  #data[[d]][data[[d]] == "FAM60A"] <- "SINHCAF"
  data[[d]][data[[d]] == "P3H4"] <- "LEPREL4"
  #data[[d]][data[[d]] == "LEPREL4"] <- "P3H4"
  data[[d]][data[[d]] == "FAM86C"] <- "FAM86C1"
  data[[d]][data[[d]] == "GPAT3"] <- "AGPAT9"
  data[[d]][data[[d]] == "HSAJ2425"] <- "HSAJ2425"  
  data[[d]][data[[d]] == "PCDHB18P"] <- "PCDHB18P"  
  data[[d]][data[[d]] == "KIAA0738"] <- "TCAF1"
  data[[d]][data[[d]] == "KRTAP2-1"] <- "KRTAP2-1"  
  data[[d]][data[[d]] == "LINC00628"] <- "LINC00628"  
  data[[d]][data[[d]] == "LINC00842"] <- "LINC00842"  
  data[[d]][data[[d]] == "LOC145945"] <- "LOC145945"  
  data[[d]][data[[d]] == "LOC163131"] <- "ZNF780B"
  data[[d]][data[[d]] == "LOC283849"] <- "EXOC3L1"
  data[[d]][data[[d]] == "LOC285181"] <- "LOC285181"  
  data[[d]][data[[d]] == "LOC390705"] <- "LOC390705"  
  data[[d]][data[[d]] == "LOC439962"] <- "LOC439962"  
  data[[d]][data[[d]] == "LOC644005"] <- "LOC644005"  
  data[[d]][data[[d]] == "LOC642441"] <- "LOC642441"  
  data[[d]][data[[d]] == "LOC645238"] <- "LOC645238"  
  data[[d]][data[[d]] == "LOC646014"] <- "LOC646014"  
  data[[d]][data[[d]] == "LOC653605"] <- "LOC653605"  
  data[[d]][data[[d]] == "LOC692247"] <- "LOC692247"  
  data[[d]][data[[d]] == "MGC20647"] <- "MGC20647"  
  data[[d]][data[[d]] == "MGC34761"] <- "CLEC18C"
  data[[d]][data[[d]] == "MGC5457"] <- "MGC5457"  
  data[[d]][data[[d]] == "MGC5566"] <- "LINC01260"
  data[[d]][data[[d]] == "MIGA2"] <- "FAM73B"
  data[[d]][data[[d]] == "PCOTH"] <- "C1QTNF9B-AS1"
  data[[d]][data[[d]] == "PRELID3A"] <- "SLMO1"
  data[[d]][data[[d]] == "PRO2852"] <- "PRO2852"  
  data[[d]][data[[d]] == "PRO2900"] <- "HDLBP"
  data[[d]][data[[d]] == "RP5-1119A7.4"] <- "FOXRED2"
  data[[d]][data[[d]] == "ZNRD2"] <- "SSSCA1"
  data[[d]][data[[d]] == "CSDA"] <- "YBX3"
  #data[[d]][data[[d]] == "B3GNT1"] <- "B4GAT1"
  data[[d]][data[[d]] == "C10orf26"] <- "WBP1L"
  data[[d]][data[[d]] == "C12orf24"] <- "FAM216A"
  data[[d]][data[[d]] == "C18orf1"] <- "LDLRAD4"
  data[[d]][data[[d]] == "CXCR7"] <- "ACKR3"
  data[[d]][data[[d]] == "PGAM2"] <- "PGAM2" 
  data[[d]][data[[d]] == "PPPDE2"] <- "DESI1"
  data[[d]][data[[d]] == "B7H6"] <- "NCR3LG1"
  data[[d]][data[[d]] == "PEO1"] <- "TWNK"
  data[[d]][data[[d]] == "TWNK"] <- "C10orf2"
  #data[[d]][data[[d]] == "C10orf2"] <- "TWNK"
  #data[[d]][data[[d]] == "C11orf48"] <- "LBHD1"
  #data[[d]][data[[d]] == "C1orf51"] <- "CIART"
  #data[[d]][data[[d]] == "DIEXF"] <- "UTP25"
  data[[d]][data[[d]] == "FAM100A"] <- "UBALD1"
  data[[d]][data[[d]] == "FLJ10661"] <- "FAM86C1"
  data[[d]][data[[d]] == "KIAA1609"] <- "TLDC1"
  data[[d]][data[[d]] == "LOC646762"] <- "LOC646762" 
  data[[d]][data[[d]] == "LOC730755"] <- "KRTAP2-3"
  data[[d]][data[[d]] == "MIR21"] <- "MIR21" 
  data[[d]][data[[d]] == "C5orf13"] <- "NREP"
  data[[d]][data[[d]] == "GPATC4"] <- "GPATCH4"
  data[[d]][data[[d]] == "SFRS10"] <- "TRA2B"
  data[[d]][data[[d]] == "DTX2P1-UPK3BP1-PMS2P11"] <- "DTX2P1-UPK3BP1-PMS2P11" 
  data[[d]][data[[d]] == "LINC00173"] <- "LINC00173" 
  data[[d]][data[[d]] == "HNRPDL"] <- "HNRNPDL"
  data[[d]][data[[d]] == "NOLA1"] <- "GAR1"
  data[[d]][data[[d]] == "HSPC111"] <- "NOP16"
  data[[d]][data[[d]] == "BAT1"] <- "DDX39B"
  data[[d]][data[[d]] == "WDSOF1"] <- "DCAF13"
  data[[d]][data[[d]] == "ZCD1"] <- "CISD1"
  data[[d]][data[[d]] == "C14orf156"] <- "SLIRP"
  data[[d]][data[[d]] == "HSP90AA2"] <- "HSP90AA2P" 
  data[[d]][data[[d]] == "IGKC"] <- "IGKC" 
  data[[d]][data[[d]] == "CHP"] <- "CHP1"
  data[[d]][data[[d]] == "TMCO7"] <- "TANGO6"
  data[[d]][data[[d]] == "VARS1"] <- "VARS"
  data[[d]][data[[d]] == "C4orf34"] <- "SMIM14"
  data[[d]][data[[d]] == "MOSC2"] <- "MARC2"
  data[[d]][data[[d]] == "KIAA0664"] <- "CLUH"
  data[[d]][data[[d]] == "CTPS"] <- "CTPS1"
  data[[d]][data[[d]] == "SNHG1"] <- "SNHG1" 
  data[[d]][data[[d]] == "SNHG13"] <- "DANCR" 
  data[[d]][data[[d]] == "SIK1B"] <- "SIK1B" 
  data[[d]][data[[d]] == "C22orf13"] <- "GUCD1"
  data[[d]][data[[d]] == "LOC645638"] <- "WFDC21P"
  data[[d]][data[[d]] == "HNRPLL"] <- "HNRNPLL"
  # Done
  data[[d]] <- unique(data[[d]])
  genes <- data[[d]]
  common <- intersect(all_genes, genes)
  missing <- genes[!(genes %in% common)]
  print(paste(d, " - Unmapped genes: ", paste(missing, collapse=","), sep=""))
}

max.length <- max(sapply(data, length))
df_data <- lapply(data, function(v) { c(v, rep(NA, max.length-length(v)))})
df_data <- t(do.call(rbind, df_data))
first_row <- data.frame(matrix(ncol=ncol(df_data),nrow=1))
colnames(first_row) <- colnames(df_data)
first_row[1,] <- "NA"
attach(ensembl_mapping)
df_data_ensembl <- matrix(with(ensembl_mapping, Ensembl[match(df_data[ ,1:ncol(df_data)], Symbol)]), nrow = nrow(df_data))
detach(ensembl_mapping)
df_data_ensembl <- as.data.frame(df_data_ensembl)
colnames(df_data_ensembl) <- colnames(df_data)
df_data_ensembl <- as.data.frame(rbind(first_row, df_data_ensembl))
#all_ensembl_IDs <- unique(unname(unlist(lapply(df_data_ensembl, function(v) { return(v[!is.na(v)][-1]) }))))
write.table(df_data_ensembl, file="myc_signature_genesets_ensembl.gmx", quote=FALSE, sep="\t", row.names=FALSE, na="")
df_data <- as.data.frame(rbind(first_row, df_data))
write.table(df_data, file="myc_signature_genesets.gmx", quote=FALSE, sep="\t", row.names=FALSE, na="")
