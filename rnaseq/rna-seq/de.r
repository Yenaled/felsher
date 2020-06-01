# Differential gene expression analysis using DESeq2

# Load R packages
pkgs <- c("DESeq2", "BiocParallel")
invisible(lapply(pkgs, function(x) suppressWarnings(suppressMessages(library(x, character.only=TRUE)))))

# Process command-line arguments
args <- commandArgs(TRUE)
numcores <- as.numeric(args[1]) # Number of CPU cores
input.path <- args[2] # Path to input files
suffix <- args[3] # Suffix of input files
alpha_value <- as.numeric(args[4]) # alpha value threshold
lfc_threshold <- as.numeric(args[5]) # LFC threshold
prefix <- args[6] # output file prefix
samples <- unlist(strsplit(args[7],",")) # Sample names separated by commas
variables <- args[8:(length(args)-3)] # All remaining args (except the last 3) contain explanatory variables
design_string <- args[length(args)-2] # The design formula
gtf_file <- args[length(args)-1] # GTF file for calculating FPKMs
coef_comparison <- args[length(args)] # Number specifying which variable in the design to compare

# Register cores for parallel processing
register(MulticoreParam(numcores))

# Prepare data matrix by going through samples
data <- NULL
identifiers <- NULL # For tximport
for (input in samples) {
    file <- paste(input.path,"/",input,suffix,sep="")
    d <- NULL
    if (suffix == "FeatureCounts.txt") { # Using featureCounts results
	d <- read.table(file, header=TRUE, sep="\t", quote="\"", row.names=1, check.names=FALSE, stringsAsFactors=FALSE)
	d <- d[,ncol(d),drop=FALSE]
    } else if (suffix == "kallisto") {
        require(rhdf5)
        file <- paste(input.path,"/",input,"/abundance.h5",sep="")
        if (is.null(identifiers)) {
            identifiers <- as.character(rhdf5::h5read(file, "aux/ids"))
        } else {
            identifiers2 <- as.character(rhdf5::h5read(file, "aux/ids"))
            if (!(setequal(identifiers, identifiers2) && length(identifiers) == length(identifiers2))) {
                print("Error: Transcript identifiers in supplied files are not identical")
                quit()
            }
        }
        if (is.null(data)) {
            data <- list()
        }
        data[[input]] <- file
        next
    } else {
	d <- read.table(file,row.names=1,header=FALSE,stringsAsFactors=FALSE)
	if (ncol(d) == 3) { # If using STAR quantMode (STAR output gives 4 columns [first column = our row names]) 
		d <- d[,3,drop=FALSE] # Use fourth (technically third) column
		d <- d[!(rownames(d) %in% c("N_unmapped","N_multimapping","N_noFeature","N_ambiguous")),,drop=FALSE] # Get rid of these rows
    	} else {
    		d <- d[,1,drop=FALSE] # Only care about the first column
    	}
    }
    if (is.null(data)) {
        data <- d
    } else {
        data <- cbind(data, d)
    }
}
if (suffix == "kallisto") { # Do tximport
    require(tximport)
    require(GenomicFeatures)
    data <- unlist(data)
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
    tx2gene <- mapping[,c("target_id","gene_name")]
    colnames(tx2gene) <- c("TXNAME", "GENEID")
    txi.kallisto <- tximport(data, type = "kallisto", tx2gene=tx2gene)
    data <- txi.kallisto$counts
    mapping <- NULL
} else {
    data <- as.matrix(data)
    colnames(data) <- samples
}
# Prepare data.info dataframe which consists of the explanatory variables
data.info <- data.frame(row.names=colnames(data))
for (variable in variables) {
        v = unlist(strsplit(variable,":")) # Parse the string which is "variable_name:values"
        variable_name = v[1]
        variable_values = unlist(strsplit(v[2],",")) # Values are comma-separated and correspond to samples (in order)
	data.info[,variable_name] <- variable_values
}

# Prepare GTF annotation file
ebg <- NULL # Exons by gene
if (suffix != "kallisto") {
    mapping <- NULL # Mapping from ENSEMBL to Gene Name
    if (toupper(gtf_file) != "NULL") {
            if (require(GenomicFeatures)) {
                    txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")
                    ebg <- exonsBy(txdb, by = "gene")
		    mapping <- read.table(gtf_file, header=FALSE, stringsAsFactors=FALSE, sep="\t")
		    mapping <- mapping[,9]
		    mapping <- lapply(mapping, function(x) {
		        a <- trimws(unlist(strsplit(x, ";")))
		        return(c("gene_id"=substring(a[grepl("gene_id", a)], 9), "gene_name"=gsub("\"", "", substring(a[grepl("gene_name", a)], 11))))
		    })
		    mapping <- as.data.frame(do.call(rbind, mapping))
		    mapping <- mapping[!duplicated(mapping$gene_id),]
            }
    }
}

print(model.matrix(eval(parse(text=design_string)), data.info))

if (suffix == "kallisto") {
    dds <- DESeqDataSetFromTximport(txi.kallisto, colData = data.info, design=eval(parse(text=design_string)))
} else {
    dds <- DESeqDataSetFromMatrix(countData = data, colData = data.info, design=eval(parse(text=design_string)))
}
dds <- DESeq(dds, parallel=TRUE)
print(resultsNames(dds))
if (toupper(coef_comparison) == "NULL") {
	variable_coef <- length(resultsNames(dds))
} else {
	variable_coef <- as.numeric(variable_coef)
}
print(paste("Doing differential gene expression on variable ",variable_coef,": ", resultsNames(dds)[variable_coef], sep=""))
res <- results(dds, parallel=TRUE, alpha=alpha_value, lfcThreshold=lfc_threshold)
res <- lfcShrink(dds=dds, coef=length(resultsNames(dds)), res=res, parallel=TRUE)
res <- res[order(res$pvalue),]

# Function to prepare object for being written out to a file

prepForOutput <- function(x, mapping) {
	x <- as.data.frame(x)
	if (!is.null(mapping)) {
		x <- merge(mapping, x, by.x="gene_id", by.y="row.names", all.x=FALSE, all.y=TRUE)
		rownames(x) <- x$gene_id
		x$gene_id <- NULL
		rownames(x) <- gsub("\\..*", "", rownames(x)) # Get rid of periods in ENSEMBL IDs
	}
	return(x)
}

# Time to output files and do some transformations

write.csv(prepForOutput(res, mapping), paste(input.path, "/", prefix, "_deseq2_results.csv", sep=""), quote=FALSE)

transformed <- rlog(dds, blind = FALSE)
write.csv(prepForOutput(assay(transformed), mapping), file=paste(input.path, "/", prefix, "_deseq2_normalized_counts.csv", sep=""), quote=FALSE)
if ("Batch" %in% colnames(data.info)) {
	if (require("sva")) { # Apply ComBat batch correction
		assay(transformed) <- ComBat(assay(transformed), transformed$Batch)
		write.csv(prepForOutput(assay(transformed), mapping), file=paste(input.path, "/", prefix, "_deseq2_normalized_combat_counts.csv", sep=""), quote=FALSE)
	}
}
write.csv(prepForOutput(counts(dds,normalized=FALSE), mapping), paste(input.path, "/", prefix, "_deseq2_raw_counts.csv", sep=""), quote=FALSE)

if (!is.null(ebg)) {
	dds <- dds[rownames(dds) %in% names(ebg),]
	# Remove: rowRanges(dds) <- GRangesList(ebg)
	gene_lengths <- as.data.frame(sum(width(reduce(ebg))))
	gene_lengths$Geneid <- rownames(gene_lengths)
	gene_lengths <- gene_lengths[match(rownames(mcols(dds, use.names=TRUE)), gene_lengths$Geneid),]
	mcols(dds)$basepairs <- gene_lengths[,1]
	write.csv(prepForOutput(fpkm(dds,robust=FALSE), mapping), paste(input.path, "/", prefix, "_deseq2_fpkm.csv", sep=""), quote=FALSE)
}
