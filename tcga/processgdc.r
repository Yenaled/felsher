###############################################################
### File: processgdc.r
### Description: Gets gene expression values
###              from organized GDC gene data,
###              among other data.
### Written by Delaney Sullivan
###############################################################

# Load R packages
pkgs <- c("XML")
invisible(lapply(pkgs, function(x) suppressWarnings(suppressMessages(library(x, character.only=TRUE)))))

# Set up some variables
args <- commandArgs(TRUE) # Extract command-line arguments
project.name <- args[1]
current_directory <- args[2]
project.control <- args[3]
project.experimental <- args[4]
project.output <- args[5]
project.path <- paste(current_directory,"/organized/",project.name,sep="")

# Iterate through folders containing the data:

organized.subfolders <- list.dirs(path=project.path, recursive=FALSE)

# Retrieves RNA-seq data (of type datatype) for a patient within group (group being Tumor, Normal, etc.).
process.data <- function(patient, datatype, group, log2=FALSE, has_header=FALSE) {
    filename <- paste(patient,"/",group,"/", datatype, ".txt",sep="")
    if (file.exists(filename)) {
        # Process the data into a dataframe
        data <- read.table(filename, sep="\t", header=has_header)
        rownames(data) <- data[,1]
        data[,1] <- NULL
        data <- data[,1,drop=FALSE] # If there's multiple columns, just take the first one
        colnames(data) <- c(basename(patient))
        if (log2) {
            data <- log2(data+1) # log2(x+1) transformation
        }
        # Transpose rows and cols for now, and return
        return(t(data))
    } else {
        return(NULL)
    }
}

# Processes biospecimen data to associate patients with batch numbers
process.batchnums <- function(patient) {
    filename <- paste(patient,"/", "biospecimen.xml",sep="")
    
    if (file.exists(filename)) {
        data <- xmlParse(filename)
        data <- xmlToList(data)
        batch.number <- unlist(data[["admin"]][["batch_number"]][["text"]])
        return(data.frame(batch.number, row.names = c(basename(patient))))
    } else {
        return(NULL)
    }
}

if (dir.exists(project.output)) {
    unlink(project.output, recursive=TRUE) # Remove output directory if it already exists
}
dir.create(project.output, FALSE) # Create directory for output

# Write basic information about control and experimental groups to file
writeLines(paste("Control: \"", project.control, "\"\nExperimental: \"", project.experimental, "\"", sep=""), paste(project.output, "/info.txt", sep=""))

# Write counts data from experimental group to file
result <- do.call(rbind, lapply(organized.subfolders, process.data, datatype="counts", group=project.experimental))
if (!all(is.null(result))) {
    write.csv(t(result), file = paste(project.output, "/counts_experimental.csv", sep=""), quote = FALSE)
}

# Write counts data from control group to file
result <- do.call(rbind, lapply(organized.subfolders, process.data, datatype="counts", group=project.control))
if (!all(is.null(result))) {
    write.csv(t(result), file = paste(project.output, "/counts_control.csv", sep=""), quote = FALSE)
}

# Write batch numbers and their associated patients to file
result <- do.call(rbind, lapply(organized.subfolders, process.batchnums))
if (!all(is.null(result))) {
    result <- t(result)
    write.table(data.frame(Patient = colnames(result), Batch = result["batch.number",]), file=paste(project.output, "/batch.txt", sep=""), row.names = FALSE, quote = FALSE)
}
