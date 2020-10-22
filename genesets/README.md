**myc_signature_genesets_ensembl.gmx** and **myc_signature_genesets.gmx** are the files containing gene sets from previous MYC studies.

**aggregated_genesets_myc_up.RData** contains an object named 'data' that can be loaded in R. This object contains the raw downloaded gene sets. The **process.r** script processes the gene sets from the RData file to clean up gene names and map gene names to ENSEMBL IDs, generating the .gmx files.
