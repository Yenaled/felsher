# Table of contents
1. [Microarray](#microarray)
2. [RNA-Seq](#rnaseq)
3. [TCGA](#tcga)
4. [CCLE](#ccle)
5. [Machine Learning](#ml)
6. [ChIP-Seq](#chipseq)
7. [Gene Symbol Mapping Tool](#mapping)
8. [Appendix](#appendix)

# Microarray<a name="microarray"></a>

The microarray is deposited at GEO accession number: <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE143254">GSE143254</a> with following two subseries:
<ul>
  <li>Tissue: <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE143253">GSE143253</a></li>
  <li>Cell lines: <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE143250">GSE143250</a></li>
</ul>

## Normalizing Raw Data

The following command will obtain the raw and normalized microarray values and put them into a sub-directory, named <i>processed_microarray</i>, created within the <i>microarray</i> directory using the Mouse-WG6 annotation file supplied (to map probe IDs to gene IDs/symbols). The normalized values deposited on GEO will be reproduced and, separately, the normalized values used in the paper (where poorly-detected probes were filtered) will also be produced.

<pre>Rscript microarray/normalize.r microarray/processed_microarray microarray/annotation/MouseWG-6v2.csv</pre>

## Downloading normalized data

The following command will download the normalized microarray values (from GEO) into a sub-directory, named <i>processed_microarray</i>, created within the <i>microarray</i> directory. (Use this if you only want to obtain the normalized data as deposited on GEO; keep in mind that the analysis on the paper doesn't use these values).

<pre>Rscript microarray/download.r microarray/processed_microarray</pre>

## Differential Gene Expression (Tissue)

The following command will perform differential gene expression analysis on normalized tissue microarray data deposited in <i>microarray/processed_microarray/tissue/</i> (described previously). The output of the analysis will be deposited in that same directory. Each tissue will have a set of files produced from the analysis. For example, for liver, we'll have a set of files prefixed with <i>liver</i> containing differential gene expression analysis and expression values. Both the differential gene expression analysis and expression values are done at probe-level and are aggregated at gene-level (for differential expression, Fisher's method is used to combine p-values and the probe with the highest absolute value of log2 fold change is used as effect size; for expression values, the mean of the probe intensities is used for gene-level aggregation).

<pre>Rscript microarray/diffexpr.r microarray/processed_microarray microarray/annotation/MouseWG-6v2.csv</pre>

## Rank Product Analysis (Cell Lines)

# RNA-Seq<a name="rnaseq"></a>

# TCGA<a name="tcga"></a>

# CCLE<a name="ccle"></a>

# Machine Learning<a name="ml"></a>

# ChIP-Seq<a name="chipseq"></a>

# Gene Symbol Mapping Tool<a name="mapping"></a>

# Appendix<a name="appendix"></a>

## Software Versions
<ul>
  <li>R version 3.6.1</li>
  <li>GEOquery 2.54.0</li>
  <li>lumi 2.38.0</li>
  <li>limma 3.42.0</li>  
  <li>aggregation 1.0.1</li>
  <li>UpSetR 1.4.0</li>
</ul>

## Loading files in Excel

When loading .csv files or files containing tab-separated values in Excel, use File -> Import to import the files into Excel and set "Column data format" for all columns containing gene names/symbols to "Text" rather than "General". This will allow you to avoid issues with gene names being reformatted (e.g. March7 becoming 7-Mar).
