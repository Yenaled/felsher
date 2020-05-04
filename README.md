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

The following command will reproduce the normalized values in a sub-directory, named <i>processed_microarray</i>, created within the <i>microarray</i> directory.

<pre>Rscript microarray/normalize.r microarray/processed_microarray</pre>

## Downloading normalized data

The following command will download the normalized microarray values (from GEO) into a sub-directory, named <i>processed_microarray</i>, created within the <i>microarray</i> directory. (Use this if you only want to obtain the normalized data as deposited on GEO; keep in mind that the analysis on the paper doesn't use these values).

<pre>Rscript microarray/download.r microarray/processed_microarray</pre>

## Differential Gene Expression (Tissue)
### Probe-Level
The following command will perform differential gene expression analysis on normalized tissue microarray data deposited in <i>microarray/processed_microarray/tissue/</i> (from the previous step). The output of the analysis will be deposited in that same directory. Each tissue will have a set of files produced from the analysis. For example, for liver, we'll have the files liver<b>_diffexpr.csv</b> (containing the differential gene expression log fold changes and p-values), liver<b>_control_values.csv</b> (containing the probe intensities for the control samples), and liver<b>_myc_values.csv</b> (containing the probe intensities for the MYC-induced tumor samples).

<pre>Rscript microarray/diffexpr.r microarray/processed_microarray</pre>

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
</ul>
