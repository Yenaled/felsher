# Table of contents
1. [Microarray](#microarray)
2. [RNA-Seq](#rnaseq)
3. [Integrating Differential Gene Expression Analyses](#mouse_integration)
4. [Downstream Analyses of Mouse Differential Gene Expression Data](#mouse_downstream_analyses)
5. [TCGA](#tcga)
6. [CCLE](#ccle)
7. [Machine Learning](#ml)
8. [ChIP-Seq](#chipseq)
9. [Gene Symbol Mapping Tool](#mapping)
10. [Figures](#figures)
11. [Appendix](#appendix)

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

# RNA-Seq<a name="rnaseq"></a>

# Integrating Differential Gene Expression Analyses<a name="mouse_integration"></a>

The following command will integrate information from multiple differential gene expression analyses: microarrays of liver, kidney, and lung cancers as well as RNA-seq of lymphomas. We control FDR at 0.05 and set the fold change threshold to 2 (i.e. a log2 Fold Change of 1). The output will be placed in output/mouse_de/.

<pre>Rscript overlap.r 2 0.05 output/mouse_de/ liver_myc,kidney_myc,lung_myc,lung_mycras,lung_ras,tall_myc,eumyc_myc "HCC,RCC,LAC (MYC),LAC (MYC+KRAS),LAC (KRAS),T-ALL,BCL" microarray/processed_microarray/tissue/liver_diffexpr_gene.csv,microarray/processed_microarray/tissue/kidney_diffexpr_gene.csv,microarray/processed_microarray/tissue/lung_diffexpr_gene.csv,microarray/processed_microarray/tissue/lung_mycras_diffexpr_gene.csv,microarray/processed_microarray/tissue/lung_ras_diffexpr_gene.csv,rnaseq/tall_rnaseq/kallisto_aligned/tallmycon_sleuth_results_genes.csv,rnaseq/eumyc_rnaseq/kallisto_aligned/eumycon_sleuth_results_genes.csv microarray/annotation/MouseWG-6v2.csv rnaseq/annotation/gencode_GRCm38_vM15.csv</pre>  

# Downstream Analyses of Mouse Differential Gene Expression Data<a name="mouse_downstream_analyses"></a>

## Mouse tissue enrichment analysis<a name="mouse_tissue_enrichment"></a>

Run the following script for mouse tissue enrichment analysis. Output folder: output/enrichr_mouse_tissue
<pre>Rscript analyze_mouse_tissue_enrichment.r</pre>

# TCGA<a name="tcga"></a>

# CCLE<a name="ccle"></a>

# Machine Learning<a name="ml"></a>

# ChIP-Seq<a name="chipseq"></a>

# Gene Symbol Mapping Tool<a name="mapping"></a>

# Figures<a name="figures"></a>

Figures in paper were generated as follows:
<ol>
  <li>
    <ol style="list-style-type: upper-alpha" type="A">
      <li>Self-created figure</li>
      <li>
        <ul>
          <li>PDF file: output/mouse_de/venn.pdf</li>
          <li>Raw data: output/mouse_de/</li>
          <li><a href="mouse_integration">Link to method</a></li>
        </ul>
      </li>
      <li>Self-created figure</li>
      <li>
        <ul>
          <li>Graphpad Prism file: prism/mouse_tissue_enrichment.pzfx</li>
          <li>Raw data: output/enrichr_mouse_tissue/</li>
          <li><a href="mouse_tissue_enrichment">Link to method</a></li>
        </ul>
      </li>
    </ol>
  </li>
</ol>
1. Test
a. test
b. test

# Appendix<a name="appendix"></a>

## Software Versions
<ul>
  <li>Graphpad Prism 8.4.1</li>
  <li>R version 3.6.1</li>
  <li>GEOquery 2.54.0</li>
  <li>lumi 2.38.0</li>
  <li>limma 3.42.0</li>  
  <li>aggregation 1.0.1</li>
  <li>UpSetR 1.4.0</li>
  <li>VennDiagram 1.6.20</li>
  <li>openxlsx 4.1.5</li>
  <li>kallisto 0.44.0</li>
  <li>sleuth 0.30.0</li>
  <li>DESeq2 1.22.2</li>
  <li>GenomicFeatures 1.38.0</li>
  <li>fgsea 1.12.0</li>
  <li>enrichR 2.1</li>
</ul>

## Loading files in Excel

When loading .csv files or files containing tab-separated values in Excel, use File -> Import to import the files into Excel and set "Column data format" for all columns containing gene names/symbols to "Text" rather than "General". This will allow you to avoid issues with gene names being reformatted (e.g. March7 becoming 7-Mar).
