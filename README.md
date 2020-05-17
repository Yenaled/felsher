# Table of contents
1. [Microarray](#microarray)
2. [RNA-Seq](#rnaseq)
3. [TCGA](#tcga)
4. [Integrating Differential Gene Expression Analyses](#mouse_integration)
5. [Downstream Analyses of Mouse Differential Gene Expression](#mouse_downstream_analyses)
6. [CCLE](#ccle)
7. [Machine Learning](#ml)
8. [ChIP-Seq](#chipseq)
9. [Figures](#figures)
10. [Appendix](#appendix)

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

# TCGA<a name="tcga"></a>

## Obtaining TCGA data from GDC

Here, we will go through an example of how to download data for the TCGA-LIHC study:
1. Go to the GDC portal, navigate to the TCGA-LIHC repository, add all the biospecimen files (bcr xml format) and the RNA-seq files (HTSeq - Counts) to your cart, then download the manifest file (save the manifest file as gdc_manifest_lihc.txt).
2. Using NCI's gdc-client software, run the following command to download the file into the path: tcga/data/lihc:

<pre>gdc-client download -m gdc_manifest_lihc.txt -d tcga/data/lihc</pre>

3. Then, run the following command to organize the downloaded data (the organized data will be stored in tcga/organized/lihc):

<pre>tcga/organize_gdcdata.sh lihc</pre>

4. Then, run the following command to further process the data (which will be outputted into tcga/processed/lihc):

<pre>tcga/processgdc.sh lihc</pre>

# Integrating Differential Gene Expression Analyses<a name="mouse_integration"></a>

The following command will integrate information from multiple differential gene expression analyses: microarrays of liver, kidney, and lung cancers as well as RNA-seq of lymphomas. We control FDR at 0.05 and set the fold change threshold to 2 (i.e. a log2 Fold Change of 1). The output will be placed in output/mouse_de/.

<pre>Rscript overlap.r 2 0.05 output/mouse_de/ liver_myc,kidney_myc,lung_myc,lung_mycras,lung_ras,tall_myc,eumyc_myc "HCC,RCC,LAC (MYC),LAC (MYC+KRAS),LAC (KRAS),T-ALL,BCL" microarray/processed_microarray/tissue/liver_diffexpr_gene.csv,microarray/processed_microarray/tissue/kidney_diffexpr_gene.csv,microarray/processed_microarray/tissue/lung_diffexpr_gene.csv,microarray/processed_microarray/tissue/lung_mycras_diffexpr_gene.csv,microarray/processed_microarray/tissue/lung_ras_diffexpr_gene.csv,rnaseq/tall_rnaseq/kallisto_aligned/tallmycon_sleuth_results_genes.csv,rnaseq/eumyc_rnaseq/kallisto_aligned/eumycon_sleuth_results_genes.csv microarray/annotation/MouseWG-6v2.csv rnaseq/annotation/gencode_GRCm38_vM15.csv</pre>  

# Downstream Analyses of Mouse Differential Gene Expression<a name="mouse_downstream_analyses"></a>

## Mouse tissue enrichment analysis<a name="mouse_tissue_enrichment"></a>

Run the following script for mouse tissue enrichment analysis. Output folder: output/enrichr_mouse_tissue
<pre>Rscript analyze_mouse_tissue_enrichment.r</pre>

# CCLE<a name="ccle"></a>

# Machine Learning<a name="ml"></a>

# ChIP-Seq<a name="chipseq"></a>

# Figures<a name="figures"></a>

Figures in paper were generated as follows:

<ol>
  <li>
    <ul>
      <li>A: Self-created figure</li>
      <li>B:
        <ul>
          <li>PDF file: output/mouse_de/venn.pdf</li>
          <li>Raw output: output/mouse_de/</li>
          <li><a href="#mouse_integration">Link to method</a></li>
        </ul>
      </li>
      <li>C: Self-created figure</li>
      <li>D:
        <ul>
          <li>Graphpad Prism file: prism/mouse_tissue_enrichment.pzfx</li>
          <li>Raw output: output/enrichr_mouse_tissue/</li>
          <li><a href="#mouse_tissue_enrichment">Link to method</a></li>
        </ul>
      </li>
    </ul>
  </li>
</ol>

Supplementary figures in paper were generated as follows:

<ol>
  <li>
    <ul>
      <li>Self-created figure</li>
    </ul>
  </li>
  <li>
    <ul>
      <li>output/mouse_de/upset_up.pdf
        output/mouse_de/upset_down.pdf
        <ul>
          <li>Raw output: output/mouse_de/</li>
          <li><a href="#mouse_integration">Link to method</a></li>
        </ul>
      </li>
    </ul>
  </li>
</ol>

Supplementary tables in paper were generated as follows:

<ol>
  <li>
    <ul>
      <li>output/mouse_de/de_table.xlsx
        <ul>
          <li><a href="#mouse_integration">Link to method</a></li>
        </ul>
      </li>
    </ul>
  </li>
  <li>
    <ul>
      <li>output/enrichr_mouse_tissue/tissue_enrichment_table.xlsx
        <ul>
          <li><a href="#mouse_tissue_enrichment">Link to method</a></li>
        </ul>
      </li>
    </ul>
  </li>
</ol>

# Appendix<a name="appendix"></a>

## Software Versions
<ul>
  <li>Graphpad Prism 8.4.1</li>
  <li>gdc-client 1.5.0</li>
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
  <li>XML 3.98-1.20</li>
  <li>sva 3.34.0</li>
  <li>BiocParallel 1.20.0</li>
</ul>

## Gene identifiers ortholog mapping

To perform mapping between gene identifiers (e.g. NCBI ID, ENSEMBL, Symbol, etc.) or orthologs for annotation purposes, please see the files distributed here:

https://ftp.ncbi.nlm.nih.gov/gene/DATA/

Here are some files of note:
<ul>
  <li>https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz - Mapping between NCBI Gene ID and ENSEMBL ID</li>
  <li>https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz - Mapping between orthologs by NCBI IDs (note: human taxonomy ID: 9606; mouse taxonomy ID: 10090)</li>
  <li>https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz - Mapping between human gene symbols and NCBI IDs, among other things</li>
    <li>https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz - Mapping between mouse gene symbols and NCBI IDs, among other things</li>
</ul>

## Loading files in Excel

When loading .csv files or files containing tab-separated values in Excel, use File -> Import to import the files into Excel and set "Column data format" for all columns containing gene names/symbols to "Text" rather than "General". This will allow you to avoid issues with gene names being reformatted (e.g. March7 becoming 7-Mar).
