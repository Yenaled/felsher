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

Pipeline readme here: https://github.com/Yenaled/felsher/blob/master/rnaseq/rna-seq/README.md

# TCGA<a name="tcga"></a>

## Obtaining TCGA data from UCSC Xena Toil

Let's download TCGA data from UCSC Xena Toil and put them into the tcga/ directory.

We obtain log2(x+1)-transformed DESeq2-normalized RSEM values (for between-samples comparisons) as follows:

<pre>wget --continue https://toil.xenahubs.net/download/TCGA-GTEx-TARGET-gene-exp-counts.deseq2-normalized.log2.gz
zcat < TCGA-GTEx-TARGET-gene-exp-counts.deseq2-normalized.log2.gz|awk -v OFS='\t' 'NR==1{for (i=1; i<=NF; i++) if (!($i ~ /TCGA/)) a[i]; 
     else printf "%s%s", $i, OFS; print ""; next}
     {for (i=1; i<=NF; i++) if (!(i in a)) printf "%s%s", $i, OFS; print "" }' | gzip > tcga/TCGA.rsem.deseq2.log2.tsv.gz</pre>

We obtain gene-level log2(x+0.001)-trransformed TPM values (for within-samples comparisons) as follows:

<pre>wget --continue https://toil.xenahubs.net/download/TcgaTargetGtex_rsem_gene_tpm.gz
zcat < TcgaTargetGtex_rsem_gene_tpm.gz|awk -v OFS='\t' 'NR==1{for (i=1; i<=NF; i++) if (!($i ~ /TCGA/)) a[i]; 
    else printf "%s%s", $i, OFS; print ""; next}
    {for (i=1; i<=NF; i++) if (!(i in a) || i==1) printf "%s%s", $i, OFS; print "" }' | gzip > tcga/TCGA.rsem.TPM.log2.tsv.gz</pre>
     
Alternately, the TCGA normalized RSEM and TPM files are posted here (which you can download and put into the tcga/ directory):
<ul>
  <li>https://github.com/Yenaled/felsher/releases/download/felsher/TCGA.rsem.deseq2.log2.tsv.gz</li>
  <li>https://github.com/Yenaled/felsher/releases/download/felsher/TCGA.rsem.TPM.log2.tsv.gz</li>
</ul>

## Obtaining TCGA batch identifier data (Tissue Source Site and Batch ID)

The final results of these steps are already in this github repository, but here is how those final results wereobtained:

1. Go to the GDC portal, navigate to the repository, add all the biospecimen files (bcr xml format) to your cart, then download the manifest file (save the manifest file as gdc_biospecimen.txt).
2. Using NCI's gdc-client software, run the following command to download the file into the path: tcga/biospecimen:

<pre>gdc-client download -m gdc_biospecimen.txt -d tcga/biospecimen</pre>

3. Then, run the following command to format the downloaded batch identifiers into tcga/tcga_batches.csv:

<pre>Rscript tcga/processgdc.r tcga/biospecimen tcga/tcga_batches.csv</pre>

4. Finally, download TCGA sample information into the tcga/ directory (more details on TCGA sample information are described below).

<pre>wget --continue https://gdc.cancer.gov/files/public/file/tcga_code_tables.zip
unzip tcga_code_tables.zip bcrBatchCode.tsv
unzip tcga_code_tables.zip tissueSourceSite.tsv
mv bcrBatchCode.tsv tcga/
mv tissueSourceSite.tsv tcga/</pre>

## Processing downloaded TCGA data

Now, we can process the TCGA data to only include samples of interest and to batch-correct batch ID and tissue source site information. Note that TCGA sample information (e.g. tissue source sites, tumor types, tumor vs. normal, etc.) is described at https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/ and how that information is encoded into TCGA barcodes can be found at https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/

Process TCGA data to get our final batch-corrected output in the tcga/ directory:

<pre>Rscript tcga/processTCGAbatches.r tcga/TCGA.rsem.deseq2.log2.tsv.gz tcga/bcrBatchCode.tsv tcga/tissueSourceSite.tsv tcga/TCGA.processed. tcga/tcga_batches.csv</pre>

The output files from the preceding command (i.e. the processed files associated with the batch-corrected data) will include:

<ul>
  <li>TCGA.processed.tumor_normal_lfc.tsv.gz - Log2 fold changes of paired tumor versus matched normal</li>
  <li>Batch correction:
    <ul>
      <li>TCGA.processed.info_corrected.csv - Information about cancer types, tumor/normal, batches, etc. for each sample</li>
      <li>TCGA.processed.batch_effects_view.pdf - PCA plots of the result of batch correction</li>
      <li>TCGA.processed.tumors_corrected.tsv.gz - Tumor sample expression data</li>
    </ul>
  </li>
  <li>Non-batch corrected:
    <ul>
      <li>TCGA.processed.info.csv - Information about cancer types, tumor/normal, batches, etc. for each sample</li>
      <li>TCGA.processed.tumors.tsv.gz - Tumor sample expression data</li>
    </ul>
  </li>
</ul>

Note: The smaller files are in this Github repository. However, the larger .tsv.gz files can be downloaded from:
<ul>
  <li>https://github.com/Yenaled/felsher/releases/download/felsher/TCGA.processed.tumor_normal_lfc.tsv.gz</li>
  <li>https://github.com/Yenaled/felsher/releases/download/felsher/TCGA.processed.tumors_corrected.tsv.gz</li>
  <li>https://github.com/Yenaled/felsher/releases/download/felsher/TCGA.processed.tumors.tsv.gz</li>
</ul>
  

## Obtaining TCGA clinical data


# Integrating Differential Gene Expression Analyses<a name="mouse_integration"></a>

The following command will integrate information from multiple differential gene expression analyses: microarrays of liver, kidney, and lung cancers as well as RNA-seq of lymphomas. We control FDR at 0.05 and set the fold change threshold to 2 (i.e. a log2 Fold Change of 1). The output will be placed in output/mouse_de/.

<pre>Rscript overlap.r 2 0.05 output/mouse_de/ liver_myc,kidney_myc,lung_myc,lung_mycras,lung_ras,tall_myc,eumyc_myc "HCC,RCC,LAC (MYC),LAC (MYC+KRAS),LAC (KRAS),T-ALL,BCL" microarray/processed_microarray/tissue/liver_diffexpr_gene.csv,microarray/processed_microarray/tissue/kidney_diffexpr_gene.csv,microarray/processed_microarray/tissue/lung_diffexpr_gene.csv,microarray/processed_microarray/tissue/lung_mycras_diffexpr_gene.csv,microarray/processed_microarray/tissue/lung_ras_diffexpr_gene.csv,rnaseq/tall_rnaseq/kallisto_aligned/tallmycon_sleuth_results_genes.csv,rnaseq/eumyc_rnaseq/kallisto_aligned/eumycon_sleuth_results_genes.csv microarray/annotation/MouseWG-6v2.csv rnaseq/annotation/gencode_GRCm38_vM15.csv</pre>  

# Downstream Analyses of Mouse Differential Gene Expression<a name="mouse_downstream_analyses"></a>

## Mouse tissue enrichment analysis<a name="mouse_tissue_enrichment"></a>

Run the following script for mouse tissue enrichment analysis. Output folder: output/enrichr_mouse_tissue
<pre>Rscript analyze_mouse_tissue_enrichment.r</pre>

## Mouse DE genes in TCGA<a name="mouse_tcga"></a>

Run the following script to get the expression log2 fold changes of mouse DE genes in TCGA human cancer samples (comparing human tumor with matched normal). Output folder: output/mouse_tcga

<pre>Rscript analyze_mouse_genes_tcga.r</pre>

## Mouse DE genes GO analysis<a name="mouse_de_go"></a>

Run the following script to perform GSEA of GO terms for mouse DE genes. For the top enriched GO terms, the script will also perform GO overrepresentation analysis for relevant mouse gene atlas tissue gene sets. This is to assess whether enriched GO terms in mouse DE genes have tissue-lineage dependence. Output folder: output/enrichr_mouse_go

The script has two parts. The first part as follows and will generate the GSEA results as well as a list of top GO terms:

<pre>Rscript analyze_mouse_go_enrichment.r</pre>

The GSEA results will be outputted in go_gsea.xlsx and the top GO terms will be outputted in topterms.txt in the output folder. We then manually curate the top GO terms and assign them to broader categories for ease of display. The category assignment should be in a .csv file with each column name being a category name and with each GO term being assigned to a specific column. One column should be named "Other" and this contains GO terms that couldn't be assigned to a particular category (and therefore won't be displayed on the heatmap). Finally, we run the second part of the script -- supplying the .xlsx GSEA results (generated from the first part) and the .csv category assignment file (that we manually created):

<pre>Rscript analyze_mouse_go_enrichment.r output/enrichr_mouse_go/go_gsea.xlsx data/go_term_collapse.csv</pre>

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
      <li>C/D:
        <ul>
          <li>Graphpad Prism file: prism/mouse_tissue_enrichment.pzfx</li>
          <li>Raw output: output/enrichr_mouse_tissue/</li>
          <li><a href="#mouse_tissue_enrichment">Link to method</a></li>
        </ul>
      </li>
      <li>E:
        <ul>
          <li>Graphpad Prism file: prism/mouse_tcga.pzfx</li>
          <li>Raw output: output/mouse_tcga/</li>
          <li><a href="#mouse_tcga">Link to method</a></li>
        </ul>
      </li>
    </ul>
  </li>
  <li>
    <ul>
      <li>Under construction</li>
    </ul>
  </li>
  <li>
    <ul>
      <li>A:
        <ul>
          <li>PDF file: output/enrichr_mouse_go/heatmap_go.pdf</li>
          <li>Raw output:
            <ul>
              <li>output/enrichr_mouse_go/go_gsea.xlsx</li>
              <li>output/enrichr_mouse_go/tissue_go.xlsx</li>
              <li>output/enrichr_mouse_go/heatmap_go_terms.csv</li>
            </ul>
          </li>
          <li><a href="#mouse_de_go">Link to method</a></li>
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
  <li>
    <ul>
      <li>output/enrichr_mouse_go/heatmap_go_terms.csv
        <ul>
          <li><a href="#mouse_de_go">Link to method</a></li>
        </ul>
      </li>
    </ul>
  </li>
  <li>
    <ul>
      <li>output/enrichr_mouse_go/go_gsea.xlsx
        <ul>
          <li><a href="#mouse_de_go">Link to method</a></li>
        </ul>
      </li>
    </ul>
  </li>
  <li>
    <ul>
      <li>output/enrichr_mouse_go/tissue_go.xlsx
        <ul>
          <li><a href="#mouse_de_go">Link to method</a></li>
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
  <li>ComplexHeatmap 2.2.0</li>
  <li>survival 3.1-12</li>
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
