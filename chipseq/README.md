# chip-seq
Pipeline for ChIP-Seq analysis. We will be using the ENCODE pipeline.
<details><summary>Click here for instructions on getting the pipeline (for Stanford's Sherlock cluster)</summary>
<p>Instructions for getting the pipeline onto Sherlock: https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/master/docs/tutorial_sherlock.md</p>
<p>We will be using Singularity so we won't be installing conda (i.e. disregard the step that specifies "For Conda users").</p>
<p>Within the new chip-seq-pipeline2 directory, edit the file <b>workflow_opts/sherlock.json</b> to remove /oak/stanford from singularity_bindpath (note: this is not necessary if the lab has paid for oak storage). Additionally, you might also want to consider having a dedicated directory within the lab group folder as a place to store all the sequencing data (e.g. /home/groups/dfelsher/delaneysequencing). In which case, you should create that directory (using mkdir) and then add that directory in singularity_bindpath. All in all, the resulting workflow_opts/sherlock.json file should look something like the following (note: v1.1.5 may be different depending on your version number):</p>
<pre>{
    "default_runtime_attributes" : {
        "slurm_partition" : "normal",
        "singularity_container" : "/home/groups/cherry/encode/pipeline_singularity_images/chip-seq-pipeline-v1.1.5.simg",
        "singularity_bindpath" : "/scratch,/lscratch,/home/groups/dfelsher/delaneysequencing,/home/groups/cherry/encode"
    }
}</pre>
</details>

## Versions

<ul>
    <li>ENCODE-DCC chip-seq-pipeline2 v1.1.5</li>
    <li>deeptools 3.1.3</li>
    <li>HOMER v4.10</li>
</ul>

# Step 1: Run the ENCODE pipeline

(Skip to step 2 to start from the processed files rather than raw data)

## Part a: Setup your raw data

Follow the instructions here for preliminary quality control and processing of your raw fastQ sequencing read files: https://github.com/Yenaled/felsher/tree/master/rnaseq/prep_fastq

## Part b: Create the input configuration file

Read <i>carefully</i> the input file specification: https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/master/docs/input.md and alter the options as needed. You can add additional replicates (e.g. fastqs_rep3_R1, fastqs_rep4_R1, etc.) if you have them available (be sure to add the corresponding _R2 if you have paired-end reads and specify the corresponding control FastQs as well).
<p>IMPORTANT:</p>
<ul>
    <li>The fastq files must be gzipped (.gz) otherwise the pipeline might fail!! </li>
    <li>The input json file must should not have a comma (,) after the very last element is specified (i.e. on the line right before the curly brace } at the end of the file) so alter the input file to make sure that final comma isn't there.</li>
 </ul>
 
Note: Edit the input file to use macs2 as the chip.peak_caller regardless of whether you are doing transcription factor ChIP and histone ChIP. (SPP is too slow as a peak caller).
 
## Part c: Execute the pipeline
 
See further instructions at https://github.com/ENCODE-DCC/chip-seq-pipeline2

# Step 2: Collect the files

Put all the processed files into the data/ directory. If working from the pipeline above, you'll need to collect the .fc.signal.bigwig fold change signal track files. The path to it might look something like call-macs2_pooled/execution or, if there are no replicates, call-macs2/shard-0/execution/. Alternately, the .fc.signal.bigwig files (the bigwig files containing the fold change signal tracks) are listed below:

<ul>
	<li><a href="https://github.com/Yenaled/felsher/releases/download/felsher/eumyc_h3k27ac_C.fc.signal.bigwig">eumyc_h3k27ac_C.fc.signal.bigwig</a></li>
	<li><a href="https://github.com/Yenaled/felsher/releases/download/felsher/eumyc_h3k27ac_T.fc.signal.bigwig">eumyc_h3k27ac_T.fc.signal.bigwig</a></li>
	<li><a href="https://github.com/Yenaled/felsher/releases/download/felsher/eumyc_h3k4me3_C.fc.signal.bigwig">eumyc_h3k4me3_C.fc.signal.bigwig</a></li>
	<li><a href="https://github.com/Yenaled/felsher/releases/download/felsher/eumyc_h3k4me3_T.fc.signal.bigwig">eumyc_h3k4me3_T.fc.signal.bigwig</a></li>
	<li><a href="https://github.com/Yenaled/felsher/releases/download/felsher/hcc_h3k27ac_C.fc.signal.bigwig">hcc_h3k27ac_C.fc.signal.bigwig</a></li>
	<li><a href="https://github.com/Yenaled/felsher/releases/download/felsher/hcc_h3k27ac_T.fc.signal.bigwig">hcc_h3k27ac_T.fc.signal.bigwig</a></li>
	<li><a href="https://github.com/Yenaled/felsher/releases/download/felsher/hcc_h3k4me3_C.fc.signal.bigwig">hcc_h3k4me3_C.fc.signal.bigwig</a></li>
	<li><a href="https://github.com/Yenaled/felsher/releases/download/felsher/hcc_h3k4me3_T.fc.signal.bigwig">hcc_h3k4me3_T.fc.signal.bigwig</a></li>
</ul>
	

# Step 3: Analysis with deeptools

<p>Install <b>deeptools</b>: https://deeptools.readthedocs.io/en/develop/content/installation.html</p>

The reference genome (e.g. for mm10 mouse) in .bed format is in the data/ folder and was obtained via:
<pre>wget -qO- http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/refGene.txt.gz | gunzip -c - | awk 'BEGIN{ OFS="\t" }{ print $3, $5, $6, $13, ".", $4  }' - > data/refGene.bed</pre>

We can create the log2 ratio (tumor vs. control) bigwig files in the data/ directory via the following:

<pre>
bigwigCompare --bigwig2 data/eumyc_h3k27ac_C.fc.signal.bigwig --bigwig1 data/eumyc_h3k27ac_T.fc.signal.bigwig --outFileName data/log2ratio_eumyc_h3k27ac.bigwig
bigwigCompare --bigwig2 data/eumyc_h3k4me3_C.fc.signal.bigwig --bigwig1 data/eumyc_h3k4me3_T.fc.signal.bigwig --outFileName data/log2ratio_eumyc_h3k4me3.bigwig
bigwigCompare --bigwig2 data/hcc_h3k27ac_C.fc.signal.bigwig --bigwig1 data/hcc_h3k27ac_T.fc.signal.bigwig --outFileName data/log2ratio_hcc_h3k27ac.bigwig
bigwigCompare --bigwig2 data/hcc_h3k4me3_C.fc.signal.bigwig --bigwig1 data/hcc_h3k4me3_T.fc.signal.bigwig --outFileName data/log2ratio_hcc_h3k4me3.bigwig
</pre>

Next, we create the gene lists files (e.g. files containing upregulated/downregulated gene symbols) in the data/ directory. We do so as follows:

<pre>
awk -F'\t' -v c="eumyc_myc" 'NR==1{for (i=1; i<=NF; i++) if ($i==c){p=i; break}; next} {print $p}' ../output/mouse_de/de_genes_down_symbols.txt > data/eumyc_down.txt
awk -F'\t' -v c="eumyc_myc" 'NR==1{for (i=1; i<=NF; i++) if ($i==c){p=i; break}; next} {print $p}' ../output/mouse_de/de_genes_up_symbols.txt > data/eumyc_up.txt
awk -F'\t' -v c="liver_myc" 'NR==1{for (i=1; i<=NF; i++) if ($i==c){p=i; break}; next} {print $p}' ../output/mouse_de/de_genes_down_symbols.txt > data/hcc_down.txt
awk -F'\t' -v c="liver_myc" 'NR==1{for (i=1; i<=NF; i++) if ($i==c){p=i; break}; next} {print $p}' ../output/mouse_de/de_genes_up_symbols.txt > data/hcc_up.txt
</pre>

Finally, run the analysis to make heatmaps and metagene plots:

<pre>./make_profile_plots.sh</pre>

# Super-Enhancer (SE) Analysis

## dbSUPER

The data was obtained from https://asntech.org/dbsuper/

* The zip file containing the BED files for the mm9 SEs can be found here: https://github.com/Yenaled/felsher/blob/master/chipseq/data/all_mm9_bed.zip -- it was originally obtained via the following link: https://asntech.org/dbsuper/data/bed/mm9/all_mm9_bed.zip
* The mapping between SE ID and gene symbol can be found here: https://github.com/Yenaled/felsher/blob/master/chipseq/data/superenhancer_gene_annotations.csv -- it was originally obtained via going to the <a href="https://asntech.org/dbsuper/adv_search.php">Detailed Search</a> feature of the dbSUPER website, selecting "Mouse (mm9)" as the genome, and clicking Search.

The following commands were used to convert further process the files above (convert the mm9 to mm10 coordinates via liftOver; combine gene symbol, tissue name, and SE coordinates into a single file; remove enhancers not mapped to specific genes):

<pre>Rscript process_bed_zip.r data/all_mm9_bed.zip data/dbSUPER_mm9.bed data/superenhancer_gene_annotations.csv</pre>

Next, feed the resulting data/dbSUPER_mm9.bed output file into UCSC's liftOver tool: https://genome.ucsc.edu/cgi-bin/hgLiftOver setting the options as follows. Original Genome: Mouse; Original Assembly: "July 2007 (NCBI37/mm9)", New Genome: Mouse, New Assembly: "Dec. 2011 (GRCm38/mm10)". Then download the resulting file (now in mm10 coordinates) as data/dbSUPER_mm10.bed

The final mm10 file of the SE records can be found here: https://github.com/Yenaled/felsher/blob/master/chipseq/data/dbSUPER_mm10.bed

## SEA

The data was obtained from http://sea.edbc.org/

* The mouse mm10 SEs BED file can be found (gzip-compressed) here: https://github.com/Yenaled/felsher/blob/master/chipseq/data/SEA00201.bed.gz -- it was originally obtained by navigating to "Download Data" via http://sea.edbc.org/ and then selecting SEA00201 under DataSourceID (Date "2018-06-12"; Species "Mouse"; Description "mouse super-enhancers by SEA").

Finally, execute the following command to finish processing the BED file (i.e. remove unnecessary columns):

<pre>zcat < data/SEA00201.bed.gz|awk -v OFS="\t" -F"\t" '{print $2,$3,$4,$16,$7}' > data/SEA00201_processed.bed</pre>

The final SE file can be found here: https://github.com/Yenaled/felsher/blob/master/chipseq/data/SEA00201_processed.bed

## Make files for graphing

<pre>Rscript merge_bed_with_genes.r data/dbSUPER_mm10.bed dbSUPER_hcc_up.bed data/hcc_up.txt
Rscript merge_bed_with_genes.r data/dbSUPER_mm10.bed dbSUPER_hcc_down.bed data/hcc_down.txt
Rscript merge_bed_with_genes.r data/dbSUPER_mm10.bed dbSUPER_eumyc_up.bed data/eumyc_up.txt
Rscript merge_bed_with_genes.r data/dbSUPER_mm10.bed dbSUPER_eumyc_down.bed data/eumyc_down.txt

Rscript merge_bed_with_genes.r data/SEA00201_processed.bed SEA_hcc_up.bed data/hcc_up.txt
Rscript merge_bed_with_genes.r data/SEA00201_processed.bed SEA_hcc_down.bed data/hcc_down.txt
Rscript merge_bed_with_genes.r data/SEA00201_processed.bed SEA_eumyc_up.bed data/eumyc_up.txt
Rscript merge_bed_with_genes.r data/SEA00201_processed.bed SEA_eumyc_down.bed data/eumyc_down.txt</pre>
