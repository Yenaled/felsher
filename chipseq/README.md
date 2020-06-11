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
bigwigCompare --bigwig1 data/eumyc_h3k27ac_C.fc.signal.bigwig --bigwig2 data/eumyc_h3k27ac_T.fc.signal.bigwig --outFileName data/log2ratio_eumyc_h3k27ac.bigwig
bigwigCompare --bigwig1 data/eumyc_h3k4me3_C.fc.signal.bigwig --bigwig2 data/eumyc_h3k4me3_T.fc.signal.bigwig --outFileName data/log2ratio_eumyc_h3k4me3.bigwig
bigwigCompare --bigwig1 data/hcc_h3k27ac_C.fc.signal.bigwig --bigwig2 data/hcc_h3k27ac_T.fc.signal.bigwig --outFileName data/log2ratio_hcc_h3k27ac.bigwig
bigwigCompare --bigwig1 data/hcc_h3k4me3_C.fc.signal.bigwig --bigwig2 data/hcc_h3k4me3_T.fc.signal.bigwig --outFileName data/log2ratio_hcc_h3k4me3.bigwig
</pre>

Next, we create the gene lists files (e.g. files containing upregulated/downregulated genes) in the data/ directory. We do so as follows:

<pre>
awk -F'\t' -v c="eumyc_myc" 'NR==1{for (i=1; i<=NF; i++) if ($i==c){p=i; break}; next} {print $p}' ../output/mouse_de/de_genes_down_symbols.txt > data/eumyc_down.txt
</pre>

Put all the genes you're interested in studying in a file (genes.txt) with <i>n</i> lines where <i>n</i> is your number of genes. The gene symbols in that file should match those in refGene.bed. Then run the following to filter the reference genome .bed file so that it only contains your genes of interest:
<pre>
fname="genes.txt"
output="filtered_refGene.bed"

while read -r line
do
        echo "$line""$(printf '\t')"
done < "$fname" > "/tmp/genes.txt"
cat "/tmp/genes.txt"|sed '/^[[:space:]]*$/d'|tr -d '\r' > "/tmp/genes_pattern.txt"
grep -f /tmp/genes_pattern.txt refGene.bed|awk '!seen[$0]++' > "$output"
</pre>

<p>Inspect the newly created file (in the example above, filtered_refGene.bed) -- it should have only your genes of interest there. Note that some genes may appear multiple times due to having multiple transcription start sites (TSS's)</p>

<p>Find the .fc.signal.bigwig fold change signal track file; the path to it might look something like call-macs2_pooled/execution.</p>

First, we want to deal with genes that have multiple TSS's. In the folder with the bigwig file, open up a python session and insert the following (obviously modifying the file names and whatnot):
<pre>
input_bigwig_file = "name_of_file.fc.signal.bigwig"
input_bed_file =  "filtered_refGene.bed"
output_file = "filtered_nodups_refGene.bed"
padding=3000

#####
import pyBigWig
import csv
bw = pyBigWig.open(input_bigwig_file)
output = {}
for line in open(input_bed_file):
        cols = line.strip().split()
        try:
                if (cols[5] == "-"):
                        tss_start = max(int(cols[2]) - padding, 0)
                        tss_end = int(cols[2]) + padding
                else:
                        tss_start = max(int(cols[1]) - padding, 0)
                        tss_end = int(cols[1]) + padding
                avg = bw.stats(cols[0], tss_start, tss_end)[0]
		if avg == None:
			avg = 0
        except:
                print("Skipping the following:")
                print(cols)
                continue
        if cols[3] in output:
                current_avg = int(output[cols[3]][3])
                if avg > current_avg:
                        output[cols[3]] = [cols[0], cols[1], cols[2], avg, cols[4], cols[5]]
        else:
                output[cols[3]] = [cols[0], cols[1], cols[2], avg, cols[4], cols[5]]

bw.close()
with open(output_file, 'w') as outfile:
	csv_writer = csv.writer(outfile, delimiter='\t')
	for k,v in output.items():
		o = csv_writer.writerow([v[0], v[1], v[2], k, v[4], v[5]])
</pre>

<p>This will create a new .bed file (in this case, filtered_nodups_refGene.bed) that does not have duplicate gene names because if a gene appears multiple times due to multiple TSS's, only the one with the highest average enrichment for the signal in its TSS (+/- padding) will be used. Use that new .bed file (instead of the old filtered_refGene.bed) in the following steps.</p>
In that folder with your bigwig signal track file, run computeMatrix to generate a matrix file. Example (where we create a matrix based on +/- 2000 bp of the TSS):</p>
<pre>computeMatrix reference-point --referencePoint TSS -S SRR3020034_pass_1_trimmed.merged.nodup_x_ctl_for_rep1.fc.signal.bigwig -R filtered_nodups_refGene.bed --beforeRegionStartLength 2000 --afterRegionStartLength 2000 --skipZeros -o matrix.mat.gz</pre>
Then plot graphs like:
<pre>plotHeatmap -m matrix.mat.gz -out heatmap.pdf --heatmapWidth 15</pre>

# Step 4: Annotate Peak Calls

The peak calls might be stored in a file ending in Peak.gz (e.g. for replicated transcription factor ChIP, you might see files like optimal_peak.narrowPeak.gz in a folder like call-reproducibility_idr/execution/, which contains peaks identified from IDR analysis). You can view more info about the various peak files here: https://www.encodeproject.org/data-standards/chip-seq/. To analyze the peak files and to assign peaks to genes, use the annotatePeaks.pl utility of HOMER (http://homer.ucsd.edu/homer/). Be sure to unzip the .gz peaks file prior to using HOMER.
