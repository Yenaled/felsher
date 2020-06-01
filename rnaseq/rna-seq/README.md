# rna-seq

Pipeline for RNA-Seq analysis. Supports both STAR alignment and kallisto pseudoalignment.

## Versions
<ul>
  <li>STAR 2.5.4b</li>
  <li>kallisto 0.44.0</li>
  <li>featureCounts 1.6.3</li>
  <li>MultiQC 1.7</li>
  <li>DESeq2 1.22.2</li>
  <li>sleuth 0.30.0</li>
  <li>fastq-dump 2.8.2</li>
</ul>

# Installation

First, clone this repository by opening a terminal, logging into your server if you're doing processing remotely (e.g. if you're using Sherlock, log into Sherlock), and then enter the following command:
<pre>git clone https://github.com/Yenaled/felsher/rnaseq/rna-seq</pre>

# Genome files and annotations

GENCODE Mouse (mm10): (release 15) - https://www.gencodegenes.org/mouse/
<ul>
  <li>Transcriptome: gencode.vM15.pc_lncRNA_combined_transcripts.fa - Merged the "Protein-coding transcript sequences" and the "Long non-coding RNA transcript sequences"</li>
  <li>Genome: GRCm38.primary_assembly.genome.fa - Genome sequence, primary assembly (GRCm38)</li>
  <li>Annotations: gencode.vM15.primary_assembly.annotation.gtf - Comprehensive gene annotation	(Primary assembly)</li>
  <li>Annotations (w/ human MYC): gencode.vM15.primary_assembly.annotation_edited.gtf - The previous file merged with hmyc_final.gtf</li>
  <li>Human MYC annotation: hmyc_final.gtf - The final annotation for human MYC (chromosome:GRCh38:8:127734834:127742034:1) with all coordinates shifted left by 127735433 bases and chromosome names changed to hMYC (instead of chr8)</li>
  <li>Human MYC transcript FASTA: hmyc_transcript.fa</li>
  <li>Human MYC genome FASTA: hmyc.fa</li>
</ul>

# Downloading RNA-seq data

There are many ways to download RNA-seq data. Below is a way to do so (via fastq-dump) [note: this is neither the most efficient nor the best way to obtain sequencing files]. We'll download 6 sequencing files from GSE106078 (first three are normal; next three are tumor) into the folder: ./tall_rnaseq:

<pre>fastq-dump --skip-technical --readids --read-filter pass --split-3 --clip --outdir ./tall_rnaseq SRR6206939
fastq-dump --skip-technical --readids --read-filter pass --split-3 --clip --outdir ./tall_rnaseq SRR6206940
fastq-dump --skip-technical --readids --read-filter pass --split-3 --clip --outdir ./tall_rnaseq SRR6206941
fastq-dump --skip-technical --readids --read-filter pass --split-3 --clip --outdir ./tall_rnaseq SRR6206942
fastq-dump --skip-technical --readids --read-filter pass --split-3 --clip --outdir ./tall_rnaseq SRR6206943
fastq-dump --skip-technical --readids --read-filter pass --split-3 --clip --outdir ./tall_rnaseq SRR6206944</pre>

At this stage, it is very important to check that the files were downloaded successfully. If there's an error, you'd see a message such as "err: incomplete while reading file within network system module - failed SRR6206939". If successful, you'd see a message such as "Read 27684182 spots for SRR6206939", "Written 27684182 spots for SRR6206939". Make sure all sequence files were obtained successfully before proceeding.

The resulting downloaded files should be named something like SRR6206939_pass.fastq (if it were paired-end reads, it would be named SRR6206939_pass_1.fastq and SRR6206939_pass_2.fastq).

# Setup your raw data

Follow the instructions here for preliminary quality control and processing of your raw fastQ sequencing read files: https://github.com/Yenaled/felsher/blob/master/rnaseq/prep_fastq/README.md

# Usage

It's best to run this program as a batch job on a server and an sample template for configuring a batch job submission for the Sherlock computing cluster (at Stanford) is provided in rnaseq.batch. Edit that file as you see fit. The bottom of that file displays an example of the main commands that are going to be run for the alignment (rnaseq.sh) and the differential gene expression & normalization (de.r and de2.r):

kallisto:
<pre>./rnaseq.sh -t 10 -p kallisto -g "gencode.vM15.pc_lncRNA_combined_transcripts.fa hmyc_transcript.fa" -s gencode.vM15.primary_assembly.annotation_edited.gtf ./tall_rnaseq</pre>
<pre>Rscript ./de2.r 10 ./tall_rnaseq/kallisto_aligned/ tallmycon "SRR6206939_pass,SRR6206940_pass,SRR6206941_pass,SRR6206942_pass,SRR6206943_pass,SRR6206944_pass" "c,c,c,t,t,t" gencode.vM15.primary_assembly.annotation_edited.gtf</pre>
STAR:
<pre>./rnaseq.sh -t 10 -p STAR -g "GRCm38.primary_assembly.genome.fa hmyc.fa" -l 50 -s gencode.vM15.primary_assembly.annotation_edited.gtf ./tall_rnaseq</pre>
<pre>Rscript ./de.r 10 ./tall_rnaseq/STAR_aligned/ FeatureCounts.txt 0.05 0 tallmycon "SRR6206939_pass,SRR6206940_pass,SRR6206941_pass,SRR6206942_pass,SRR6206943_pass,SRR6206944_pass" treatment:c,c,c,t,t,t "~treatment" gencode.vM15.primary_assembly.annotation_edited.gtf NULL</pre>

You can view the help menu (for the alignment step) by running <b>./rnaseq.sh -h</b>

Further details on how to figure the commands are as follows:

<pre>rnaseq.sh -t &lt;number of processes&gt; -p &lt;kallisto or STAR&gt; -g &lt;genome/transcriptome FASTA files&gt; -l &lt;read length (if using STAR)&gt; -s &lt;GTF annotation file&gt; &lt;output folder&gt;</pre>
<pre>Rscript de2.r &lt;number of processes&gt; &lt;kallisto alignment directory&gt; &lt;output files prefix&gt; &lt;list of alignment file prefixes&gt; &lt;list of control (c) and treatment (t) matching the order of the supplied alignment files&gt; &lt;GTF annotation file&gt;</pre>
<pre>Rscript de.r &lt;number of processes&gt; &lt;STAR alignment directory&gt; &lt;suffix of gene quantification file which will be either FeatureCounts.txt or ReadsPerGene.out.tab&gt; &lt;FDR&gt; &lt;log2FC threshold&gt; &lt;output files prefix&gt; &lt;explanatory variables (see example)&gt; &lt;design formula&gt; &lt;GTF annotation file&gt; &lt;Number specifying which variable in the design to compare (or NULL for default)&gt;</pre>

## Job submission

Finally, when you're done editing rnaseq.batch, you can submit the job by typing in the command as follows:
<pre>sbatch rnaseq.batch</pre>

## Output files

All new files/folders will be created in the directory you specified to be your input directory. Within that directory, the directory <b>multiqc_alignment</b> contains the QC for the alignments while <b>STAR_aligned</b> and <b>kallisto_aligned</b> will contain the actual alignments, quantifications, and DESeq2/sleuth analyses. 

## Seeing which commands were run

The output of the program will have all commands that were executed highlighted in purple. On Sherlock, the output will generally be stored in a file, e.g. slurm-36996399.out (where 36996399 is the assigned job ID). If you have such a file, to view a list of all commands that were executed, run the following (replacing JobID with the numerical assigned job ID):
<pre>grep $'\033\[0;35m' slurm-JobID.out</pre>
