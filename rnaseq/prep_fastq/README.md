# prep_fastq
Here, we will generate quality control reports of FASTQ files as well as do additional processing of FASTQ files (e.g. adapter trimming, reads filtering, etc.). We will be leveraging fastQC, multiQC, and cutadapt for these tasks so make sure they are installed.
<details><summary>For Sherlock users: Click here to check installation and for installation instructions</summary>
Cutadapt and multiQC should be in $PI_HOME/python, therefore if you run the following:
<pre>
module load python/3.6.1
$PI_HOME/python/bin/cutadapt --help
$PI_HOME/python/bin/multiqc --help
</pre>
You should be able to view the cutadapt help page and the multiQC help page.
However, if they are not installed, install them by running the following commands:
<pre>
mkdir $PI_HOME/python
module load python/3.6.1
export PYTHONUSERBASE=$PI_HOME/python
pip3 install --user cutadapt
pip3 install --user multiqc
</pre>
</details>

## Versions
<ul>
  <li>Python 3.6.1</li>
  <li>FastQC 0.11.8</li>
  <li>TrimGalore 0.5.0</li>
  <li>Cutadapt 1.18</li>
  <li>MultiQC 1.7</li>
</ul>

# Installation

First, clone this repository by opening a terminal, logging into your server if you're doing processing remotely (e.g. if you're using Sherlock, log into Sherlock), and then enter the following command:
<pre>git clone https://github.com/Yenaled/felsher/rnaseq/prep_fastq</pre>

# Usage

It's best to run this program as a batch job on a server and an sample template for configuring a batch job submission for Sherlock is provided in <b>prep_fastq.batch</b>. Edit that file as you see fit. The bottom of that file displays an example of the main command that's going to be run.
<pre>./processfq.sh -t 4 -g /path/to/sequencing/reads --stringency 3</pre>
You can view the help menu by running ./processfq.sh -h.

## trim_galore

The above command searches for fastQ files recursively in /path/to/sequencing/reads and generates QC reports using 4 threads (i.e. multithreading) as well as trims adapter sequences using <b>trim_galore</b> with all its default options (except with stringency set to 3). More details about the usage of trim_galore can be found here: https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md
<p>You can view trim_galore's help menu by entering: trim_galore --help. (Also, trim_galore should already be installed as a module on Sherlock -- but if not, visit here: https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/).</p>

## cutadapt

<p>If you want to use <b>cutadapt</b> for adapter trimming, you can remove the <b>-g</b> option and enter cutadapt options like the following:</p>
<pre>./processfq.sh -t 4 /path/to/sequencing/reads --cores=4 -g GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -e 0.05 -O 8 -q 10</pre>
Essentially, the above command searches for fastQ files recursively in /path/to/sequencing/reads and generates QC reports using 4 threads (i.e. multithreading) as well as trims (using 4 cores) the 5' adapter sequence, GATCGGAAGAGCACACGTCTGAACTCCAGTCAC, with an error tolerance of 5% and a minimum overlap of 8; additionally, the bases on the 3' end will be trimmed with a read quality threshold of 10. You can view all the possible options for adapter trimming in the cutadapt program's official documentation: https://cutadapt.readthedocs.io/en/stable/guide.html.

## QC analysis only

<p>If you only want to do QC analysis (and not adapter/quality trimming via cutadapt or trim_galore), you can simply omit the cutadapt options and omit the -g option as well, e.g.:</p>
<pre>./processfq.sh -t 4 /path/to/sequencing/reads</pre>

## Job submission

Finally, when you're done editing prep_fastq.batch, you can submit the job by typing in the command as follows:
<pre>sbatch prep_fastq.batch</pre>
The integrated final QC report will be located in a folder called <b>multiqc</b> in the path where you specified your fastQ files are contained. You can view the report in that folder via the file: <b>multiqc_report.html</b>. All the newly trimmed sequencing files will have the same filename as the original sequencing file but will be suffixed with _trimmed (e.g. file_trimmed.fq.gz). Individual HTML and Zip files for each individual sample's QC report are also provided.

## Note for paired-end reads

Paired-end reads (two fastQ files) should be placed in the same directory and named with the suffix <b>_1.fq</b> and <b>_2.fq</b> (.fq can be .fastq, .fq.gz, or fastq.gz but the _1 and _2 labels are absolutely required to denote that the reads are paired end). 

## Seeing which commands were run

The output of the program will have all commands that were executed highlighted in purple. On Sherlock, the output will generally be stored in a file, e.g. slurm-36996399.out (where 36996399 is the assigned job ID). If you have such a file, to view a list of all commands that were executed, run the following (replacing JobID with the numerical assigned job ID):
<pre>grep $'\033\[0;35m' slurm-JobID.out</pre>
