#!/bin/bash

###############################################################
### Description: Generates indices for alignment, searches 
###              for FastQ files and aligns them, and 
###              performs QC and quant.
### Written by Delaney Sullivan
###############################################################

num_threads=1
align_prog="STAR"
index=""
genome_file=""
secondpass=0
readlen=0
use_trimmed=0
kallisto_num_bootstraps=100
while getopts ":t:hi:p:g:l:sr" opt; do
	case $opt in
		h) 
		   echo "Usage: $0 [-h] [-g genome_file(s)] [-p align_prog] [-i index] [-l read_len] [-t num_threads] [-s] <annot_file> <input_path>"
		   echo "Options and arguments:"
		   echo "  -h                : Help"
                   echo "  -g                : The reference genome (or transcriptome) fasta file(s) to be used."
                   echo "                      Multiple files should be quoted and space-separated, e.g. \"A.fa B.fa\""
                   echo "  -p                : The alignment program (options: STAR or kallisto)"
                   echo "  -i                : The index to be used (a file for kallisto; a directory for STAR)"
                   echo "                    if not supplied, will be generated."
                   echo "  -l                : Sequencing read length (needed for generating index with STAR)"
                   echo "  -t                : Number of threads (Default: 1)"
                   echo "  -s                : Use this option if you want to use two-pass "
                   echo "                    alignment for STAR."
                   echo "  -r                : Use this option if you did adapter trimming and want to do "
                   echo "                    mapping using the trimmed fastq files."
                   echo "  annot_file        : The GTF annotation file to be used"
                   echo "  input_path        : Directory path containing the FastQ files"
		   exit 1
		   ;;
		g) genome_file="$OPTARG"
		  ;;
		t) num_threads="$OPTARG"
		  ;;
		p) align_prog="$OPTARG"
		  ;;
		i) index="$OPTARG"
		  ;;
		l) readlen="$OPTARG"
		  ;;
                s) secondpass=1
                  ;;
		r) use_trimmed=1
		  ;;
		\?)
		   echo "Invalid option: -$OPTARG" >&2 # Fold output from stdout to stderr (file_descriptor=2)
		   exit 1
		   ;;
	esac
done

shift $((OPTIND-1))

# Some preliminary checks

if [ ! "$#" -eq 2 ]; then
	echo "Usage: $0 [-h] [-g genome_file] [-p align_prog] [-i index] [-t num_threads] [-s] <annot_file> <input_path>" >&2
	exit 1
fi

if [ ! "$align_prog" == "STAR" -a ! "$align_prog" == "kallisto" ]; then
	echo "Must choose either STAR or kallisto" >&2
	exit 1
fi

annot_file="$1"
input_path="$2"

if [ ! -f "$annot_file" ]; then
        echo "Annotation file does not exist: $annot_file" >&2
        exit 1
fi

if [ ! -d "$input_path" ]; then
	echo "The following directory, specified for <input_path>, does not exist: $input_path" >&2
        exit 1
fi

# Function to output and execute commands

cmdexec() {
	cmd="$1"
        BEGINHIGHLIGHT='\033[0;35m'
	ENDHIGHLIGHT='\033[0m'
	printf "${BEGINHIGHLIGHT}$cmd${ENDHIGHLIGHT}\n"
	eval "$cmd"
}


# Now, build the index

## For STAR:
if [ "$align_prog" == "STAR" -a "$index" == "" ]; then
	declare -i overhang
	overhang="$readlen"
	if [ $overhang -eq 0 ]; then
		echo "Please specify a valid sequencing read length (read_len) for STAR index generation" >&2
		exit 1
	fi
	if [ "$genome_file" == "" ]; then
                echo "Please specify a valid reference genome file that exists" >&2
                exit 1
	fi
	overhang=$overhang-1
	cmdexec "mkdir -p '${input_path}/genomedir'"
	cmdexec "STAR --runThreadN '$num_threads' --runMode genomeGenerate --sjdbOverhang '${overhang}' --genomeDir '${input_path}/genomedir' --genomeFastaFiles ${genome_file} --sjdbGTFfile '${annot_file}'"
	index="${input_path}/genomedir"
fi

## For kallisto:
if [ "$align_prog" == "kallisto" -a "$index" == "" ]; then
	transcriptome_files=($genome_file)
	if [ ${#transcriptome_files[@]} -gt 1 ]; then
		cmdexec "tmpfile=$(mktemp)"
		for i in "${transcriptome_files[@]}"
		do
			cmdexec "cat '$i' >> '$tmpfile'"
		done
		cmdexec "kallisto index -i '${input_path}/kallisto_transcripts.idx' '${tmpfile}' &> '${input_path}/index_kallisto.log'"
		cmdexec "rm '${tmpfile}'"
	else
                cmdexec "kallisto index -i '${input_path}/kallisto_transcripts.idx' '${genome_file}' &> '${input_path}/index_kallisto.log'"
	fi
	index="${input_path}/kallisto_transcripts.idx"
fi

# Do the alignment

do_alignment() {
	align_prog="$1"
	input_path="$2"
	extension="$3"
	use_trimmed="$4"
	output_folder_name="$5"
	options="$6"
	annot_file="$7"
	num_threads="$8"

	files=""
	if [ $use_trimmed -eq 1 ]; then
		extension="_trimmed$extension"
		files=$(find "$input_path" -name "*$extension" -exec echo "{}" \;|sort)
	else
		files=$(find "$input_path" -name "*$extension" -not -name "*_trimmed$extension" -exec echo "{}" \;|sort)
	fi
	files=($files)
	numfiles=${#files[@]}
	suffix_trim=$(echo "$extension"|wc -m)
        output_folder="$input_path/$output_folder_name/"
	bam_suffix="Aligned.sortedByCoord.out.bam"

	files_processed=() # Keep track of which files we have already processed
	for (( i=0; i<${numfiles}; i++ )); do
		file=${files[$i]}
		if [[ $file == *_1"$extension" ]]; then
			suffix=_1"$extension"
			suffix_length=${#suffix}
			prefix="${file::-$suffix_length}"
			file_2="${prefix}"_2"$extension"
			prefix=$(basename "$prefix")
			if [ -f "$file_2" ]; then # Paired end
				files_processed+=($file)
				files_processed+=($file_2)
				file_1=$(echo "$file"|rev|cut -c ${suffix_trim}-|rev)
				file_2=$(echo "$file_2"|rev|cut -c ${suffix_trim}-|rev)
        	                if [ "$align_prog" == "STAR" ]; then
                	                cmdexec "STAR ${options} --readFilesIn '${file_1}$extension' '${file_2}$extension' --outFileNamePrefix '${output_folder}${prefix}'"
					if [ -f "${output_folder}${prefix}${bam_suffix}" ]; then
						cmdexec "featureCounts -p -T $num_threads -a '${annot_file}' -o '${output_folder}${prefix}FeatureCounts.txt' '${output_folder}${prefix}${bam_suffix}'"
					fi
	                        fi
				if [ "$align_prog" == "kallisto" ]; then
                                	cmdexec "mkdir -p '${output_folder}/${prefix}'"
					cmdexec "kallisto quant ${options} -o '${output_folder}/${prefix}' '${file_1}$extension' '${file_2}$extension' &> '${output_folder}/${prefix}/${prefix}_kallisto.log'"
				fi
			fi
		fi
		if [[ ! " ${files_processed[@]} " =~ " ${file} " ]]; then
			files_processed+=($file)
			suffix="$extension"
                        suffix_length=${#suffix}
                        prefix="${file::-$suffix_length}"
                        prefix=$(basename "$prefix")
                        file=$(echo "$file"|rev|cut -c ${suffix_trim}-|rev)
                        if [ "$align_prog" == "STAR" ]; then
                	        cmdexec "STAR ${options} --readFilesIn '${file}$extension' --outFileNamePrefix '${output_folder}${prefix}'"
				if [ -f "${output_folder}${prefix}${bam_suffix}" ]; then
					cmdexec "featureCounts -T $num_threads -a '${annot_file}' -o '${output_folder}${prefix}FeatureCounts.txt' '${output_folder}${prefix}${bam_suffix}'"
				fi
			fi
			if [ "$align_prog" == "kallisto" ]; then
				### NOTE: Change -l and -s estimates for single end reads
				cmdexec "mkdir -p '${output_folder}/${prefix}'"
                                cmdexec "kallisto quant ${options} --single -l 200 -s 20 -o '${output_folder}/${prefix}' '${file}$extension' &> '${output_folder}/${prefix}/${prefix}_kallisto.log'"
			fi
		fi
	done
}

## STAR
if [ "$align_prog" == "STAR" ]; then
        output_folder_name="STAR_aligned"
	if [ $secondpass -eq 0 ]; then
		cmdexec "mkdir -p '$input_path/$output_folder_name/'"
		options_to_include="--quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --runThreadN $num_threads --genomeDir '$index'"
		do_alignment "$align_prog" "$input_path" ".fq" $use_trimmed "$output_folder_name" "$options_to_include" "$annot_file" "$num_threads"
		do_alignment "$align_prog" "$input_path" ".fq.gz" $use_trimmed "$output_folder_name" "$options_to_include --readFilesCommand zcat" "$annot_file" "$num_threads"
		do_alignment "$align_prog" "$input_path" ".fastq" $use_trimmed "$output_folder_name" "$options_to_include" "$annot_file" "$num_threads"
		do_alignment "$align_prog" "$input_path" ".fastq.gz" $use_trimmed "$output_folder_name" "$options_to_include --readFilesCommand zcat" "$annot_file" "$num_threads"
	fi
	if [ $secondpass -eq 1 ]; then
		output_folder_name="STAR_firstpass"
                cmdexec "mkdir -p '$input_path/$output_folder_name/'"
		options_to_include="--runThreadN $num_threads --genomeDir '$index'"
                do_alignment "$align_prog" "$input_path" ".fq" $use_trimmed "$output_folder_name" "$options_to_include" "$annot_file" "$num_threads"
                do_alignment "$align_prog" "$input_path" ".fq.gz" $use_trimmed "$output_folder_name" "$options_to_include --readFilesCommand zcat" "$annot_file" "$num_threads"
                do_alignment "$align_prog" "$input_path" ".fastq" $use_trimmed "$output_folder_name" "$options_to_include" "$annot_file" "$num_threads"
                do_alignment "$align_prog" "$input_path" ".fastq.gz" $use_trimmed "$output_folder_name" "$options_to_include --readFilesCommand zcat" "$annot_file" "$num_threads"

		# For two pass, get the tab files
		tab_files=$(find "$input_path/STAR_firstpass" -name "*.out.tab" -exec printf "'{}' " \;)

		output_folder_name="STAR_aligned"
                options_to_include="--quantMode GeneCounts --limitSjdbInsertNsj 3000000 --outSAMtype BAM SortedByCoordinate --runThreadN $num_threads --genomeDir '$index' --sjdbFileChrStartEnd ${tab_files}"
                output_folder_name="STAR_aligned"
                cmdexec "mkdir -p '$input_path/$output_folder_name/'"
                do_alignment "$align_prog" "$input_path" ".fq" $use_trimmed "$output_folder_name" "$options_to_include" "$annot_file" "$num_threads"
                do_alignment "$align_prog" "$input_path" ".fq.gz" $use_trimmed "$output_folder_name" "$options_to_include --readFilesCommand zcat" "$annot_file" "$num_threads"
                do_alignment "$align_prog" "$input_path" ".fastq" $use_trimmed "$output_folder_name" "$options_to_include" "$annot_file" "$num_threads"
                do_alignment "$align_prog" "$input_path" ".fastq.gz" $use_trimmed "$output_folder_name" "$options_to_include --readFilesCommand zcat" "$annot_file" "$num_threads"
	fi

	# QC:
	cmdexec "multiqc -f -o '$input_path/multiqc_alignment' -c ./rnaseq_config.yaml '$input_path'"
fi

## Kallisto
if [ "$align_prog" == "kallisto" ]; then
	output_folder_name="kallisto_aligned"
	cmdexec "mkdir -p '$input_path/$output_folder_name/'"
	options_to_include="-i $index -b $kallisto_num_bootstraps -t $num_threads" # --genomebam -g '${annot_file}'
	do_alignment "$align_prog" "$input_path" ".fq" $use_trimmed "$output_folder_name" "$options_to_include" "$annot_file" "$num_threads"
        do_alignment "$align_prog" "$input_path" ".fq.gz" $use_trimmed "$output_folder_name" "$options_to_include" "$annot_file" "$num_threads"
        do_alignment "$align_prog" "$input_path" ".fastq" $use_trimmed "$output_folder_name" "$options_to_include" "$annot_file" "$num_threads"
        do_alignment "$align_prog" "$input_path" ".fastq.gz" $use_trimmed "$output_folder_name" "$options_to_include" "$annot_file" "$num_threads"

	# QC:
	cmdexec "multiqc -f -o '$input_path/multiqc_alignment' -c ./rnaseq_config.yaml '$input_path'"
fi
