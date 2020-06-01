#!/bin/bash

###############################################################
### Description: Searches for FastQ files and processes them
###              by generating QC metrics, trimming adapters,
###              and doing read quality filtering.
### Written by Delaney Sullivan
###############################################################

# Process command line arguments
num_threads=1
trim_prog="cutadapt"
while getopts ":t:hg" opt; do
	case $opt in
		h) 
		   echo "Usage: $0 [-h] [-g] [-t num_threads] <input_path> [cutadapt_options]"
		   echo "Options and arguments:"
		   echo "  -h                : Help"
                   echo "  -t                : Number of threads (number of files that "
                   echo "                    will be process simultaneously) for FastQC."
                   echo "                    (Default: 1)."
                   echo "  -g                : Use this option if you want to use trim_galore "
                   echo "                    instead of cutadapt."
		   echo "  cutadapt_options  : Options that will be passed to cutadapt."
                   echo "                    For more info, see: cutadapt --help."
                   echo "                    Do not use the option -o as we will be "
                   echo "                    generating output file names automatically."
		   echo "  input_path        : Directory path containing the FastQ files"
		   exit 1
		   ;;
		t) num_threads="$OPTARG"
		  ;;
		g) trim_prog="trim_galore"
		  ;;
		\?)
		   echo "Invalid option: -$OPTARG" >&2 # Fold output from stdout to stderr (file_descriptor=2)
		   exit 1
		   ;;
	esac
done

shift $((OPTIND-1))

if [ "$#" -eq 0 ]; then
	echo "Usage: $0 [-h] [-g] [-t num_threads] <input_path> [cutadapt_options]" >&2
	exit 1
fi

# Make sure the options -o and -p are not included for cutadapt
cutadapt_options="${@:2}"
for item in "$cutadapt_options"; do
	if [[ "-o" == "$item" ]]; then
		   echo "Invalid option for cutadapt: -o" >&2
                   exit 1
	fi
        if [[ "-p" == "$item" ]]; then
                   echo "Invalid option for cutadapt: -p" >&2
                   exit 1
        fi

done

# Make sure input_path exists
input_path="$1"
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

# Trimming function

trimfiles() {
	trim_prog="$1"
	input_path="$2"
	extension="$3"
	suffix_trim=$(echo "$3"|wc -m)
	cutadapt_options="$4"
	files=$(find "$input_path" -name "*$extension" -not -name "*_trimmed$extension" -exec echo "{}" \;|sort)
	files=($files)
	numfiles=${#files[@]}

	files_processed=() # Keep track of which files we have already processed
	for (( i=0; i<${numfiles}; i++ )); do
		file=${files[$i]}
		folder=$(dirname "$file")
		if [[ $file == *_1"$extension" ]]; then
			suffix=_1"$extension"
			suffix_length=${#suffix}
			file_2="${file::-$suffix_length}"_2"$extension"
			if [ -f "$file_2" ]; then # Paired end
				files_processed+=($file)
				files_processed+=($file_2)
				file_1=$(echo "$file"|rev|cut -c ${suffix_trim}-|rev)
				file_2=$(echo "$file_2"|rev|cut -c ${suffix_trim}-|rev)
        	                if [ "$trim_prog" == "cutadapt" ]; then
                	                cmdexec "cutadapt $cutadapt_options -o '${file_1}_trimmed$extension' -p '${file_2}_trimmed$extension' '$file_1$extension' '$file_2$extension'"
	                        fi
				if [ "$trim_prog" == "trim_galore" ]; then
                                	cmdexec "trim_galore --paired $cutadapt_options --output_dir '$folder' '$file_1$extension' '$file_2$extension'"
                                        cmdexec "mv '${file_1}_val_1$extension' '${file_1}_trimmed$extension'"
                                        cmdexec "mv '${file_2}_val_2$extension' '${file_2}_trimmed$extension'"
				fi 
			fi
		fi
		if [[ ! " ${files_processed[@]} " =~ " ${file} " ]]; then
			files_processed+=($file)
			file=$(echo "$file"|rev|cut -c ${suffix_trim}-|rev)
                        if [ "$trim_prog" == "cutadapt" ]; then
                	        cmdexec "$cutadapt $cutadapt_options -o '${file}_trimmed$extension' '$file$extension'"
			fi
                        if [ "$trim_prog" == "trim_galore" ]; then
				cmdexec "trim_galore $cutadapt_options --output_dir '$folder' '$file$extension'"
			fi
		fi
	done
}

# Do cutadapt (if user entered options for cutadapt)

if [ ! -z "$cutadapt_options" ] || [ "$trim_prog" == "trim_galore" ]; then
	trimfiles "$trim_prog" "$input_path" ".fastq" "$cutadapt_options"
	trimfiles "$trim_prog" "$input_path" ".fq" "$cutadapt_options"
	trimfiles "$trim_prog" "$input_path" ".fastq.gz" "$cutadapt_options"
	trimfiles "$trim_prog" "$input_path" ".fq.gz" "$cutadapt_options"
fi

# Generate FastQC reports

cmdexec "find '$input_path' -name '*.fq' | xargs --no-run-if-empty --verbose fastqc -t $num_threads"
cmdexec "find '$input_path' -name '*.fastq' | xargs --no-run-if-empty --verbose fastqc -t $num_threads"
cmdexec "find '$input_path' -name '*.fq.gz' | xargs --no-run-if-empty --verbose fastqc -t $num_threads"
cmdexec "find '$input_path' -name '*.fastq.gz' | xargs --no-run-if-empty --verbose fastqc -t $num_threads"

# Generate an integrated multiQC report

cmdexec "mkdir -p '$input_path/multiqc'"
cmdexec "multiqc -f -o '$input_path/multiqc' -c ./prep_fastq_config.yaml '$input_path'"
