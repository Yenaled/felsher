#!/bin/sh

###############################################################
### File: process_gdc.sh
### Description: Processes GDC gene data that was already
###              organized by organize_gdcdata.sh by calling
###              the R script: processgdc.r to get gene
###              expression values, among other types of
###              data.
### Written by Delaney Sullivan
###############################################################

SCRIPT_PATH=$(dirname $(realpath -s "$0"))
processgdc_program="$SCRIPT_PATH/processgdc.r"
control="Solid Tissue Normal"
experimental="Primary Tumor"

while getopts ":h" opt; do
	case $opt in
		h) 
		   echo "Usage: $0 [-h] <project> [control] [experimental]"
		   echo "Options and arguments:"
		   echo "  -h              : Help"
		   echo "  project         : The directory within the organized folder,"
                   echo "                    containing the organized GDC data"
		   echo "  control         : Folder name for the control samples"
                   echo "                    (default: 'Solid Tissue Normal')"
		   echo "  experimental    : Folder name for the experimental samples"
                   echo "                    (default: 'Primary Tumor')"
		   exit 1
		   ;;
		\?)
		   echo "Invalid option: -$OPTARG" >&2 # Fold output from stdout to stderr (file_descriptor=2)
		   exit 1
		   ;;
	esac
done

shift $((OPTIND-1))

if [ "$#" -eq 0 -o "$#" -ge 4 ]; then
	echo "Usage: $0 [-h] <project> [control] [experimental]" >&2
	exit 1
fi

processgdc_project="$1"

if [ "$#" -ge 2 ]; then
	control="$2"
fi

if [ "$#" -eq 3 ]; then
	experimental="$3"
fi

output_dir="$SCRIPT_PATH/processed/""$processgdc_project"

mkdir -p "$output_dir"
Rscript $processgdc_program $processgdc_project $SCRIPT_PATH "$control" "$experimental" $output_dir
