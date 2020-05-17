#!/bin/sh

###############################################################
### File: normalize.sh
### Description: Takes GDC/TCGA data, that have been processed
###              and formatted by processgdc.r, and obtains
###              normalized RNA-seq count values via
###              the script: normalize.r.
### Written by Delaney Sullivan
###############################################################

SCRIPT_PATH=$(dirname $(realpath -s "$0"))
norm_program="$SCRIPT_PATH/normalize.r"

while getopts ":h" opt; do
case $opt in
h)
echo "Usage: $0 [-h] <project>"
echo "Options and arguments:"
echo "  -h              : Help"
echo "  project         : The directory within the processed folder,"
echo "                    containing the processed GDC data"
exit 1
;;
\?)
echo "Invalid option: -$OPTARG" >&2 # Fold output from stdout to stderr (file_descriptor=2)
exit 1
;;
esac
done

shift $((OPTIND-1))

if [ "$#" -ne 1 ]; then
echo "Usage: $0 [-h] <project>"
exit 1
fi

Rscript $norm_program "$1" $SCRIPT_PATH
