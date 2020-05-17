#!/bin/sh

###############################################################
### File: organize_gdcdata.sh
### Description: Organizes GDC data that was already downloaded
###              by putting all samples from a single
###              a patient into one directory.
### Written by Delaney Sullivan
###############################################################

while getopts ":h" opt; do
	case $opt in
		h) 
		   echo "Usage: $0 [-h] <data_directory>"
		   echo "Options and arguments:"
		   echo "  -h              : Help"
		   echo "  data_directory  : Directory (within tcga/data/) containing the downloaded GDC data"
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
	echo "Usage: $0 [-h] <data_directory>" >&2
	exit 1
fi

SCRIPT_PATH=$(dirname $(realpath -s "$0"))

if [ ! -d "$SCRIPT_PATH/data/$1" ]; then
	echo "Data directory: $1 does not exist" >&2
	exit 1
fi

# Extract UUID values by listing directories that were downloaded
# UUID values are then surrounded by quotes and comma-separated for JSON

uuid_values=$(for file in $SCRIPT_PATH/data/"$1"/*; do
  basename "$file"
done | sed "$ ! s/\(.*\)/\"\1\",/; $ s/\(.*\)/\"\1\"/")

# Json payload:

payload='{
	"filters":{
		"op":"in",
		"content":{
			"field":"files.file_id",
			"value":['$uuid_values']
		}
	},
	"format":"CSV",
	"fields":"file_id,file_name,cases.submitter_id,cases.samples.sample_type",
	"size":"100000"
}'

# Define what the output directory will be

output_dir="$SCRIPT_PATH/organized/""$1"
input_dir="$SCRIPT_PATH/data/$1"

# Retrieve metadata information from GDC

data=$(curl --request POST --header "Content-Type: application/json" --data "$payload" 'https://api.gdc.cancer.gov/files')
data=`echo "$data" | tr -d '\r'`

file_name_col=$(echo "$data"|head -1|tr "," "\n"|grep -n file_name|cut -d":" -f1)
file_id_col=$(echo "$data"|head -1|tr "," "\n"|grep -n file_id|cut -d":" -f1)
tcga_id_col=$(echo "$data"|head -1|tr "," "\n"|grep -n "submitter_id"|cut -d":" -f1)
sample_type_col=$(echo "$data"|head -1|tr "," "\n"|grep -n "sample_type"|cut -d":" -f1)
data=$(echo "$data"|awk -F"," -v a=$file_name_col -v b=$tcga_id_col -v c=$file_id_col -v d=$sample_type_col '{print($a","$b","$c","$d)}')

# Unzips downloaded flies and puts them into their appropriate directories by patient ID

echo "$data" | awk -v i="$input_dir" -v o="$output_dir" -F',' 'BEGIN { RS="\n" } NR>1 {
dir = o "/" ""$2"/"$4
command = "mkdir -p \""dir"\""
if ($1 ~ /\.gz$/) {
    filename=substr($1, 1, length($1)-3)
    if ($1 ~ /\.htseq\.counts\.gz$/) {
        filename="counts.txt"
        unzip_command = "gunzip -c \"" i "/" ""$3"/"$1"\" > \""dir"/"filename"\""
        command = command" && "unzip_command
    }
} else {
    filename=$1
    if ($1 ~ /\.xml$/) {
        dir=o "/" ""$2
        if ($1 ~ /biospecimen/) {
            filename="biospecimen.xml"
            cp_command = "cp \"" i "/" ""$3"/"$1"\" \""dir"/"filename"\""
            command = command" && "cp_command
        }
    }
}
system(command)
}'

