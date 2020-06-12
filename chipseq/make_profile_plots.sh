#!/bin/sh

r="data/refGene.bed"
padding=3000

make_bed_file () {
    refbedfile="$1"
    fname="$2"
    tmpfile=$(mktemp)
    trap "rm -f $tmpfile" 0 2 3 15
    tmpfile_pattern=$(mktemp)
    trap "rm -f $tmpfile_pattern" 0 2 3 15
    genes_list=$(while read -r line
    do
        echo "$line""$(printf '\t')"
    done < "$fname")
    echo "$genes_list"|sed '/^[[:space:]]*$/d'|tr -d '\r' > "$tmpfile_pattern"
    grep -f "$tmpfile_pattern" "$refbedfile"|awk '!seen[$0]++' > "$tmpfile"
    echo "$tmpfile"
}

# Make temporary bed files containing for our lists of genes
eumyc_up_tmp=$(make_bed_file $r data/eumyc_up.txt)
eumyc_down_tmp=$(make_bed_file $r data/eumyc_down.txt)
hcc_up_tmp=$(make_bed_file $r data/hcc_up.txt)
hcc_down_tmp=$(make_bed_file $r data/hcc_down.txt)

eumyc_down=$(mktemp)
trap "rm -f $eumyc_down" 0 2 3 15
hcc_up=$(mktemp)
trap "rm -f $hcc_up" 0 2 3 15
hcc_down=$(mktemp)
trap "rm -f $hcc_down" 0 2 3 15

# Do the main deepTools analysis
do_analysis() {
    prefix="$1"
    bed_file_up="$2"
    bed_file_down="$3"
    c_file="$4"
    t_file="$5"
    fc_file="$6"
    padding="$7"

    # Make our "final" bed files
    # Because some genes may appear multiple times due to having multiple transcription start sites (TSS's), these final files will aggregate the duplicates
    # Only the gene with the highest average enrichment for the signal in its TSS (+/- padding) will be used.
    final_bed_file_up=$(mktemp)
    trap "rm -f $final_bed_file_up" 0 2 3 15
    final_bed_file_down=$(mktemp)
    trap "rm -f $final_bed_file_down" 0 2 3 15
    python aggregate_dups.py "$fc_file" "$bed_file_up" "$final_bed_file_up" $padding
    python aggregate_dups.py "$fc_file" "$bed_file_down" "$final_bed_file_down" $padding

    # Plot with deepTools
    computeMatrix reference-point --referencePoint TSS -S "$fc_file" -R "$final_bed_file_up" "$final_bed_file_down" --beforeRegionStartLength $padding --afterRegionStartLength $padding --skipZeros -o "./output/mat/${prefix}_log2FC.mat.gz"
    computeMatrix reference-point --referencePoint TSS -S "$c_file" -R "$final_bed_file_up" "$final_bed_file_down" --beforeRegionStartLength $padding --afterRegionStartLength $padding --skipZeros -o "./output/mat/${prefix}_C.mat.gz
    plotProfile -m "./output/mat/${prefix}_log2FC.mat.gz" -out "./output/figures/profile_${prefix}.pdf" --colors red blue --regionsLabel Upregulated Downregulated
    tmp_outfile_sorted=$(mktemp)
    trap "rm -f $tmp_outfile_sorted" 0 2 3 15
    plotHeatmap -m "./output/mat/${prefix}_C.mat.gz" -out "./output/figures/heatmap_${prefix}_C.pdf" --colorMap Reds --whatToShow "heatmap and colorbar" --zMax 10 --outFileSortedRegions "$tmp_outfile_sorted"
    computeMatrix reference-point --referencePoint TSS -S "$t_file" -R "$tmp_outfile_sorted" --beforeRegionStartLength $padding --afterRegionStartLength $padding --skipZeros -o "./output/mat/${prefix}_T.mat.gz
    plotHeatmap -m "./output/mat/${prefix}_T.mat.gz" -out "./output/figures/heatmap_${prefix}_T.pdf" --colorMap Reds --whatToShow "heatmap and colorbar" --zMax 10 --sortRegions no
}

do_analysis eumyc_h3k27ac "$eumyc_up_tmp" "$eumyc_down_tmp" data/eumyc_h3k27ac_C.fc.signal.bigwig data/eumyc_h3k27ac_T.fc.signal.bigwig data/log2ratio_eumyc_h3k27ac.bigwig $padding
do_analysis eumyc_h3k4me3 "$eumyc_up_tmp" "$eumyc_down_tmp" data/eumyc_h3k4me3_C.fc.signal.bigwig data/eumyc_h3k4me3_T.fc.signal.bigwig data/log2ratio_eumyc_h3k4me3.bigwig $padding
do_analysis hcc_h3k27ac "$hcc_up_tmp" "$hcc_down_tmp" data/hcc_h3k27ac_C.fc.signal.bigwig data/hcc_h3k27ac_T.fc.signal.bigwig data/log2ratio_hcc_h3k27ac.bigwig $padding
do_analysis hcc_h3k4me3 "$hcc_up_tmp" "$hcc_down_tmp" data/hcc_h3k4me3_C.fc.signal.bigwig data/hcc_h3k4me3_T.fc.signal.bigwig data/log2ratio_hcc_h3k4me3.bigwig $padding
