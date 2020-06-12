# Make our "final" bed files
# Because some genes may appear multiple times due to having multiple transcription start sites (TSS's), this script will aggregate the duplicates
# Only the gene with the highest average enrichment for the signal in its TSS (+/- padding) will be used.

import pyBigWig
import csv
import sys

input_bigwig_file = sys.argv[1]
input_bed_file =  sys.argv[2]
output_file = sys.argv[3]
padding = int(sys.argv[4])

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
                if abs(avg) > abs(current_avg):
                        output[cols[3]] = [cols[0], cols[1], cols[2], avg, cols[4], cols[5]]
        else:
                output[cols[3]] = [cols[0], cols[1], cols[2], avg, cols[4], cols[5]]

bw.close()
with open(output_file, 'w') as outfile:
	csv_writer = csv.writer(outfile, delimiter='\t')
	for k,v in output.items():
		o = csv_writer.writerow([v[0], v[1], v[2], k, v[4], v[5]])
