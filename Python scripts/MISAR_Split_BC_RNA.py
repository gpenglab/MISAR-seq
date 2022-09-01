from Bio.SeqIO.QualityIO import FastqGeneralIterator
from gzip import open as gzopen

import argparse

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required=True, help="input file")
ap.add_argument("-o", "--output", required=True, help="output file")

args = vars(ap.parse_args())

input_file = args["input"]
output_file = args["output"]

seq_end=0

umi_start=0
umi_end=10

bc2_start=10
bc2_end=18

bc1_start=50
bc1_end=58

with gzopen(input_file, "rt") as in_handle, open(output_file, "w") as out_handle:
    for title, seq, qual in FastqGeneralIterator(in_handle):
        barcode = seq[bc2_start:bc2_end] + seq[bc1_start:bc1_end] + seq[umi_start:umi_end] # !!! BC2 + BC1 + UMI
        new_qual = qual[bc2_start:bc2_end] + qual[bc1_start:bc1_end] + qual[umi_start:umi_end]       
        out_handle.write("@%s\n%s\n+\n%s\n" % (title, barcode, new_qual))
