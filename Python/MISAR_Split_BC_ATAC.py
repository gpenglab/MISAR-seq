from Bio.SeqIO.QualityIO import FastqGeneralIterator
from gzip import open as gzopen

import argparse

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required=True, help="input file")
ap.add_argument("-o1", "--output_R1", required=True, help="output file R1")
ap.add_argument("-o2", "--output_R2", required=True, help="output file R2")
args = vars(ap.parse_args())

input_file_R2 = args["input"]
output_file_R1 = args["output_R1"]
output_file_R2 = args["output_R2"]

seq_end=75

bc2_start=75
bc2_end=83

bc1_start=115
bc1_end=123

with gzopen(input_file_R2, "rt") as in_handle_R1, open(output_file_R1, "w") as out_handle_R1, open(output_file_R2, "w") as out_handle_R2:
    for title, seq, qual in FastqGeneralIterator(in_handle_R1):
        new_seq_R1 = seq[:seq_end]
        new_qual_R1 = qual[:seq_end]
        barcode = seq[bc2_start:bc2_end] + seq[bc1_start:bc1_end] # !!! BC2 + BC1
        new_qual_R2 = qual[bc2_start:bc2_end] + qual[bc1_start:bc1_end]        
        out_handle_R1.write("@%s\n%s\n+\n%s\n" % (title, new_seq_R1, new_qual_R1))
        out_handle_R2.write("@%s\n%s\n+\n%s\n" % (title, barcode, new_qual_R2))
