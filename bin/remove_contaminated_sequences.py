#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("-cont", type=str)
parser.add_argument("-lib", type=str)
parser.add_argument("-o", type=str)
args = parser.parse_args()
#args = parser.parse_args(["-cont", "../work/30/fac29bc1d78b96617032822ecdf4ef/primer_contamination.csv",
#                          "-lib", "../work/30/fac29bc1d78b96617032822ecdf4ef/derepped.mea.fasta",
#                          "-o", "derp"])


contaminated_seqs = []
with open(args.cont, "r") as contamination_file:
    for line in contamination_file:
        line = line.rstrip().split('\t')
        contaminated_seqs.append(line[1])


with open(args.lib, "r") as library_file, open(args.o, "w") as out_file:

    while True:
        seq_id = library_file.readline().rstrip()
        seq = library_file.readline().rstrip()

        if seq_id is None or seq_id == "":
            break

        if seq.upper() in contaminated_seqs:
            continue

        out_file.write(seq_id)
        out_file.write("\n")
        out_file.write(seq)
        out_file.write("\n")
