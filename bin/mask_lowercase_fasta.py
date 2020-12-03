#!/usr/bin/env python3
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("-i", type=str)
parser.add_argument("-c", action="store_true", help="cut primers")
args = parser.parse_args()
#args = parser.parse_args(["-i", "../output_EF05.excl.nomin.all/fold/aptamer.unique.mea_masked.fasta", "-c"])

with open(args.i, "r") as input_fasta:
    while True:
        id_str = input_fasta.readline().strip()
        seq_str = input_fasta.readline().strip()

        if id_str is None or id_str == "":
            break

        seq_str = seq_str.replace('a', 'N').replace('g', 'N').replace('c', 'N').replace('t', 'N')
        if args.c:
            seq_str = seq_str[23:len(seq_str)-23]

        print(id_str, end='\n')
        print(seq_str, end='\n')

