#!/usr/bin/env python3
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("-i", type=str)
parser.add_argument("-l", type=int)
args = parser.parse_args()

seqs = set()

with open(args.i, "r") as input_fasta:
    while True:
        id_str = input_fasta.readline().strip()
        seq_str = input_fasta.readline().strip()

        if id_str is None or id_str == "":
            break

        seq_str = seq_str.replace('a', 'N').replace('g', 'N').replace('c', 'N').replace('t', 'N')

        if len(seq_str) == args.l and 'CCCCC' not in seq_str and seq_str not in seqs:
            print(id_str, end='\n')
            print(seq_str, end='\n')
            seqs.add(seq_str)

