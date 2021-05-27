#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, required=True)
parser.add_argument("-p1", type=int, required=True)
parser.add_argument("-p2", type=int, required=True)
args = parser.parse_args()

with open(args.i, 'r') as fasta_file:
    for seq_id in fasta_file:
        seq = fasta_file.readline().rstrip()
        print(seq_id, end="")
        print(seq[args.p1:-args.p2])
