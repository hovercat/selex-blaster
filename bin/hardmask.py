#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, required=True)
args = parser.parse_args()

with open(args.i, 'r') as fasta_file:
    for seq_id in fasta_file:
        seq = fasta_file.readline().rstrip()

        seq_lc = []
        for c in seq:
            if c.islower():
                seq_lc.append('N')
            else:
                seq_lc.append(c)

        print(seq_id, end="")
        print("{}".format(''.join(seq_lc)))
