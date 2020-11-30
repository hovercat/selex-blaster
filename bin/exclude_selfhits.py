#!/usr/bin/env python3
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--blast-file", type=str, help="BLAST CSV file", required=True)
args = parser.parse_args()

with open(args.blast_file, "r") as blast_file:
    i = 0
    for line in blast_file:
        line_split = line.split('\t')

        if line_split[0] == line_split[1]:
            continue

        print(line.strip())

