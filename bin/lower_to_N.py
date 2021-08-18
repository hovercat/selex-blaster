#!/usr/bin/env python
import sys

if len(sys.argv) < 1:
    print("Provide fold.csv.")
    sys.exit(1)

with open(sys.argv[1], "r") as fold_file:
    for id in fold_file:
        id = id.strip()
        seq = fold_file.read().strip()

        seq = seq.replace('a', 'N').replace('g', 'N').replace('c', 'N').replace('t', 'N')

        print(id)
        print(seq)
