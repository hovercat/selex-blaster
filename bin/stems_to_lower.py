#!/usr/bin/env python

"""
This script takes a folded sequence file in csv format as outputted by the
clustamer workflow and converts every stem element to lowercase.
This is done for the fasta search, to exclude non-loop areas.
Also the primer areas are lowercase masked.
"""
import sys

if len(sys.argv) < 2:
    print("Provide fold.csv and out.fasta.")
    sys.exit(1)

FOLD_FILE_PATH = sys.argv[1]
OUT_FILE_PATH = sys.argv[2]

start_length = 23
end_length = 23

lowercase_fasta = open(OUT_FILE_PATH, "w")
fold_file = open(FOLD_FILE_PATH, "r")
for line in fold_file:
    if len(line) < 50: # hacky workarond for AWK error
        continue

    line_spl = line.split('\t')

    #fold = line_spl[2].strip()[start:(start + length)]
    fold = line_spl[2].strip()
    #seq = line_spl[1].strip()[start:(start + length)]
    seq = line_spl[1].strip()

    seq_lc = []
    for i in range(0, len(seq)):
        if i < start_length or i > len(seq) - end_length:
            seq_lc.append(seq[i].lower())
            continue
        if fold[i] == '(' or fold[i] == ')':
            seq_lc.append(seq[i].lower())
        else:
            seq_lc.append(seq[i])

    lowercase_fasta.write(
        ">{}\n".format(line_spl[0])  # NAME
    )

    lowercase_fasta.write(
        "{}\n".format(''.join(seq_lc)) # LC SEQ
    )
fold_file.close()
lowercase_fasta.close()
