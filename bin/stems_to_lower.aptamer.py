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

while True:
	line = fold_file.readline().strip()
	if line is None or line == '':
		break
	
	if line.startswith('>'):
		seq_name = line[1:]
		seq_seq = fold_file.readline().strip()
		j = 0
	else:
		structure = line.strip()
		seq_lc = []
		for i in range(0, len(seq_seq)):
		    if i < start_length or i > len(seq_seq) - end_length:
		        seq_lc.append(seq_seq[i].lower())
		        continue
		    if structure[i] == '(' or structure[i] == ')':
		        seq_lc.append(seq_seq[i].lower())
		    else:
		        seq_lc.append(seq_seq[i])

		lowercase_fasta.write(
		    ">{}_{}\n".format(j, seq_name)  # NAME
		)

		lowercase_fasta.write(
		    "{}\n".format(''.join(seq_lc)) # LC SEQ
		)
		
		j+=1
fold_file.close()
lowercase_fasta.close()
