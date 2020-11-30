#!/usr/bin/env python3
import os
import argparse
from collections import Counter
from pandas import DataFrame

args_parser = argparse.ArgumentParser(
    prog="SELEXderep",
    description="SELEXderep combines FASTA files acquired from a SELEX process. "
                "The FASTA files should be preprocessed (filtered, merged and primers cut). "
                "SELEXderep writes to STDOUT or writes to FASTA files with the primerless sequences as identifier "
                "and the corresponding read counts per SELEX round.\n"
                "\n"
                "Author: Ulrich Aschl (ulrich.aschl@tuwien.ac.at)",
    formatter_class=argparse.RawTextHelpFormatter
)
args_parser.add_argument("-o", "--out-fasta", type=str, default=None,
                         help="Output fasta file location. Default: stdout")
args_parser.add_argument("-f", "--force", action="store_true", default=False,
                         help="Overwrite output files if they already exists")
args_parser.add_argument("-c", "--out-csv", type=str, default=None,
                         help="Output csv file location. Default: no csv output")
args_parser.add_argument("SELEX_Rounds", metavar="SELEX_Round_FASTA", type=str, nargs='+',
                         help="FASTA files from SELEX process. Files have to be in sequential order")


def get_round_name(fasta_file_path):
    return os.path.splitext(os.path.basename(fasta_file_path))[0]


def read_file(fasta_file_path, seq_hash_map: dict):
    fasta_file = open(fasta_file_path, "r")

    round_hash_list = []

    for seq_id in fasta_file:
        seq = fasta_file.readline().rstrip('\n')
        seq_hash = hash(seq)
        if seq_hash_map.get(seq_hash, 0) == 0:
            seq_hash_map[seq_hash] = seq

        round_hash_list.append(seq_hash)
    return round_hash_list


def read_files(fasta_files):
    seq_hash_map = dict()
    selex_dict = dict()
    # Reading of sequences
    for fasta_file in fasta_files:
        round_name = get_round_name(fasta_file)
        selex_dict[round_name] = Counter(read_file(fasta_file, seq_hash_map))

    counts_df = DataFrame.from_dict(selex_dict)
    index_df = DataFrame.from_dict(seq_hash_map, orient='index', columns=["seq"])
    selex_df = counts_df.join(index_df)
    selex_df = selex_df.set_index("seq")
    selex_df[selex_df.isna()] = 0
    return selex_df


def print_fasta_file(out_file, round_names, sequences):
    # Print counted sequences to out_file
    for sequence, round_counts in sequences.iterrows():
        round_counts = '-'.join(str(int(x)) for x in list(round_counts))

        # WRITE TO FILE:
        out_file.write(">{} {}\n".format(sequence, round_counts))
        out_file.write("{}\n".format(sequence))

    out_file.flush()


def print_csv_file(out_file, round_names, sequences):
    # CSV Header
    out_file.write("sequence\t{}\n".format('\t'.join(round_names)))
    for sequence, round_counts in sequences.iterrows():
        # WRITE TO FILE:
        out_file.write("{}\t{}\n".format(
            sequence,
            '\t'.join(str(int(x)) for x in list(round_counts))
        ))
    out_file.flush()


def main():
    args = args_parser.parse_args()

    # ======= Parameter verification =======
    out_fasta = args.out_fasta
    if (out_fasta is not None) and (not args.force) and os.path.exists(out_fasta):
        print("Output fasta file {} already exists. Add --force to overwrite the output file".format(
            out_fasta))
        print("Halting.")
        return 1

    out_csv = args.out_csv
    if (out_csv is not None) and (not args.force) and os.path.exists(out_csv):
        print("Did not start! Output csv file {} already exists. Add --force to overwrite the output file".format(
            out_csv))
        print("Halting.")
        return 1

    fasta_files = args.SELEX_Rounds
    # check if fasta files exist
    missing_file_flag = False
    for fasta_file in fasta_files:
        if not os.path.exists(fasta_file):
            print("Input file {} does not seem to be at the specified path.".format(fasta_file))
            missing_file_flag = True
    if missing_file_flag:
        print("Halting.")
        return 1

    # ======= Read fasta files =======
    sequences = read_files(fasta_files)

    # ======= Print fasta file =======
    round_names = [get_round_name(fasta_file) for fasta_file in fasta_files]
    if out_fasta is not None:
        out_file = open(out_fasta, "w")
        print_fasta_file(out_file, round_names, sequences)
        out_file.close()
    else:
        print_fasta_file(sys.stdout, round_names, sequences)

    # ======= Print csv file =======
    if out_csv is not None:
        out_file = open(out_csv, "w")
        print_csv_file(out_file, round_names, sequences)
        out_file.close()


main()
