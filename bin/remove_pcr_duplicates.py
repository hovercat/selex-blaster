#!/usr/bin/env python3

import argparse

from Sequence_Node import Sequence_Node

args_parser = argparse.ArgumentParser()
args_parser.add_argument("-f", type=str, required=True)
args_parser.add_argument("-b", type=str, required=True)
args_parser.add_argument("-o", type=str, required=True)

import pandas as pd
import numpy as np
import networkx as nx


def main():
    args = args_parser.parse_args()

    print("Reading Blast file")
    blast = pd.read_csv(args.b,
                        delimiter="\t",
                        header=0,
                        names=[
                            "seq_x",
                            "seq_y",
                            "identity",
                            "length",
                            "mismatches",
                            "gaps",
                            "start_x",
                            "end_x",
                            "start_y",
                            "end_y",
                            "e_val",
                            "bitscore"
                        ])


    # read_sequence_file
    print("Reading FASTA file")
    g = nx.Graph()
    sequence_nodes = dict()
    friends_list = list()
    remove_nodes = list()
    with open(args.f, "r") as fasta_file:
        for id in fasta_file:
            id = id.rstrip()
            seq = fasta_file.readline().rstrip()
            seq_node = Sequence_Node(seq, id)

           # if seq_node.total_count >= 2:
           #     sequence_nodes[seq] = seq_node
           # else:
           #     remove_nodes.append(seq_node)
            sequence_nodes[seq] = seq_node

    #filtered_results = blast[blast.e_val <= 0.05]
    filtered_results = blast
    filtered_results = filtered_results[filtered_results.seq_x.isin(sequence_nodes.keys()) & filtered_results.seq_y.isin(sequence_nodes.keys())]
    seq_x = filtered_results['seq_x'].to_list()
    seq_y = filtered_results['seq_y'].to_list()
    del filtered_results
    del blast

    g.add_edges_from(zip(seq_x, seq_y))

    print("Finding best friends and print to fasta")
    representative_sequences = []
    with open(args.o, 'w') as out:
        for component in nx.connected_components(g):
            if len(component) == 1:
                seq = component.pop()
                seq_node = sequence_nodes[seq]

                out.write('>{} {}\n'.format(seq_node.seq, '-'.join([str(x) for x in seq_node.count])))
                out.write('{}\n'.format(seq_node.seq))
                out.flush()
            else:
                seq_nodes = [sequence_nodes[x] for x in component]
                max_count = max(seq_node.total_count for seq_node in seq_nodes)
                indices = [i for i, x in enumerate(seq_nodes) if x.total_count == max_count]
                counts = np.sum([sequence_nodes[x].count for x in component], axis=0)

                best_node = seq_nodes[indices[0]]
                out.write('>{} {}\n'.format(best_node.seq, '-'.join([str(x) for x in counts])))
                out.write('{}\n'.format(best_node.seq))
                out.flush()



main()
