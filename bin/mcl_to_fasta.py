#!/usr/bin/env python3
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--min", type=int, help="Minimum number of cluster sequences", required=True)
parser.add_argument("--mcl-file", type=str, help="MCL output file", required=True)
args = parser.parse_args()

mcl_file = args.mcl_file

os.mkdir("clusters")
with open(mcl_file, "r") as mcl_cluster_file:
    i = 0
    for cluster in mcl_cluster_file:
        cluster_size = 0
        for aptamer in cluster.split():
            cluster_size = cluster_size + 1
        if cluster_size < args.min:
            continue

        with open("clusters/c{:06d}.fasta".format(i), "w") as single_cluster_file:
            for aptamer in cluster.split():
                single_cluster_file.write(">{}\n".format(aptamer))
                single_cluster_file.write("{}\n".format(aptamer))
        i = i+1
