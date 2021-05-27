#!/usr/bin/env python3

import os
import time

import pandas as pd
from multiprocessing import Pool
from itertools import repeat

import glob
import xml.etree.ElementTree as etree

# wrap your csv importer in a function that can be mapped
from sphinx.addnodes import only
import re


tags_cluster = set()
def analyse_cluster(filename, selex, tags, out_file):
    start = time.time()
    cluster_name = filename.split('/')[-2]
    cluster_xml = etree.parse(filename)
    # motifs = {}
    # for motif in cluster_xml.find("motifs"):
    #     motif_str = motif.attrib['id']
    #     e_val = motif.attrib['evalue']
    #     motifs[motif_str] = []
    #
    #     for nt_pos in motif:
    #         motifs[motif_str].append({
    #             "A": nt_pos.attrib["A"],
    #             "C": nt_pos.attrib["C"],
    #             "G": nt_pos.attrib["G"],
    #             "T": nt_pos.attrib["T"],
    #         })

   # sequences = set()
   # mcl_fasta = filename[:filename.rindex('/')] + ".fasta"
   # with open(mcl_fasta) as mcl_fasta_file:
   #     while True:
   #         seq_id = mcl_fasta_file.readline().strip()
   #         seq_fold = mcl_fasta_file.readline().strip()
   #         seq_id_seq = seq_id.split(" ")[0][1:]
   #         if seq_id is None or seq_id == "":
   #             break
   #         sequences.add(seq_id_seq)

    cluster = {}
    for sequence in cluster_xml.find("sequences"):
        seq_str = sequence.attrib['name']
        combined_p_val = sequence[0].attrib['combined_pvalue']
        e_val = sequence[0].attrib['evalue']
    #for index, row in cluster.iterrows():
        out_file.write('\t'.join([seq_str, combined_p_val, e_val, cluster_name]))
        out_file.write('\n')


       # if seq_str in sequences:
        cluster[seq_str] = [combined_p_val, e_val, cluster_name]

    #cluster = pd.DataFrame(cluster, columns=["sequence_name", "combined_p_val", "ev"])
    for index, row in tags.iterrows():
        if row["sequence_name"] in cluster.keys():
            seq_data = cluster[row["sequence_name"]]
            tags_cluster.add((row["sequence_name"], row["tag"], seq_data[0], seq_data[1], seq_data[2]))

    # #tags_cluster = set()
    # # look for tag
    # if len(cluster) > 0:
    #     cluster_tags = cluster.merge(tags, how="inner", on="sequence_name")
    #     for index, row in cluster_tags.iterrows():
    #         tag = row["tag"]
    #         seq = row["sequence_name"]
    #         ev = row["ev"]
    #         p_val = row["combined_p_val"]
    #         tags_cluster.add((seq, tag, ev, p_val, cluster_name)) # todo this is not thread safe, but whatever
    #
    # cluster['cluster'] = cluster_name
   # print("{} read: {}s".format(cluster_name, time.time() - start))


    return cluster


def main():
    # read aptamer dereplicated csv
    selex = pd.read_csv(
        '/home/ulrichaschl/Bioinformatics/Diplomarbeit/github/selex-blaster/output_EF05.excl.nomin.all/selex.aptamers.raw.csv',
        delimiter='\t', comment='#')
    selex = pd.wide_to_long(selex, stubnames="R", i="sequence", j="round_id")
    selex = selex.reset_index()
    selex.columns = ["sequence_name", "round_id", "count"]

    # read folded selex file
    folded_sequences = []
  #  with open("/home/ulrichaschl/Bioinformatics/Diplomarbeit/github/selex-blaster/output_EF05.excl.nomin.all/fold/aptamer.unique.mea_masked.fasta", "r") as folded_file:
  #      while True:
  #          seq_id = folded_file.readline().strip()
  #          seq_fold = folded_file.readline().strip()
  #          seq_id_seq = seq_id.split(" ")[0][1:]
  #          if seq_id is None or seq_id == "":
  #              break
#
#            folded_sequences.append((seq_id_seq, seq_fold[23:len(seq_fold)-23]))

    #dt_folded_sequences = pd.DataFrame(folded_sequences)
    #dt_folded_sequences.columns = ["sequence_name", "folded"]

    #selex = selex.merge(dt_folded_sequences)

    # read tag file
    tags = pd.read_csv('/home/ulrichaschl/Bioinformatics/Diplomarbeit/github/selex-blaster/tagfile.txt', delimiter='\t',
                       names=["sequence_name", "tag"])
    tags = tags.reset_index()

    # get a list of file names
    files = glob.glob(
        '/home/ulrichaschl/Bioinformatics/Diplomarbeit/github/selex-blaster/output_EF05.excl.nomin.all/fimo/c*/mast.xml')
    files = sorted(files)#[0:400]

    print("start")

    combined_df = None
    with open("/home/ulrichaschl/Bioinformatics/Diplomarbeit/github/selex-blaster/output_EF05.excl.nomin.all/mast_cat.csv", "w") as out_file:

        out_file.write('\t'.join(["sequence_name", "combined_p_val", "ev", "cluster"]))
        out_file.write('\n')
        for cluster in files:
            cluster_data = analyse_cluster(cluster, selex, tags, out_file)

    #combined_df.to_csv("/home/ulrichaschl/Bioinformatics/Diplomarbeit/github/selex-blaster/output_EF05.excl.nomin.all/mast_cat.csv")

    with open(
            "/home/ulrichaschl/Bioinformatics/Diplomarbeit/github/selex-blaster/output_EF05.excl.nomin.all/motif_analysis.tags.csv",
            "w") as motif_tags:
        motif_tags.write('\t'.join(["tag", "seq", "ev", "combined_p_val", "cluster_id"]))
        motif_tags.write('\n')
        for row in tags_cluster:
            seq = row[0]
            tag = row[1]
            ev = row[2]
            combined_p_val = row[3]
            cluster = row[4]
            motif_tags.write('\t'.join([tag, seq, ev, combined_p_val, cluster]))
            motif_tags.write('\n')

    print("end")

multiproc = False

if __name__ == '__main__':
    main()
