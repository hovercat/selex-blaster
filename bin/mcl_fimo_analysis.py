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

tags_cluster = set()
def analyse_cluster(filename, selex, tags):
    'converts a filename to a pandas dataframe'
    start = time.time()
    try:
        #cluster = pd.read_csv(filename, delimiter='\t', comment='#')
        cluster_xml = etree.parse(filename)
        cluster_dict = {}

        for sequence in cluster:
            cluster.append()
    except Exception:
        return None

    cluster_name = filename.split('/')[-2]

    #cluster = cluster[cluster['score'] > 0]


    #tags_cluster = set()
    # look for tag
    cluster_tags = cluster[["motif_id", "sequence_name"]].merge(tags, how="inner", on="sequence_name")
    for index, row in cluster_tags.iterrows():
        print("jaa")
        tag = row[3]
        motif_id = row[0]
        seq = row[1]
        tags_cluster.add((seq, tag, cluster_name, motif_id))

    cluster_selex = cluster.merge(selex, on="sequence_name")

    cluster_data = cluster_selex.groupby(["motif_id", "round_id"]).agg(
        count_total=('count', sum),
        n=('count', lambda x: sum(x > 0)),
        med_fimo_score=('score', 'mean'),
        sd_fimo_score=('score', 'std')
    )

    # e.g. from '/home/ulrichaschl/Bioinformatics/Diplomarbeit/github/selex-blaster/output_EF05.excl.nomin.all/fimo/c000003/fimo.tsv'
    # we get c000003 as value
    cluster_data['cluster'] = cluster_name

    print("{} read: {}s".format(cluster_name, time.time() - start))

    return cluster_data.reset_index()


def main():
    # read aptamer dereplicated csv
    selex = pd.read_csv(
        '/home/ulrichaschl/Bioinformatics/Diplomarbeit/github/selex-blaster/output_EF05.excl.nomin.all/selex.aptamers.raw.csv',
        delimiter='\t', comment='#')
    selex = pd.wide_to_long(selex, stubnames="R", i="sequence", j="round_id")
    selex = selex.reset_index()
    selex.columns = ["sequence_name", "round_id", "count"]

    # read tag file
    tags = pd.read_csv('/home/ulrichaschl/Bioinformatics/Diplomarbeit/github/selex-blaster/tagfile.txt', delimiter='\t',
                       names=["sequence_name", "tag"])
    tags = tags.reset_index()

    # get a list of file names
    files = glob.glob(
        '/home/ulrichaschl/Bioinformatics/Diplomarbeit/github/selex-blaster/output_EF05.excl.nomin.all/fimo/c*/mast.xml')
    files = sorted(files)[1:300]

    print("start")
    if multiproc:
        # set up your pool
        with Pool(processes=20) as pool:  # or whatever your hardware can support

            # have your pool map the file names to dataframes
            df_list = pool.starmap(analyse_cluster, zip(files, repeat(selex), repeat(tags)))

            # reduce the list of dataframes to a single dataframe
            combined_df = pd.concat(df_list, ignore_index=True)
    else:
        combined_df = None
        for cluster in files:
            cluster_data = analyse_cluster(cluster, selex, tags)
            if combined_df is None:
                combined_df = cluster_data
            else:
                combined_df = combined_df.append(cluster_data)

    avg_round_counts = combined_df.merge(selex.groupby("round_id").agg({"count": "sum"}), on="round_id")
    avg_round_counts["count_p"] = avg_round_counts["count_total"] / avg_round_counts["count"]
    avg_round_counts["counts_per_seq"] = avg_round_counts["count_total"] / avg_round_counts["n"]

    avg_round_counts_p = avg_round_counts.reset_index().pivot(index=["cluster", "motif_id"], columns="round_id",
                                                              values="count_p")
    avg_round_counts_p.to_csv(
        "/home/ulrichaschl/Bioinformatics/Diplomarbeit/github/selex-blaster/output_EF05.excl.nomin.all/motif_analysis.csv")

    avg_round_counts_cps = avg_round_counts.reset_index().pivot(index=["cluster", "motif_id"], columns="round_id",
                                                                values="counts_per_seq")
    avg_round_counts_cps.to_csv(
        "/home/ulrichaschl/Bioinformatics/Diplomarbeit/github/selex-blaster/output_EF05.excl.nomin.all/motif_analysis.cps.csv")

    avg_round_counts.to_csv(
        "/home/ulrichaschl/Bioinformatics/Diplomarbeit/github/selex-blaster/output_EF05.excl.nomin.all/motif_analysis.long.csv")

    with open(
            "/home/ulrichaschl/Bioinformatics/Diplomarbeit/github/selex-blaster/output_EF05.excl.nomin.all/motif_analysis.tags.csv",
            "w") as motif_tags:
        motif_tags.write('\t'.join(["tag", "seq", "cluster_id", "motif_id"]))
        motif_tags.write('\n')
        #for index, clusters in tags_cluster.items():
        for row in tags_cluster:
            seq = row[0]
            tag = row[1]
            cluster = row[2]
            motif_id = row[3]
            #seq = index[0]
            #tag = index[1]
            #for cluster in clusters:
            #    cluster_id = cluster[0]
            #    motif_id = cluster[1]
            motif_tags.write('\t'.join([tag, seq, cluster, motif_id]))
            motif_tags.write('\n')

    print("end")
    print(len(avg_round_counts))


multiproc = False

if __name__ == '__main__':
    main()
