Disclaimer: Do not use, as this does not work properly. The k-mer based graph generation currently leads to a much too fragmented graph to gain any real insights.

# NextFlow Pipeline selex-blaster
The NextFlow pipeline *selex-blaster* is performing motif-search and clustering on aptamer pools gathered in SELEX.

*selex-blaster* folds all aptamers using RNAfold, to find all regions which directly interact with another region of the aptamer itself.
These regions are then masked in the sequences, which are used to build a BLAST database.
The only regions unmasked are loops, which have been shown to interact with targets.

The aptamer pool is then searched for in the database.
Sequences with similar motifs are fetched and returned with an E-score, a bit score and an alignment length.
The result is a network of similarity scores between the aptamers, based on short (4-10nt) unfolded sequences.
<!-- todo e-score? -->
The BLAST results are then clustered using the MCL algorithm (Markov Cluster Algorithm).
DREME, a discriminative derivate of the MEME-algorithm, is then used on finding the motif sequence of every cluster.
The found motifs are searched again in the aptamer pool.
These found sequences are then used for enrichment analysis.

For further instructions please see:
'U. Aschl:
"Data Analysis of HT-SELEX against Complex Targets";
Betreuer/in(nen): A.H. Farnleitner, G. Reischer, C. Kolm; Institut für Verfahrenstechnik, Umwelttechnik und Technische Biowissenschaften, 2021; Abschlussprüfung: 22.11.2021.'
