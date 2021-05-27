#!/usr/bin/env nextflow
"""
========================================================
Groovy Helper Functions
========================================================    
"""
def reverse_complement(String s) {
    complement(s.reverse());
}

def complement(String s) {
    def acgt_map = [
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "a": "t",
        "c": "g",
        "g": "c",
        "t": "a"
    ];

    char[] sc = new char[s.length()];
    for (int i = 0; i < s.length(); i++) {
        sc[i] = acgt_map[s[i]];
    }
    new String(sc);
}

def remove_all_extensions(String s) {
    s.substring(0, s.indexOf("."));
}

def get_round_id(String round_name) {
    if (params.round_order == null || params.round_order == "" || params.round_order.size() == 0) return 0;
    else return params.round_order.indexOf(round_name);
}

    

"""
========================================================
Make output directory if it doesn't exist
========================================================    
"""
dir_output = file(params.output_dir)
if (!dir_output.exists()) {
    if (!dir_output.mkdir()) {
        println("Couldn't create output directory.")
    }
}



"""
========================================================
Data import
========================================================
"""    
fasta_input = Channel.fromPath(params.input_dir + "/" + params.fasta_pattern, checkIfExists:true, type: "file")
    .view()
    .map { it -> [get_round_id(it.baseName), it].flatten() }
    .filter { it[0] >= 0 }
    .toSortedList( { a -> a[0] } )
    .transpose()
    .last()
    .collect()
    
process dereplicate_sequences {
	conda 'pandas'
    publishDir "${params.output_dir}",
        mode: "copy"
        
    input:
        file(fasta) from fasta_input
    output:
        file("randomregion.unique.fasta") into dereplicated_selex_db, dereplicated_selex_q, dereplicated_selex_asv
    script:
    """ 
        selex_dereplicate_fasta.py -o randomregion.unique.fasta -c selex.aptamers.raw.csv ${fasta}
        
        head -n 100000 randomregion.unique.fasta > randomregion.unique.fasta_
        mv randomregion.unique.fasta_ randomregion.unique.fasta
    """
}
library_ch = Channel.fromPath(params.input_dir + "/" + params.library + params.fasta_pattern , checkIfExists:true, type: "file").collect()

"""
=========
ASV Removal
==========
"""
process asv_removal_blastdb {
	conda 'bioconda::blast'

    input:
		file(derep_fasta) from dereplicated_selex_db
    output:
        file("blast_db*") into asv_blastdb
    script:
    """
    	makeblastdb -in $derep_fasta -dbtype nucl -out blast_db
    """
}

dereplicated_selex_q
	.splitFasta(by: 1000, file:true)
	.set { dereplicated_selex_q_split }

process asv_removal_blastn {
	conda 'bioconda::blast'
	publishDir "${params.output_dir}/",
	mode: "copy"

    input:
		file(blast_db) from asv_blastdb.collect()
		each file(f) from dereplicated_selex_q_split
    output:
        file("blast.${f}.csv") into asv_removal_blastn_results
    script:
    """
    	blastn -task blastn-short -query $f -num_threads 1 -evalue 0.05  -strand plus -db blast_db -outfmt 6 -gapopen 0 -gapextend 4 -penalty -3 -reward 2 > blast.${f}.csv 
    """
}

process asv_concat_blastn_results {
	conda 'bioconda::blast'
	publishDir "${params.output_dir}/",
	mode: "copy"

    input:
		file(blast_csv) from asv_removal_blastn_results.collect()
    output:
        file("blastn_results.csv") into asv_removal_blastn_results_concat
    script:
    """
		find -L . -wholename './blast.*.csv' | sort | xargs cat > blastn_results.csv
    """
}




process asv_extraction {
	conda 'pandas networkx'
	publishDir "${params.output_dir}/",
	mode: "copy"

    echo true
    input:
    	file(derep) from dereplicated_selex_asv
		file(blastn_csv) from asv_removal_blastn_results_concat
    output:
         file("aptamers.clean.fasta") into dereplicated_asv_removed
    script:
    """
		remove_pcr_duplicates.py -f $derep -b $blastn_csv -o aptamers.clean.fasta
    """
}




"""
========================================================
Folding Aptamers
========================================================    
"""


process folding {
    publishDir "${params.output_dir}/fold"
    cpus params.cpus

    input:
	    file(random_region_fasta) from dereplicated_asv_removed
    output:
        file("aptamers.mfe.masked.fasta") into folded, folded_blast_query
    script:
    """
         sed '2~2s/^/${params.primer5}/;2~2s/\$/${params.primer3}/' ${random_region_fasta} > aptamers.fasta
         
         RNAfold --paramFile=$params.RNAfold_mathews2004_dna \
            --noconv --noPS --noDP --partfunc=1 \
            --jobs=${params.cpus} \
            -T ${params.folding_temp} \
            -i aptamers.fasta \
            --MEA \
            --outfile=aptamers.vienna
            
		awk 'BEGIN { RS = ">"; FS = "\\n|( \\\\( *)|)\\n" } { print \$1, "\\t", \$2, "\\t", \$3, "\\t", \$4 }' aptamers.vienna > aptamers.vienna.mfe
       
        stems_to_lower.py aptamers.vienna.mfe aptamers.mfe.masked.fasta
    """
}

"""
========================================================
BLAST 
========================================================
"""
process makeblastdb {
    conda 'bioconda::blast'
    input:
        file(aptamers_masked_fasta) from folded
    output:
        file("blast_db*") into blast_db
    script:
    """
		 convert2blastmask -in ${aptamers_masked_fasta} -masking_algorithm repeat -masking_options "repeatmasker, default" -outfmt maskinfo_asn1_bin -out blast_db.mask
         makeblastdb -in ${aptamers_masked_fasta} -dbtype nucl -mask_data blast_db.mask -out blast_db
    """
}

folded_blast_query
    .splitFasta(by: 1000, file: "query.fasta")
    .set { folded_blast_query_split }

process blast_motifs {
    conda 'bioconda::blast'
    cpus 1
    maxForks params.cpus
    
    input:
        file(blast_db) from blast_db.collect()
        each file(f) from folded_blast_query_split
    output:
        file("blast.${f}.csv") into blast_motifs_results
    script:
    """
        blastn -task blastn-short -db blast_db -query $f \
            -num_threads 1 -strand plus \
            -reward 1 -penalty -4 -gapopen 1 -gapextend 2 \
            -outfmt 6 -max_target_seqs 200 -db_hard_mask 40 -evalue 10000 \
            -lcase_masking > blast.${f}.csv
    """
}


"""
========================================================
Clustering of resulting network
========================================================
"""
process prepare_cluster_motifs {
    publishDir "debug/prepare_cluster_motifs"
    
    conda 'bioconda::mcl'
    input:
        file(f) from blast_motifs_results.collect()
    output:
        tuple file("seq.mci"), file("seq.abc"), file("seq.tab") into mcl_preparation
    script:
    """
    	# combine all blast results
        find -L . -type f -name 'blast.*.csv' -print0 | xargs -0 cat >> collect.all_blast.csv
        exclude_selfhits.py --blast-file collect.all_blast.csv > collect.blast.csv
        
        cut -f 1,2,4 collect.blast.csv > seq.abc
        mcxload -abc seq.abc --stream-mirror  -o seq.mci -write-tab seq.tab
    """
}


process do_cluster_motifs {
    publishDir "${params.output_dir}/I${params.inflation}.MCL/",
        mode: "copy"
    conda 'bioconda::mcl'
    cpus params.cpus
    maxForks 1
    
    input:
        tuple file("seq.mci"), file("seq.abc"), file("seq.tab") from mcl_preparation
    output:
        file("c*.fasta") into clusters_fasta
    script:
    """
        mcl seq.mci -I ${params.inflation} -te ${params.cpus} -o seq.mcl
        mcxdump -icl seq.mcl -tabr seq.tab -o dump.seq.mci
        mcl_to_fasta.py --out-dir . --min 1 --mcl-file dump.seq.mci
        
        # move from ./clusters to ./
        # find -L clusters -type f -name 'c*.fasta' -print0 | xargs -r0 mv -t ${params.inflation}/.
    """
}


"""
========================================================
Finding Motifs in their clusters
========================================================
"""
process fold_library {
	maxForks 1
	cpus params.cpus
	
	input:
		file(library) from library_ch
	output:
		file("library.mfe.masked.fasta") into library_folded
	
	script:
	"""
		sed '2~2s/^/${params.primer5}/;2~2s/\$/${params.primer3}/' ${library} > library.fasta
         
         RNAfold --paramFile=$params.RNAfold_mathews2004_dna \
            --noconv --noPS --noDP --partfunc=1 \
            --jobs=${params.cpus} \
            -T ${params.folding_temp} \
            -i library.fasta \
            --MEA \
            --outfile=library.vienna
            
		awk 'BEGIN { RS = ">"; FS = "\\n|( \\\\( *)|)\\n" } { print \$1, "\\t", \$2, "\\t", \$3, "\\t", \$4 }' library.vienna > library.vienna.mfe
       
        stems_to_lower.py library.vienna.mfe library.mfe.masked.fasta
	"""
}

process extract_cluster_motifs {
    conda 'bioconda::meme'
    maxForks params.cpus
    publishDir "${params.output_dir}/streme",
        pattern: "*streme.html",
        mode: "copy"
    
    input:
        file(library) from library_folded
        each file(cluster) from clusters_fasta
    output:
        tuple file("${cluster.baseName}.streme.html"), file("library_10000.fasta"), file(cluster) into streme
    script:
    """
       echo '
ALPHABET "DNA" RNA-LIKE

# Core symbols
A "Adenine" CC0000
C "Cytosine" 0000CC
T "Thymine" 008000
G "Guanine" FFB300
U "Uracil" 008000

# Ambiguous symbols
#U = T # alias Uracil to Thymine (permit U in input sequences)
R = AG
Y = CT
#K = GT
#M = AC
#S = CG
#W = AT
B = CGT
D = GAT
H = ACT
V = ACG
N = ACGT # wildcard symbol
#X = ACGT # wildcard symbol
' > alphabet.file
    
    
    	hardmask.py -i ${cluster} > cluster.masked.fasta
    	hardmask.py -i ${library} > library.masked.fasta
    	cutprimers.py -i library.masked.fasta -p1 ${params.primer5.length()} -p2 ${params.primer3.length()}  > library.rr.masked.fasta
    
        fasta-subsample library.rr.masked.fasta 10000  > library_10000.fasta
        streme --alph alphabet.file --minw 5 --maxw 16 --o streme --p cluster.masked.fasta --n library_10000.fasta 
        mv streme/streme.html ${cluster.baseName}.streme.html
    """
}


/*
"""
fimo

"""
process fimo_ {
    input:
        file("folded.mfe.fasta") from folded_mast
    output:
        file("hardmasked.mfe.fasta") into folded_mast_
    script:
    """
        mask_lowercase_fasta.py -i folded.mfe.fasta -c > hardmasked.mfe.fasta
    """
}

process fimo {
    errorStrategy 'ignore'
    conda 'python=2 bioconda::meme'
    maxForks params.cpus
    publishDir "${params.output_dir}/fimo/",
        mode: "copy"
    
    input:
        file("hardmasked.mfe.fasta") from folded_mast_
        tuple file(dreme), file(control_fasta), file(cluster), file("randomregion.unique.fasta"), file("aptamer.unique.fasta") from dreme
    output:
        file("${cluster.baseName}/") into fimo
        file(cluster) into fimo_fasta
    script:
    """
        #fimo --thresh 0.0002 --norc -o ${cluster.baseName} --bfile --motif-- $dreme randomregion.unique.fasta
        mast -norc -bfile --uniform-- -remcorr -c 5 -ev 100000 -o ${cluster.baseName} $dreme hardmasked.mfe.fasta
        #rm ${cluster.baseName}/*.xml
        #rm ${cluster.baseName}/*.gff
    """
}*/


