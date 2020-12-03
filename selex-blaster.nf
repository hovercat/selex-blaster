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
    publishDir "${params.output_dir}",
        mode: "copy"
        
    input:
        file(fasta) from fasta_input
    output:
        file("randomregion.unique.fasta") into dereplicated_selex
        file("selex.aptamers.raw.csv") into dereplicated_csv
    script:
    """ 
        selex_dereplicate_fasta.py -o randomregion.unique.fasta -c selex.aptamers.raw.csv ${fasta}
    """
}
process draw_sample {
    publishDir "${params.output_dir}",
        mode: "copy"
    conda 'python=2 bioconda::meme'
    input:
        file("randomregion.unique.fasta") from dereplicated_selex
    output:
        tuple file("randomregion.unique.resmpl.fasta"), file("aptamer.unique.fasta") into dereplicated_selex_smpl
        tuple file("randomregion.unique.resmpl.fasta"), file("aptamer.unique.fasta") into dereplicated_selex_meme
    script:
    """ 
        fasta-subsample randomregion.unique.fasta ${params.sample} > randomregion.unique.resmpl.fasta
        if [ ${params.sample} -lt 1 ]
        then
            cp randomregion.unique.fasta randomregion.unique.resmpl.fasta
        fi
        sed '2~2s/^/${params.primer5}/;2~2s/\$/${params.primer3}/' randomregion.unique.resmpl.fasta > aptamer.unique.fasta
    """
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
Folding Aptamers
========================================================    
"""
process folding {
    publishDir "${params.output_dir}/fold"
    cpus params.cpus

    input:
        tuple file("randomregion.unique.fasta"), file("aptamer.unique.fasta") from dereplicated_selex_smpl
    output:
        tuple file("randomregion.unique.fasta"), file("aptamer.unique.fasta"), file("aptamer.unique.mfe_masked.fasta"), file("aptamer.unique.mea_masked.fasta") into folded_to_blastdb
        file("aptamer.unique.mea_masked.fasta") into folded_to_blastquery
        file("aptamer.unique.mea_masked.fasta") into folded_mast
    script:
    """
         RNAfold --paramFile=$params.RNAfold_mathews2004_dna \
            --noPS --noconv --noDP --partfunc=1 \
            --jobs=${params.cpus} \
            -T ${params.folding_temp} \
            --MEA \
            -i aptamer.unique.fasta \
            --outfile=aptamer.unique.vienna
            
        awk 'BEGIN { RS = ">"; FS = "\\n|( \\\\( *)|)\\n" } { print \$1, "\\t", \$2, "\\t", \$3, "\\t", \$4 }' aptamer.unique.vienna > aptamer.unique.vienna.mfe
        awk 'BEGIN { RS = ">"; FS = "\\n|(\\\\{(\\-|\\ |[0-9]+)|MEA=)" } { print \$1, "\\t", \$2, "\\t", \$7, "\\t", \$8 }' aptamer.unique.vienna > aptamer.unique.vienna.mea

        stems_to_lower.py aptamer.unique.vienna.mfe aptamer.unique.mfe_masked.fasta
        stems_to_lower.py aptamer.unique.vienna.mea aptamer.unique.mea_masked.fasta
    """
}


"""
========================================================
Making a BLAST db from masked sequences
========================================================
"""

process make_blast_db {
    conda 'bioconda::blast'
    input:
        tuple file("randomregion.unique.fasta"), file("aptamer.unique.fasta"), file("aptamer.unique.mfe_masked.fasta"), file("aptamer.unique.mea_masked.fasta") from folded_to_blastdb
    output:
        file("blast_db_decontaminated*") into blast_db
    script:
    """
        # Search for primer sequences contaminated with primers and remove them (TODO, necessary?)
        # Make BLAST db with all sequences
        convert2blastmask -in aptamer.unique.mea_masked.fasta -masking_algorithm repeat -masking_options "repeatmasker, default" -outfmt maskinfo_asn1_bin -out blast_db.mask
        makeblastdb -in aptamer.unique.mea_masked.fasta -dbtype nucl -mask_data blast_db.mask -out blast_db
        # Dummy primer fasta file
        touch primers.fasta
        echo '>${params.primer5}' >> primers.fasta
        echo '${params.primer5}' >> primers.fasta
        echo '>${params.primer3}' >> primers.fasta
        echo '${params.primer3}' >> primers.fasta
        # BLASTn search
        blastn -task blastn-short -db blast_db -query primers.fasta -num_threads 1 -evalue 1000 -strand plus -reward 1 -penalty -4 -gapopen 1 -gapextend 2 -word_size 6 -db_hard_mask 40  -outfmt 6 > primer_contamination.csv
        remove_contaminated_sequences.py -cont primer_contamination.csv -lib aptamer.unique.mea_masked.fasta -o aptamer.unique.mea_masked.decontaminated.fasta
        
        # Make BLAST db with decontaminated input
        convert2blastmask -in aptamer.unique.mea_masked.decontaminated.fasta -masking_algorithm repeat -masking_options "repeatmasker, default" -outfmt maskinfo_asn1_bin -out blast_db_decontaminated.mask
        makeblastdb -in aptamer.unique.mea_masked.decontaminated.fasta -dbtype nucl -mask_data blast_db_decontaminated.mask -out blast_db_decontaminated
    """
}



"""
========================================================
Searching the BLAST db for similarities
========================================================
"""
folded_to_blastquery
    .splitFasta(by: 100, file: "query.fasta")
    //.take(30)
    .set { folded_to_blastquery_split }

process blast_motifs {
    conda 'bioconda::blast'
    cpus 1
    maxForks params.cpus
    
    input:
        file(blast_db) from blast_db.collect()
        each file(f) from folded_to_blastquery_split
    output:
        file("blast.${f}.csv") into blast_motifs_results
    script:
    """
        blastn -task blastn-short -db blast_db_decontaminated -query $f \
            -num_threads 1 -strand plus \
            -reward 1 -penalty -4 -gapopen 1 -gapextend 2 \
            -outfmt 6 -max_target_seqs 100 -db_hard_mask 40 -evalue 10000 \
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
        find -L . -type f -name 'blast.*.csv' -print0 | xargs -0 cat >> collect.blast.csv
        exclude_selfhits.py --blast-file collect.blast.csv > collect.blast.exclselfhits.csv
        # cp collect.blast.csv collect.blast.exclselfhits.csv
        
        # min_val=-\$(cut -f 4 collect.blast.exclselfhits.csv | sort -n | head -1 | tr -d '[:space:]')    
        min_val=0
        cut -f 1,2,4 collect.blast.exclselfhits.csv > seq.abc
        mcxload -abc seq.abc --stream-mirror -stream-tf "add(\$min_val)" -o seq.mci -write-tab seq.tab
#        mcxload -abc seq.abc --stream-mirror  -o seq.mci -write-tab seq.tab
    """
}


process do_cluster_motifs {
    publishDir "${params.output_dir}/mcl",
        mode: "copy"
    conda 'bioconda::mcl'
    cpus params.cpus
    maxForks 1
    
    input:
        tuple file("seq.mci"), file("seq.abc"), file("seq.tab") from mcl_preparation
        val inflation from(1.4) // ,2, 2.6, 4)
    output:
        file("${inflation}/c*.fasta") into clusters_fasta
    script:
    """
        mcl seq.mci -I $inflation -te ${params.cpus} -o seq.mcl
        mcxdump -icl seq.mcl -tabr seq.tab -o dump.seq.mci
        mcl_to_fasta.py --min 1 --mcl-file dump.seq.mci
        
        # move from ./clusters to ./
        mkdir ${inflation}
        find -L clusters -type f -name 'c*.fasta' -print0 | xargs -r0 mv -t ${inflation}/.
    """
}


"""
========================================================
Finding Motifs in their clusters
========================================================
"""
library_ch = Channel.fromPath(params.input_dir + "/R0.fasta" , checkIfExists:true, type: "file").collect()
process extract_cluster_motifs {
    conda 'python=2 bioconda::meme'
    maxForks params.cpus
    publishDir "${params.output_dir}/dreme",
        pattern: "*dreme.html",
        mode: "copy"
    
    input:
        tuple file("randomregion.unique.fasta"), file("aptamer.unique.fasta") from dereplicated_selex_meme
        file(library) from library_ch
        each file(cluster) from clusters_fasta
    output:
        tuple file("${cluster.baseName}.dreme.html"), file("library_1000.fasta"), file(cluster), file("randomregion.unique.fasta"), file("aptamer.unique.fasta") into dreme
    script:
    """
        fasta-subsample $library 1000  > library_1000.fasta
        dreme -o dreme -norc -dna -p ${cluster} -n library_1000.fasta -mink 6
        mv dreme/dreme.html ${cluster.baseName}.dreme.html
    """
}

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
}


