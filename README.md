Disclaimer: This pipeline ist still under construction and not working properly yet. Please come back in January 2020. :)

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
The found motifs are searched again in the aptamer pool. (I am not clear if this is good or bad)
These found sequences are then used for enrichment analysis.

## Output

The workflow generates a vast number of different files.
## Directory Structure
The workflow is embedded in a directory structure.
By default new directories are created containing resulting files and a cache.

    .
    ├── bin/                # R and Python scripts for analysis
    ├── output_YOUREXPERIMENT/             
    │   ├── dreme/     # Plots and charts on aptamers enrichment
    │   ├── fimo/
    │   ├── mcl/
    │   ├── randomregions.derep.csv
    │   ├── randomregions.derep.fasta
    │   ├── randomregions.derep.resample.fasta
    │   ├── aptamers.derep.fasta
    │   └── selex.aptamers.rpm.csv      # CSV chart containing all aptamers in the experiment, counted in reads per million

### Output 1


## Usage
 
### Execution
Before executing the workflow a config file has to be created.
This config file stores details about the SELEX experiment and where the workflow can find the FASTQ-files.

To create such a config file, you can either write it on your own, or use the **create_config.py** script, which is recommended.
```bash
python create_config.py
# or
./create_config.py
```
The script **create_config.py** will create a new file, ending with '.config' named after your experiment.
So, if you experiment name is 'SELEX 24.Nov.2020', the config will be called 'SELEX24.Nov.2020.config'.

Below is the command which can be used to execute the workflow.
```bash
nextflow run selex-ngs-prep.nf -c YOUR_SELEX.config
```

If you are working on a cluster you probably don't want NextFlow to place intermediate files in the working directory.
The working directory can be changed as seen below:
```bash
nextflow run selex-ngs-prep.nf -c YOUR_SELEX.config -w /scratch/xyz/work
```

### Default Config


## Software Used

- Some R-Packages by Hadley Wickham: ggplot2, tidyr, dplyr, stringr
- benjjneb's dada2
- More R-Packages: here, BiocManager, rmarkdown, knitr, latex2exp
- Python Packages: pandas
- NextFlow
- cutadapt
- conda with bioconda (?)
- fastp


## License

None here yet. 
Please cite my thesis if you use this in academia.
Please get in touch if you are a commercial user.
