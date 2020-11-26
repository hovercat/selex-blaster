# selex-ngs-prep
*Selex-ngs-prep* is a NextFlow workflow made for data preparation and quality assessment of next generation sequencing (NGS) files resulting from SELEX experiments.

The pipeline works with demultiplexed files in FASTQ format.
Sequences in the FASTQ files are expected to consist of a 5'-primer, random region, and 3'-primer.

The workflow localizes and trims off flanking regions, filters the sequences by quality and merges paired-end reads.
It also provides extensive information about the SELEX experiment's success and the quality of the NGS run.

FASTA files resulting from *selex-ngs-prep* are required by other SELEX analysis pipelines we developed, it can be used as a standalone analysis workflow to gain first insights in the experiment's data though.


## Prerequisites and Installation
### Conda Environment
A conda-like python environment manager is required.
The workflow was tested using conda 4.9.0, though any conda-like environment manager should work.

Info on the installation of conda can be found here: [conda.io](https://docs.conda.io/)

New packages can be installed to the base environment, which is activated by default.
We suggest to create a new environment to avoid conflicting packages.
#### Optional: Creating a new environment
Create a new environment in which the required packages will be installed.
Be sure to always activate this environment before using the workflow, as seen below.
```bash
# Creating a new environment
conda create -n selex-pipelines
# Activateing the new environment
conda activate selex-pipelines
```

### NextFlow Workflow Manager
Install the latest version of the bioinfomatic workflow manager NextFlow on your device in the environment of your choice.
NextFlow can be obtained from the channel Bioconda.
```bash
# Installing NextFlow
conda install -c bioconda nextflow
```

## Output
<!--- The workflow will put created new files into the folder *output/*.
Preprocessed FASTA files are in put into *output/preprocessed*.
Discarded sequences are put into *output/discarded*.
Plots and charts concerning the success of the SELEX experiment (selex round enrichment, nucleotide distribution along aptamers, nucleotide distribution over SELEX rounds) are put into *output/selex_analysis*.
Plots and charts concerning the quality of the next generation sequencing is put into *output/ngs_quality*. -->

The workflow generates a vast number of different files.


### NGS Quality

Quality plots of sequencing files are created using the plotQualityProfile function from the DADA2 Amplicon Sequencing Pipeline. <!-- TODO: Cite -->
Plots are created for raw files, as well as the fully preprocessed files.

<p float="left">
  <img src="https://github.com/hovercat/selex-ngs-prep/blob/main/example_output/analysis.ngs_quality/raw_plots/ngs_quality_forward.png?raw=true" width="49%" />
  <img src="https://github.com/hovercat/selex-ngs-prep/blob/main/example_output/analysis.ngs_quality/prepped_plots/ngs_quality_forward.png?raw=true" width="49%" />
</p>

On the left side is the quality profile plot before preprocessing. On the right side after preprocessing.

These plots are also created for every single SELEX round (before and after preprocessing) and are also written as HTML files.


### Preprocessing Performance

Plots and CSV-files are created to show the performance of the *selex-ngs-prep* preprocessing workflow.
It shows the how many reads were discarded in every step of the workflow for every SELEX round.

Here's an example plot:

<img src="https://github.com/hovercat/selex-ngs-prep/blob/main/example_output/analysis.preprocessing/preprocessing.perc.png?raw=true" alt="Preprocessing Performance Plot" width="90%"/>

### Preprocessed FASTA/FASTQ Files

The most important output files are the preprocessed (trimmed, filtered, merged) FASTA and FASTQ files, containing the aptamer random regions.
These files should be used for further analysis.

### Aptamer CSV Lists

Two files called *selex.aptamers.csv* and *selex.aptamers.rpm.csv* are created.
They contain one row for every unique random region encountered, along with read counts for every SELEX round.

Values in the rpm-file (rpm is reads per million) are based on the percentage a sequences takes of a SELEX round, multiplicated by 10^6. 

### SELEX Enrichment

Plots and CSV files are created, in which random regions are binned logarithmically by base 2 and base 10.

An example plot: 

<img src="https://github.com/hovercat/selex-ngs-prep/blob/main/example_output/analysis.selex_success/selex_success.log2.csv.png?raw=true" alt="SELEX Enrichment Plot on log2 base" width="60%" height="60%"/>

### Nucleotide Distribution

Barplots are created for the whole SELEX experiment, as well as every SELEX round on its own.
CSV files which are used for the barplots are outputted as well.

The plots are also summarized in an HTML-file.
As the images are stored directly in the HTML, showing it on GitHub is not possible due to file size limitations.
Please have a look into the example_output directory if you are interested.


<p float="left">
    <img src="https://github.com/hovercat/selex-ngs-prep/blob/main/example_output/analysis.nt_distribution/nt_distribution.png?raw=true" alt="Nucleotide Distribution Plot for SELEX Experiment" width="49%" height="49%"/>
    <img src="https://github.com/hovercat/selex-ngs-prep/blob/main/example_output/analysis.nt_distribution/R9.nt_distribution.png?raw=true" alt="Nucleotide Distribution Plot for Round R9" width="49%" height="49%"/>
</p>
The left plot shows the nt distribution over all rounds for a SELEX experiment.
The right plots shows the nt distribution along the random region of a specific SELEX round.


## Usage
### Directory Structure
The workflow is embedded in a directory structure.
By default new directories are created containing resulting files and a cache.

    .
    ├── bin/                # R and Python scripts for analysis
    ├── output_YOUREXPERIMENT/             
    │   ├── analysis.ngs_quality/       # HTML-Files and plots on sequencing quality
    │   ├── analysis.nt_distribution/   # HTML-Files and plots on nucleotide distribution on aptamers
    │   ├── analysis.preprocessing/     # Plots and charts on loss during preprocessing
    │   ├── analysis.selex_success/     # Plots and charts on aptamers enrichment
    │   ├── fastx*                      # Directories with preprocessed FASTQ/FASTA files
    │   ├── selex.aptamers.csv          # CSV chart containing all aptamers in the experiment
    │   └── selex.aptamers.rpm.csv      # CSV chart containing all aptamers in the experiment, counted in reads per million
    ├── work/               # cache directory containing intermediate files and temporary conda environments (may be deleted after)
    |   └── ...             
    ├── selex-ngs-prep.nf   # The data preparation workflow
    ├── nextflow.config     # default config file (should not be changed)
    ├── create_config.py    # Config creation wizard
    ├── LICENSE
    └── README.md
    
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
The default config-file with all options is shown here.
Please use the config creation wizard **create_config.py**, it will cover all of these options.

```groovy
// nextflow.config (default fallback config file)
params {
    experiment = "SELEX"                // Experiment's name

    input_dir = "./input_data"          // optional
    output_dir = "./output"             // optional
    fastq_pattern = null                // e.g. '*_{fwd,rev}.fastq'
    selex_rounds = null                 // list of selex rounds to take: ["R0", ...]

    trim_delimiter = null               // delimiter to be used to trim file names
    
    random_region = null                // length of random region
    random_region_min = null            // min length of random region
    random_region_max = null            // max length of random region
    
    primer5 = "ACGT"                    // 5'Primer aka Forward Primer
    primer3 = "TTTT"                    // 3'Primer aka Reverse Primer

    // further settings
    trim.max_error_rate = 0.2           // maximum error rate in the primers to be recognized
    filter.min_phred_quality = 30       // default average sequence read quality to include sequence
    merge.base_correction = true        // enables base correctio if paired-end reads mismatch and there are major quality differences
    merge.max_mismatches = 1            // allowed number of mismatches for merging
}
```

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
