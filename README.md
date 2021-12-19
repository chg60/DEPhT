# *D*etection and *E*xtraction of *Ph*ages *T*ool (DEPhT)

DEPhT is a new tool for identifying prophages in bacteria, with a particular focus on Mycobacteria. It uses two cheap features to identify regions likely to contain prophages:
1.  Local average length of genes + intergenic regions, in 55-gene windows
2.  Local number of strand changes, in 55-gene windows

Each gene is assigned a probability of belonging to a prophage.  In general if a gene occurs in a 55-gene window where the average gene size is <800 bp and there are fewer than 10 strand changes, there is a very good chance it belongs to a prophage.  

This approach is not perfect, so another cheap approach is taken to improve the specificity.  MMseqs2 is used to cluster each gene against a database of clade-specific Mycobacterial core genes, to 50% identity, 80% coverage, e-value 0.001.  Regions of high likelihood prophage genes designated as non-core are taken as high probability prophages.

Depending on the selected runmode, these prophage regions are further scrutinized by functionally annotating them against HMMs of manually annotated mycobacteriophage phamilies from the [Actino_Draft Phamerator database](http://databases.hatfull.org/Actino_Draft).  Predicted prophages with too few high-probability hits into these HMMs may be culled as unlikely prophages.

Finally, the remaining prophages are subjected to a blastn-based attL/attR detection scheme that gives DEPhT superior edge detection than any tool we are aware of.

# General Information

- Current version is 1.0.0
- Most recent stable version is 1.0.0
- Our current supported bacterial genera are the following.  If you wish to contact us regarding expanding this set, please reach out to either chg60@pitt.edu or laa89@pitt.edu.
    - Mycobacterium
    - Gordonia
    - Pseudomonas 

# Installation

DEPhT requires several dependancies, and while these can be installed and run by manually compiling each of its dependancies, by far the easiest approach is to use Anaconda:

    conda create -n depht python=3.9 -y && conda activate depht
    conda install -c bioconda -c conda-forge prodigal aragorn mmseqs2=13.45111 hhsuite=3 blast=2.9 -y

After installing dependancies, DEPhT can be installed with PyPI
    
    pip install depht

   
# Setup

As a part of its components, DEPhT requires several databases to be locally installed.

These databases can be manually installed at the command line.  The first step is to download these databases from DEPhT's Open Science Framework repository at https://osf.io/zt4n3/

From there, downloaded compressed databases should be moved to DEPhT's default data directory and de-compressed.  Example:

    mv ~/Downloads/Mycobacterium.zip ~/.depht/models/
    unzip ~/.depht/models/Mycobacterium.zip

# Running DEPhT

## Basics

After installation and setup, check that DEPhT can be run on the command line.  Typing 'depht' at the command line should display something similar to the following:

    usage: __main__.py [-h] [--model] [-c] [-n] [-m {fast,normal,strict}] [-s] [-d] [-v] [-t] [-p] [-l]
                       infile [infile ...] outdir


    DEPhT scans bacterial genomes looking for prophages. Regions identified as prophage candidates are further scrutinized, and
    attachment sites identified as accurately as possible before prophage extraction and generating the final report.

    optional arguments:
      -h, --help            show this help message and exit
      --model {Mycobacterium}
                            which local model should be used [default: Mycobacterium]
      -c , --cpus           number of CPU cores to use [default: 4]
      -n, --no-draw         don't draw genome diagram for identified prophage(s)
      -m {fast,normal,sensitive}, --mode {fast,normal,sensitive}
                            select a runmode that favors speed or accuracy
      -s , --att_sens       sensitivity parameter for att site detection.
      -d, --dump-data       dump all support data to outdir
      -v, --verbose         print progress messages as the program runs
      -t , --tmp-dir temporary directory to use for file I/O [default: ~/.depht.tmp]
      -p , --products       minimum number of phage homologs to report a prophage
      -l , --length         select a minimum length for prophages [default: 20000]

DEPhT requires two arguments:
1. A fasta/genbank-formatted sequence file
2. A location to store results
    
    depht path/to/sequence_file path/to/results/directory

In the event that a prophage region is discovered, or if `-d` is specified, DEPhT will create a directory at the specified location with the name of the inputted sequence.  DEPhT outputs an .html file with visualization of the prophage regions discovered as well as fasta and genbank-formatted sequence files for each of these regions.  See [below](#output) for more details about individual files in DEPhT's output.

Progress updates during DEPhT's runtime can be toggled with `-v`.

    depht path/to/sequence_file path/to/results/directory -v

The amount of resources (virtual cores) that DEPhT is allowed to utilize can be specified with `-c`.

    depht path/to/sequence_file path/to/results/directory -c 6

## Options

DEPhT was originally designed for the precise and efficient discovery and extraction of *Mycobacterium* prophages, but can be adapted for other genera with the `--model` flag.  See [above](#general-information) for genera currently supported, as well [here](#setup) for how to download databases for our other supported genera.

    depht path/to/sequence_file path/to/results/directory --model Pseudomonas

DEPhT is a multi-modal program, and can incoorporate more/less data into its discernment of prophage regions directly resulting in greater/lesser average runtime.  DEPhT's other runtime modes can be selected with `-m`, with options of:
-fast: DEPhT discovers prophage regions as fast as possible using gene size and transctiptional strand changes.  Regions are trimmed/dampened using the identified shell genome content of the selected genera.
-normal:  DEPhT discovers prophage regions as described above, and then seeks to differentiate between active and decrepit prophages by identifying phage homologs essential to phage viability.
-sensitive: DEPhT discovers prophage regions as described above, and then seeks to differentiate between active and decrepit prophages by identifying phage homologs with a clear function.

    depht path/to/sequence_file path/to/results/directory -m fast

DEPhT attempts to differentiate between active and decrepit prophages based on the number of identified prophage homologs in a region.  The number of products can be lowered or raised with the `-p` flag and is set at 5 by default.

    depht path/to/sequence_file> <path/to/results/directory -p 2

DEPhT employs a multi-feature scoring algorithm and a library of reference sequences to determine the best possible attachment sequence (or if there is no appropriate sequence).  This runtime of this is heavily influenced by the runtime of the BLASTn algorithm, and so the more/less sequence searched directly results in greater/lesser average runtime.  The amount of sequence that is searched for an attachment sequence can be controlled by the `-s` flag which defines the sequence space searched in units of 5000 base pairs and is set at 7 by default.

    depht path/to/sequence_file path/to/results/directory -s 2

DEPhT mandates a minimum length for prophage regions reported for output quality assurance.  This minimum length threshold can be lowered or raised with the `-l` flag, and is set at 20000 base pairs by default.

    depht path/to/sequence_file path/to/results/directory -l 10000

DEPhT utilizes various software that require outputs and data intermediates that are written to file.  These files are stored at a temporary directory that can be moved using the `-t` flag, and is set at ~/.depht/tmp/ by default.

    depht path/to/sequence_file path/to/results/directory -t /path/to/desired/directory

## Output

DEPhT's output is comprised of three main files
1. A .html file with a visualization of the discovered prophage regions
2. A .csv spreadsheet with the primary data used to discern prophage regions
3. A genbank-formatted sequence file of the inputted sequence file after analysis/annotation by DEPhT

DEPhT's graphical .html output displays a cirular input genome map and linear phage region genome map with [DnaFeaturesViewer](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer) as well as the coordinates of the regions discovered in a colored table with [pretty-html-table](https://github.com/sbi-rviot/ph_table).

![DEPhT's graphical output for prophages identified in M. *abscessus* strain GD43A](/resources/images/result_visualization_example.png)

In each of these genome maps and coordinate tables, prophage and/or protein-coding sequence features are colored green for forward-oriented features and colored red for reverse-oriented features.  Above those prophage features in the circular genome map is annotated the prophage region name as given by DEPhT.  Above those protein-coding features in the linear genome map(s) is annotated phage products as identified by DEPhT.

DEPhT's data .csv output contains data for each protein-coding feature in the inputted sequence file.

![DEPhT's data output for prophages identified in M. *abscessus* strain GD43A](/resources/images/data_spreadsheet_example.png)

The columns in this output are the following:
- Gene ID: A protein-coding feature ID assigned by DEPhT
- Start: The start coordinate of a feature in the input sequence
- End: The end coodinate of a feature in the input sequence
- Prediction: The probability of a feature belonging to a prophage as analyzed by DEPhT 
- Bacterial Homology: The identity of a feature as shell genome content as analyzed by DEPhT
- Phage Homology: The probability given by an alignment of a feature to a HMM profile of phage amino acid sequences
