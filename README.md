# *D*etection and *E*xtraction of *Ph*ages *T*ool (DEPhT)

DEPhT is a new tool for identifying prophages in bacteria, and was developed with a particular interest in being 
able to rapidly scan hundreds to thousands of genomes and accurately extract complete (likely active) prophages 
from them.

A detailed manuscript has been submitted to Nucleic Acids Research, but in brief DEPhT works by using genome 
architecture (rather than homology) to identify genomic regions likely to contain a prophage. Any regions with 
phage-like architecture (characterized as regions with high gene density and few transcription direction changes) 
are then further scrutinized using two passes of homology detection. The first pass identifies genes on putative 
prophages that are homologs of (species/clade/genus-level) conserved bacterial genes, and uses any such genes to 
disrupt the prophage prediction. The second pass (disabled in the 'fast' runmode) identifies genes on putative 
prophages that are homologs of conserved, functionally annotated phage genes. Finally, prophage regions that got 
through the previous filters are subjected to a BLASTN-based attL/attR detection scheme that gives DEPhT better 
boundary detection than any tool we are aware of.


# Installation

DEPhT has several dependencies, and as a result by far the easiest way to install it is to use 
[Anaconda](https://www.anaconda.com/products/individual) or the lightweight 
[Miniconda](https://docs.conda.io/en/latest/miniconda.html) with this single command:

    conda create -n depht -c laa89 -c bioconda -c conda-forge depht -y

It may take up to a couple of minutes to complete.

   
# Setup

DEPhT requires at least one genus-specific model to be installed before it will be able to run. At present, there 
are a few models available in [our repository at the Open Science Framework](https://osf.io/zt4n3). Creating new 
models is currently non-trivial, but we have a script nearly finished that should make the process much simpler. 
This script (and relevant documentation) should be available by early January 2022.

Once a model has been downloaded (easiest way is through a web browser), it needs to be decompressed and moved into 
a directory for DEPhT. For example, if you downloaded the Mycobacterium model:

    if ! [[ -d ~/.depht/models ]]; then
        mkdir -p ~/.depht/models
    fi

    mv ~/Downloads/Mycobacterium.zip ~/.depht/models/
    unzip ~/.depht/models/Mycobacterium.zip -d ~/.depht/models/Mycobacterium

Models trained using `depht_utils` will be put in this directory by default. We are generally amenable to aiding in 
the construction of new models - the easiest way to accomplish this is by emailing either chg60@pitt.edu or 
laa89@pitt.edu. Note that some genera are better suited than others for DEPhT model creation.


# Running DEPhT

## Basics

After installation and setup, check that DEPhT can be run on the command line. Typing `depht` at the commandline 
should display something similar to the following (number of CPUs and models available will vary):

    usage: depht [-h] [--model] [-c] [-n] [-m {fast,normal,strict}] [-s] [-d] [-v] [-t] [-p] [-l]
                 infile [infile ...] outdir

    DEPhT scans bacterial genomes looking for prophages. Regions identified as prophage 
    candidates are further scrutinized, and attachment sites identified as accurately as 
    possible before prophage extraction and generating the final report.

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
      -t , --tmp-dir        temporary directory to use for file I/O [default: ~/.depht/tmp]
      -p , --products       minimum number of phage homologs to report a prophage
      -l , --length         select a minimum length for prophages [default: 20000]

In order to run DEPhT, you will need to provide two arguments:
1. One or more genome sequences in either FASTA or Genbank flatfile format
2. A desired output directory 

DEPhT will infer the input file type(s) when it parses the files. This makes DEPhT somewhat unusual among 
prophage-detection tools, as in a single run you can provide a set of files with multiple file formats. FASTA files 
will be treated as un-annotated and the sequences parsed from these input files will be auto-annotated prior to 
prophage detection. Genbank flatfiles will be treated as annotated genomes, and will therefore bypass the 
auto-annotation step and run ~20-30 seconds faster than their FASTA counterparts.

Run DEPhT on a single FASTA file like this:

    depht /path/to/my/sequence.fna /path/to/my/output/directory

Run DEPhT on a directory of FASTA files like this:

    depht /path/to/my/directory/*.fna /path/to/my/output/directory

A large set of mixed FASTA and Genbank flatfiles can be run like this:

    depht /path/to/my/directory/*.fna /path/to/my/directory/*.gb /path/to/my/output/directory

In theory, you're limited only by the number of files your Terminal will let you expand by using `*`.

For Mac users who are uncomfortable with entering paths at the commandline, modern versions of MacOS let you drag files 
from a Finder window into the Terminal and will automatically populate the path in the Terminal for you. Some Linux 
distributions may also support this kind of action.

In the event that a prophage region is discovered, or if the `-d` argument is specified, DEPhT will create a 
directory at the specified output directory for each of the input sequences. For those sequences that have predicted 
prophages, DEPhT will write an .html file with a visualization of the discovered prophage region(s). It will also 
output a FASTA (sequence) file and a Genbank (annotation) file for each extracted prophage sequence. See 
[below](#output) for more details DEPhT's output files.

Progress updates during DEPhT's runtime can be toggled with `-v`.

    depht /path/to/my/sequence.fna /path/to/my/output/directory -v

The amount of resources (CPU cores) that DEPhT is allowed to utilize can be specified with `-c`. Note that some of 
DEPhT's dependencies utilize hyper-threading, so on most modern computers DEPhT will utilize 2 threads per specified 
CPU core.

    depht /path/to/my/sequence.fna /path/to/my/output/directory -c 6


## Other Options

What follows is a description of DEPhT's optional arguments. These are described in isolation, but can be mixed and 
matched using different values to specifically tune the behavior of DEPhT to suit your needs. Default parameters 
were all set to optimize performance in _Mycobacterium_ genomes.

### Model Selection

DEPhT was originally designed for the precise and efficient discovery and extraction of *Mycobacterium* prophages, 
but can be adapted for other genera with the `--model` flag.  See [above](#setup) for instructions to download 
models that we have already trained, and [below](#general-information) for the list of currently available models.

If you have more than one model installed locally, you will need to tell DEPhT which model you'd like to use. 
Otherwise, it will choose one more-or-less at random, which will likely result in unexpectedly low-quality outputs.

    depht /path/to/my/sequence.fna /path/to/my/output/directory --model Pseudomonas

### Runmode Selection

DEPhT has multiple runmodes, intended to serve as a dial tuning the trade-off between runtime and accuracy. The `-m` 
argument lets you select one of the available runmodes:

- fast: DEPhT discovers prophage regions as fast as possible using gene size and transcription direction changes. 
  Regions are trimmed using the identified shell genome content of the selected genera, and an effort is made to 
  identify attL/attR, but are likely not as accurate as in the other runmodes.
- normal: DEPhT discovers prophage regions as in fast mode, then tries to differentiate between active and defective 
  prophages by identifying homologs of phage genes essential for viability.
- sensitive: DEPhT discovers prophage regions as in normal mode, and then tries to further differentiate between active 
  and defective prophages by identifying homologs of phage genes with a consensus annotated function.

DEPhT will run in normal mode by default (e.g. if `-m` is not given), but if one is interested in getting an 
estimate of the number of prophages as quickly as possible, they may run DEPhT like this:

    depht /path/to/my/sequence.fna /path/to/my/output/directory -m fast

Alternatively, if one wants only the most likely prophages, with as many detailed functional annotations as possible,
they might run:

    depht /path/to/my/sequence.fna /path/to/my/output/directory -m sensitive

### Product Threshold

In normal and sensitive runmodes, DEPhT attempts to differentiate between active and defective prophages based on the 
number of identified prophage homologs in a region. This number of phage products can be raised or lowered by using 
the `-p` argument. In normal mode, the default value is 5; in sensitive mode, it is 10. If one feels that the 
default value is too high and would rather use 2 for example, this can be done by running:

    depht /path/to/my/sequence.fna /path/to/my/output/directory -p 2

### Attachment Site Tuning

DEPhT employs a multi-feature scoring algorithm and a library of reference sequences to determine the best possible 
attachment site core (or if there is no appropriate sequence). The runtime of this component is heavily influenced by 
the runtime of the BLASTN algorithm, so runtime scales with the amount of sequence that is searched. However, the 
precision of extraction may benefit from searching a larger sequence space, particularly in genera where few 
high-quality reference sequences are available. The sequence space that is searched for an attachment site core 
can be controlled by using the `-s` flag, which acts as a multiplier against 5000 bp. By default, DEPhT uses `-s 7`, 
which corresponds to a search space of up to 7 x 5,000 = 35,000 bp at the left and right ends of each identified 
prophage. This can be raised to 50,000 bp by setting `-s` to 10, at the expense of some additional runtime:

    depht /path/to/my/sequence.fna /path/to/my/output/directory -s 10

### Prophage Size

DEPhT mandates a minimum length for prophage regions reported for output quality assurance. This minimum length 
threshold can be lowered or raised with the `-l` flag, and is set at 20,000 base pairs by default - just over half 
the length of the shortest known Mycobacterium prophage. Reduce this threshold to 10,000 bases like this:

    depht /path/to/my/sequence.fna /path/to/my/output/directory -l 10000

### Temporary Directory

DEPhT utilizes various software that require outputs and data intermediates that are written to file. These files 
are stored in a temporary directory, and removed once DEPhT finishes running. By default, DEPhT will use `~/.depht/tmp` 
can use any other directory that your user account has read/write permissions in, by using the `-t` argument.

    depht /path/to/my/sequence.fna /path/to/my/output/directory -t /path/to/temporary/directory


# Output

DEPhT's output consists of three main files:
1. An `.html` file with a visualization of the discovered prophage regions
2. A `.csv` spreadsheet with the primary data used to discern prophage regions - one file per contig
3. A `.gbk` Genbank flatfile with DEPhT's annotation of the inputted sequence - one file per contig

DEPhT's graphical `.html` output displays a cirular input genome map and linear phage region genome map with 
[DnaFeaturesViewer](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer) as well as the coordinates of the 
regions discovered in a colored table with [pretty-html-table](https://github.com/sbi-rviot/ph_table).

![DEPhT's graphical output for prophages identified in M. *abscessus* strain GD43A](/resources/images/result_visualization_example.png)

In each of these genome maps and coordinate tables, prophage and/or protein-coding sequence features are colored 
green for forward-oriented features and colored red for reverse-oriented features. Above those prophage features in 
the circular genome map is annotated the prophage region name as given by DEPhT. Above those protein-coding features 
in the linear genome map(s) is annotated phage products as identified by DEPhT.

DEPhT's data `.csv` output contains data for each protein-coding feature in the inputted sequence file.

![DEPhT's data output for prophages identified in M. *abscessus* strain GD43A](/resources/images/data_spreadsheet_example.png)

The columns in this output are the following:
- Gene ID: A protein-coding feature ID assigned by DEPhT
- Start: The start coordinate of a feature in the input sequence
- End: The end coodinate of a feature in the input sequence
- Prediction: The probability of a feature belonging to a prophage as analyzed by DEPhT 
- Bacterial Homology: The identity of a feature as shell genome content as analyzed by DEPhT
- Phage Homology: The probability given by an alignment of a feature to a HMM profile of phage amino acid sequences


# General Information

- Current version is 1.0.1
- Most recent stable version is 1.0.1
- We currently have models available for these bacterial genera:
    - Mycobacterium
    - Gordonia
    - Pseudomonas
- If you'd like to contact us about expanding this set, please email either chg60@pitt.edu or laa89@pitt.edu
