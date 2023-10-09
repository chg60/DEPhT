# Detection and Extraction of Phages Tool (DEPhT)

DEPhT is a new tool for identifying prophages in bacteria, and was developed with a particular interest in being 
able to rapidly scan hundreds to thousands of genomes and accurately extract complete (likely active) prophages 
from them.

In brief DEPhT works by using genome architecture (rather than homology) to identify genomic regions likely to 
contain a prophage. Any regions with phage-like architecture (characterized as regions with high gene density 
and few transcription direction changes) are then further scrutinized using two passes of homology detection. 
The first pass identifies genes on putative prophages that are homologs of (species/clade/genus-level) conserved 
bacterial genes, and uses any such genes to disrupt the prophage prediction. The second pass (disabled in the 
'fast' runmode) identifies genes on putative prophages that are homologs of conserved, functionally annotated 
phage genes. Finally, prophage regions that got through the previous filters are subjected to a BLASTN-based 
attL/attR detection scheme that gives DEPhT better boundary detection than any tool we are aware of.

# Please cite

Gauthier CH, Abad L, Venbakkam AK, Malnak J, Russell DA, Hatfull GF. DEPhT: a novel approach for efficient prophage discovery and precise extraction. 
Nucleic Acids Research, Volume 50, Issue 13, 22 July 2022, Page e75, doi: [10.1093/nar/gkac273](https://doi.org/10.1093/nar/gkac273). PMID: 35451479.

# Installation

DEPhT runs natively on macOS and Linux operating systems, and in theory should work on Windows using 
[WSL](https://docs.microsoft.com/en-us/windows/wsl/install).

## Conda install

DEPhT has several dependencies, and as a result by far the easiest way to install it is to use 
[Anaconda](https://www.anaconda.com/products/individual) or the lightweight 
[Miniconda](https://docs.conda.io/en/latest/miniconda.html) with this single command:

    conda create -n depht -c laa89 -c bioconda -c conda-forge depht -y

It may take up to a couple of minutes to complete.

## Manual install

For users that would prefer to manage their own dependencies, you'll need to install each of the following:
- [BLAST](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/) 2.9 or higher
- [HHsuite3](https://github.com/soedinglab/hh-suite)
- [MMseqs2](https://github.com/soedinglab/mmseqs2) 13.45111
- [Prodigal](https://github.com/hyattpd/Prodigal)
- [Aragorn](http://www.ansikte.se/ARAGORN/Downloads/)
- [Python](https://www.python.org/downloads/) 3.7 or higher
- [DEPhT](https://pypi.org/project/depht/) 1.2.0 or higher
- [ClustalO](http://www.clustal.org/omega/#Download)            (for model-training)

All Python dependencies will be installed automatically when using pip to install the DEPhT package.
   
# Setup

DEPhT requires at least one genus-specific model to be installed before it will be able to run. At present, there 
are a few models available in [our repository at the Open Science Framework](https://osf.io/zt4n3). New models can
also be trained locally ([instructions below](#training-new-models)).

We now provide a helper script to facilitate installation of models we have trained and made available at OSF.io.

    usage: depht_fetch_model [-h] [-f] [-v] {Mycobacterium,Gordonia}

    Helper script to facilitate downloading models from OSF.io.
    
    positional arguments:
      {Mycobacterium,Gordonia}
                            choose a model to download
    
    optional arguments:
      -h, --help            show this help message and exit
      -f, --force           overwrite existing model files if they exist [default: False]
      -v, --verbose         print verbose output to stdout [default: False]

To download the Mycobacterium model, for example, run:

    depht_fetch_model Mycobacterium -v

Which should produce output similar to the following:

    Downloading model files from https://osf.io/download/aw4up/...
    Unzipping model files to ~/.depht/models/Mycobacterium...
    Removing zip file at ~/.depht/models/Mycobacterium.zip...
    Done.

Overwriting an existing model can be done by supplying the `-f`/`--force` argument.

Models trained using `depht_train` will be put in `~/.depht/models` by default. We are generally amenable to providing guidance
or assistance in the construction of new models - the easiest way to accomplish this is by [creating an issue](https://github.com/chg60/DEPhT/issues). Note that per Figure 8 of Gauthier et. al 2022, some genera are likely better suited than others for DEPhT model creation.


# Running DEPhT

## Basics

After installation and setup, check that DEPhT can be run on the command line. NOTE: If you installed using conda, 
you'll need to activate your environment first (e.g., `conda activate depht`). Typing `depht` at the commandline 
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

DEPhT will infer the input file type(s) when it parses the files, *not using the file extensions*. As far as we are 
aware, this makes DEPhT somewhat unusual among prophage-detection tools, as in a single run you can provide a set 
of files with multiple file formats. FASTA files will be treated as un-annotated and the sequences parsed from these 
input files will be auto-annotated prior to prophage detection. Genbank flatfiles will be treated as annotated genomes, 
and will therefore bypass the auto-annotation step and run ~20-30 seconds faster than their FASTA counterparts.

Run DEPhT on a single FASTA file like this (use your own file paths/extensions):

    depht /path/to/my/sequence.fasta /path/to/my/output/directory

Run DEPhT on a directory of FASTA files like this:

    depht /path/to/my/directory/*.fasta /path/to/my/output/directory

A large set of mixed FASTA (here using .fasta extension) and Genbank (here using .gbk extension) flatfiles can be run like this:

    depht /path/to/my/directory/*.fasta /path/to/my/directory/*.gbk /path/to/my/output/directory

In theory, you're limited only by the number of files your Terminal will let you expand by using `*`.

For Mac users who are uncomfortable with entering paths at the commandline, modern versions of macOS let you drag files 
from a Finder window into the Terminal and will automatically populate the path in the Terminal for you. Some Linux 
distributions may also support this kind of action.

In the event that a prophage region is discovered, or if the `-d` argument is specified, DEPhT will create a 
directory at the specified output directory for each of the input sequences. For those sequences that have predicted 
prophages, DEPhT will write an .html file with a visualization of the discovered prophage region(s). It will also 
output a FASTA (sequence) file and a Genbank (annotation) file for each extracted prophage sequence. See 
[below](#output) for more details DEPhT's output files.

Progress updates during DEPhT's runtime can be toggled with `-v`.

    depht /path/to/my/sequence.fasta /path/to/my/output/directory -v

DEPhT will use all locally available CPU cores by default. You can limit the number of allowed CPU cores with the 
`-c`/`--cpus` argument:

    depht /path/to/my/sequence.fasta /path/to/my/output/directory -c 4


## Other Options

What follows is a description of DEPhT's optional arguments. These are described in isolation, but can be mixed and 
matched using different values to specifically tune the behavior of DEPhT to suit your needs. Default parameters 
were all set to optimize performance in _Mycobacterium_ genomes.

### Model Selection

DEPhT was originally designed for the precise and efficient discovery and extraction of *Mycobacterium* prophages, 
but can be adapted for other genera with the `--model` flag.  See [above](#setup) for instructions to download 
models that we have already trained, and [below](#general-information) for the list of currently available models.

If you have more than one model installed locally, you will need to tell DEPhT which model you'd like to use. 
Otherwise, it will choose one more-or-less at random, which may result in unexpectedly low-quality outputs.

    depht /path/to/my/sequence.fasta /path/to/my/output/directory --model Pseudomonas

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

    depht /path/to/my/sequence.fasta /path/to/my/output/directory -m fast

Alternatively, if one wants only the most likely prophages, with as many detailed functional annotations as possible,
they might run:

    depht /path/to/my/sequence.fasta /path/to/my/output/directory -m sensitive

### Product Threshold

In normal and sensitive runmodes, DEPhT attempts to differentiate between active and defective prophages based on the 
number of identified prophage homologs in a region. This number of phage products can be raised or lowered by using 
the `-p` argument. In normal mode, the default value is 5; in sensitive mode, it is 10. If one feels that the 
default value is too high and would rather use 2 for example, this can be done by running:

    depht /path/to/my/sequence.fasta /path/to/my/output/directory -p 2

### Attachment Site Tuning

DEPhT employs a multi-feature scoring algorithm and a library of reference sequences to determine the best possible 
attachment site core (or if there is no appropriate sequence). The runtime of this component is heavily influenced by 
the runtime of the BLASTN algorithm, so runtime scales with the amount of sequence that is searched. However, the 
precision of extraction may benefit from searching a larger sequence space, particularly in genera where few 
high-quality reference sequences are available. The sequence space that is searched for an attachment site core 
can be controlled by using the `-s` flag, which acts as a multiplier against 5000 bp. By default, DEPhT uses `-s 7`, 
which corresponds to a search space of up to 7 x 5,000 = 35,000 bp at the left and right ends of each identified 
prophage. This can be raised to 50,000 bp by setting `-s` to 10, at the expense of some additional runtime:

    depht /path/to/my/sequence.fasta /path/to/my/output/directory -s 10

### Prophage Size

DEPhT mandates a minimum length for prophage regions reported for output quality assurance. This minimum length 
threshold can be lowered or raised with the `-l` flag, and is set at 20,000 base pairs by default - just over half 
the length of the shortest known Mycobacterium prophage. Reduce this threshold to 10,000 bases like this:

    depht /path/to/my/sequence.fasta /path/to/my/output/directory -l 10000

### Temporary Directory

DEPhT (and its dependencies) use many temporary data files. These files are stored in a temporary directory, 
and removed once DEPhT finishes running. By default, DEPhT will use `~/.depht/tmp`, but can use any other 
directory that your user account has read/write permissions in, by using the `-t` argument:

    depht /path/to/my/sequence.fasta /path/to/my/output/directory -t /path/to/temporary/directory


# Output

DEPhT's output consists of three main files:
1. An `.html` file with a visualization of the discovered prophage regions
2. A `.csv` spreadsheet with the primary data used to discern prophage regions - one file per contig
3. A `.gbk` Genbank flatfile with DEPhT's annotation of the inputted sequence - one file per contig

DEPhT's graphical `.html` output displays a circular input genome map and linear phage region genome map with 
[DnaFeaturesViewer](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer) as well as the coordinates of the 
regions discovered in a colored table with [pretty-html-table](https://github.com/sbi-rviot/ph_table).

![DEPhT's graphical output for prophages identified in *M. abscessus* strain GD43A]
(/resources/images/result_visualization_example.png)

In each of these genome maps and coordinate tables, prophage and/or protein-coding sequence features are colored 
green for forward-oriented features and colored red for reverse-oriented features. Above those prophage features in 
the circular genome map is annotated the prophage region name as given by DEPhT. Above those protein-coding features 
in the linear genome map(s) is annotated phage products as identified by DEPhT.

DEPhT's data `.csv` output contains data for each protein-coding feature in the inputted sequence file.

![DEPhT's data output for prophages identified in *M. abscessus* strain GD43A]
(/resources/images/data_spreadsheet_example.png)

The columns in this output are the following:
- Gene ID: A protein-coding feature ID assigned by DEPhT
- Start: The start coordinate of a feature in the input sequence
- End: The end coordinate of a feature in the input sequence
- Prediction: Probability that this gene is part of a prophage, from the Gene Size/TDC Classifier
- Bacterial Homology: 1 if DEPhT identifies this gene as part of the bacterial accessory genome, else 0
- Phage Homology: HHSearch probability for genes with high-confidence hits to phage HMMs, else 0


# Training New Models

Models can be trained using the `depht_train` package, installed as part of DEPhT.
    
What follows will describe the workflow for training new models, as well as explain the thought process.
    
## Selection of Training Genomes

This is by far the highest hurdle for training new models. The better the training genomes are selected, the better
the model will perform. We *highly* recommend only training against completely sequenced bacteria and manually
annotated phages.

There's an important tradeoff you'll need to make when training models: volume of data versus quality of data. A
relatively small dataset (~100 phages and 30-45 bacteria) can yield incredibly high-quality models if the genomes are
chosen well and especially if the phage genomes are well-annotated. Assuming all the training data is high-quality,
increasing the amount of training data will likely improve the quality of predictions made by DEPhT, with the caveat
that larger models will necessarily increase the DEPhT runtime, which will be most noticeable in the fast runmode.

Ok so let's suppose you want to train a new model for Mycobacteria. A good start would be to head to 
[PATRIC](https://www.patricbrc.org/view/Taxonomy/2#view_tab=taxontree) and navigate to the Mycobacteriaceae.

### Retrieve Bacterial Genomes

In the taxonomy tree, the steps to get here are: 

Terrabacteria group >> Actinobacteria >> Actinomycetia >> Corynebacteriales >> Mycobacteriaceae

The red box below shows where to click to get to the home page for the family or genus of interest.

![patric mycobacteriaceae](/resources/images/patric_mycobacteriaceae.png)

From there, navigate to the "Genomes" tab to see all the available genomes in the chosen taxon. Click "Filters", and
a good choice might be to select only those genomes where "Genome Status" is "Complete", and "Reference Genome" is
either "Representative" or "Reference", and "Genome Quality" is "Good". Hit "Apply" to apply those filters. You can 
download FASTA files for these genomes by selecting all the genomes in the table, and clicking the "DWNLD" button.

![patric download](/resources/images/patric_download.png) 

Click "More Options", and in the popup dialog box, check the box next to "Genomic Sequences in FASTA (\*.fna)" 
before pressing "Download". 

![patric dialog box](/resources/images/patric_download_options.png)

Of course, you are free to add any additional genomes you'd like to better 
populate the spectrum of diversity in the genus. In our case, we added several _Mycobacterium abscessus_ strains to 
fill in the so-called _Mycobacterium abscessus_ complex (MAC).

### Check Bacteria for Prophages

Ideally, you'll run these genomes through PHASTER or some other prophage prediction tool to get the approximate
coordinates of any complete prophages in these strains, and recording them in a CSV file that you'll pass to the 
training module. The coordinates don't have to be perfect, though the better they are the better the resultant model 
will perform. This step will give the model an idea what integrated prophages are supposed to look architecturally, as 
opposed to only knowing what extracted phages/prophages look like.

![example csv](/resources/images/prophage_csv.png)

### Retrieve Phage Genomes

Lastly, you'll need to retrieve functionally annotated phages from Genbank or elsewhere. Like the bacteria, it's
important that these phages represent the spectrum of diversity of phages infecting hosts in the genus. Ideally there
will also be groups of at least somewhat-related phages in this dataset. Phage annotations need to be in GenBank flatfile
format.

## Running the Training Pipeline

DEPhT models are comprised of four main components, which can be built from curated phage and bacterial sequences.
In the `depht_train` package, several pipelines exist to convert and process sequence data as well as store these data
in formatted files and databases recognized by DEPhT.  The overall workflow is illustrated below.

![depht train workflow](/resources/images/github_training_schema.png)

This training workflow is available in the `depht_train` package as a single pipeline. The only required arguments are:
1. a name for the new model
2. path to a directory containing functionally annotated phage genomes for the genus of interest (GenBank flatfiles)
3. path to a directory containing bacterial genomes for the genus of interest (FASTA nucleotide or GenBank flatfiles)

Run the pipeline like this:

    depht_train create_model model_name /path/to/annotated/phage/genomes /path/to/bacterial/genomes
    
If you're trying to create a new model with the same name as an existing one, `depht_train` will not overwrite the
existing model by default, but it will then force you to pick a new name. If you'd like to overwrite the existing model,
you can do so with the `-f`/`--force` argument:

    depht_train create_model model_name /path/to/annotated/phage/genomes /path/to/bacterial/genomes -f
    
If one or more of your bacterial genomes has one or more known (or probable) prophage(s) in it, you can provide a CSV file
formatted [as above](#check-bacteria-for-prophages), using the `--prophage-coords` argument:

    depht_train create_model model_name /path/to/annotated/phage/genomes /path/to/bacterial/genomes --prophage-coords /path/to/prophage_coords.csv

When much is known about the taxonomy of a set of bacteria provided for the creation of a DEPhT model, formal clades or
taxa assigned to the bacteria may be much more informative and biologically relevant than an auto-generated cluster
schema. To assert a certain clade/cluster schema on the inputted bacterial sequences for the purposes of defining shell genome
content, you can create a CSV table mapping the given bacteria to a clade name or identifier. The table should be formatted
with a 'Name' and 'Cluster' header like the following:

![bacterial clusters table example](/resources/images/bacterial_clusters_table.png)

This CSV file can be provided to the `create_model` pipeline using the `--bacteria-clusters` argument:

    depht_train create_model model_name /path/to/annotated/phage/genomes /path/to/bacterial/genomes 
    --bacteria-clusters /path/to/bacteria_clusters.csv
 
Training a model consists of several computationally expensive steps, and as such the amount of time it takes to train
a model is highly variable, but generally influenced in these ways:
1. more genomes --> longer training time (and likely `depht` runtime)
2. more CPU cores --> shorter training time (and likely `depht` runtime)

Most new models will likely take somewhere between 15 minutes and an hour to train.

# Examples

We have some example genomes and expected outputs to make it easy for end-users to verify that DEPhT (and the Mycobacterium model)
are installed correctly. These genomes and their expected outputs can be found in the GitHub repository in resources > examples.

# General Information

- Current version is 1.2.0
- Most recent stable version is 1.1.8
- We currently have models available for these bacterial genera:
    - [Mycobacterium](https://osf.io/aw4up/download)
    - [Gordonia](https://osf.io/djwsb/download)
    - [Pseudomonas](https://osf.io/5puze/download)
- If you'd like to contact us about expanding this set, please email either chg60@pitt.edu or laa89@pitt.edu
