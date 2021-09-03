# *D*etection and *E*xtraction of *Ph*ages *T*ool (DEPhT)

DEPhT is a new tool for identifying prophages in bacteria, with a particular focus on Mycobacteria. It uses two cheap features to identify regions likely to contain prophages:
1.  Local average length of genes + intergenic regions, in 55-gene windows
2.  Local number of strand changes, in 55-gene windows

Each gene is assigned a probability of belonging to a prophage.  In general if a gene occurs in a 55-gene window where the average gene size is <800 bp and there are fewer than 10 strand changes, there is a very good chance it belongs to a prophage.  

This approach is not perfect, so another cheap approach is taken to improve the specificity.  MMseqs2 is used to cluster each gene against a database of clade-specific Mycobacterial core genes, to 50% identity, 80% coverage, e-value 0.001.  Regions of high likelihood prophage genes designated as non-core are taken as high probability prophages.

Depending on the selected runmode, these prophage regions are further scrutinized by functionally annotating them against HMMs of manually annotated mycobacteriophage phamilies from the [Actino_Draft Phamerator database](http://databases.hatfull.org/Actino_Draft).  Predicted prophages with too few high-probability hits into these HMMs may be culled as unlikely prophages.

Finally, the remaining prophages are subjected to a blastn-based attL/attR detection scheme that gives DEPhT superior edge detection than any tool we are aware of.

# Installation

While DEPhT can be installed and run by manually compiling each of its dependencies, by far the easiest approach is to use Anaconda:

    conda create -n depht python=3.9 -y && conda activate depht
    conda install -c bioconda -c conda-forge prodigal aragorn mmseqs2=13.45111 hhsuite=3 blast=2.9 -y
    git clone https://github.com/chg60/DEPhT && cd DEPhT
    pip install -r requirements.txt
    cd src
    
From there, print DEPhT's help menu by running it as a python module without arguments:

    python3 -m depht

Which shows something like this (default # CPUs will vary from computer to computer):

    usage: __main__.py [-h] -i INFILE [INFILE ...] [-f {fasta,genbank}] -o OUTDIR [-c CPUS] [-n] [-m {fast,normal,strict}]
                       [-s ATT_SENSITIVITY] [-d] [-v] [-t TMP_DIR] [-p PRODUCT_THRESHOLD] [-l LENGTH_THRESHOLD]

    DEPhT scans bacterial genomes looking for prophages. Regions identified as prophage candidates are further scrutinized, and
    attachment sites identified as accurately as possible before prophage extraction and generating the final report.

    optional arguments:
      -h, --help            show this help message and exit
      -i INFILE [INFILE ...], --infile INFILE [INFILE ...]
                            path to genome file(s) to scan for prophages
      -f {fasta,genbank}, --input-format {fasta,genbank}
                            input which format your input file(s) are in
      -o OUTDIR, --outdir OUTDIR
                            path where outputs can be written
      -c CPUS, --cpus CPUS  number of CPU cores to use [default: 4]
      -n, --no-draw         don't draw genome diagram for identified prophage(s)
      -m {fast,normal,strict}, --mode {fast,normal,strict}
                            select a runmode that favors speed or accuracy
      -s ATT_SENSITIVITY, --att_sensitivity ATT_SENSITIVITY
                            sensitivity parameter for att site detection.
      -d, --dump-data       dump all support data to outdir
      -v, --verbose         print progress messages as the program runs
      -t TMP_DIR, --tmp-dir TMP_DIR
                            temporary directory to use for file I/O [default: /tmp/prophicient]
      -p PRODUCT_THRESHOLD, --product-threshold PRODUCT_THRESHOLD
                            select a phage homolog product lower threshold
      -l LENGTH_THRESHOLD, --length-threshold LENGTH_THRESHOLD
                            select a minimum length for prophages [default: 20000]
