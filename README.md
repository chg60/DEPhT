# *D*etection and *E*xtraction of Pro*ph*ages *T*ool (DEPhT)

DEPHT is a new tool for identifying prophages in bacteria, with a particular focus on Mycobacteria.

Popular extant tools primarily rely on homology detection techniques for identifying likely prophage regions:
  * [PHASTER](https://phaster.ca)
  * [ProphET](https://github.com/jaumlrc/ProphET)
  * [VirSorter2](https://github.com/jiarong/VirSorter2)

These techniques typically have good sensitivity, because every gene is tested for homology to known (pro)phage genes.  However, this sensitivity comes with a high computational cost (i.e. slow, without a high end server).  Additionally, these tools struggle with phage boundary detection - only rarely are the attL/attR properly called.

DEPHT uses two cheap features to identify regions likely to contain prophages:
1.  Local average length of genes + intergenic regions, in 55-gene windows
2.  Local number of strand changes, in 55-gene windows

Each gene is assigned a probability of belonging to a prophage.  In general if a gene occurs in a 55-gene window where the average gene size is <800 bp and there are fewer than 10 strand changes, there is a very good chance it belongs to a prophage.  

This approach is not perfect, so another cheap approach is taken to improve the specificity.  MMseqs2 is used to cluster each gene against a database of clade-specific Mycobacterial core genes, to 50% identity, 80% coverage, e-value 0.001.  Regions of high likelihood prophage genes designated as non-core are taken as high probability prophages.

Depending on the selected runmode, these prophage regions are further scrutinized by functionally annotating them against HMMs of manually annotated mycobacteriophage phamilies from the [Actino_Draft Phamerator database](http://databases.hatfull.org/Actino_Draft).  Predicted prophages with too few high-probability hits into these HMMs may be culled as unlikely prophages.

Finally, the remaining prophages are subjected to a blastn-based attL/attR detection scheme that gives DEPHT superior edge detection than any tool we are aware of.
