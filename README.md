# Prophicient

Prophicient (name not fully settled) will be a new tool for identifying prophages in bacteria.

Some of the extant tools are:
  * [PHASTER](https://phaster.ca)
  * [ProphET](https://github.com/jaumlrc/ProphET)
  * [VirSorter2](https://github.com/jiarong/VirSorter2)

These generally perform well at identifying core prophage regions, but struggle to accurately call
phage boundaries.  Additionally, the effective runtime (server queue + program runtime) for some of
these programs can exceed 24 hours.

The following workflow is our attempt at improving prophage boundary detection by anchoring the predictions to a reference genome:
  * Auto-annotate protein-coding genes with [Prodigal](https://github.com/hyattpd/Prodigal)
  * Auto-annotate tRNA genes with [Aragorn](http://www.ansikte.se/ARAGORN)
  * Identify phage-like genes using [HHsuite3](https://github.com/soedinglab/hh-suite) against a custom database
  * Identify regions with high phage character
  * Use k-mer counting / blastn vs. reference genome / EMBOSS tool to find most likely _attL_/_attR_
