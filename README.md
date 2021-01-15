# ProphageHunter

ProphageHunter ~is~ will be a new tool for identifying prophages in bacteria.

Extant tools (e.g. [ProphET](https://github.com/jaumlrc/ProphET), [PHASTER](https://phaster.ca)) struggle to accurately identify the boundaries of predicted prophages.

The following workflow is our attempt at improving prophage boundary detection by anchoring the predictions to a reference genome:
  * Auto-annotate the new genome, using [prokka](https://github.com/tseemann/prokka)
  * Identify "phage-y" genes using HHsuite3 (https://github.com/soedinglab/hh-suite)
  * Identify regions with higher-than-background proportion of "phage-y" genes
  * Use blastn to align "phage-y" regions to the reference genome for enhanced boundary detection
  * (optional) run Aragorn/tRNAscan-SE?
  * (optional) EMBOSS tool to find terminal repeats?
