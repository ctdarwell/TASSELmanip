# TASSELmanip
Python scripts for dealing with TASSEL outputs and other aspects of GWAS analyses

1. `tassel2geneTable_msu7.v1.0.py`

Plots Table S1 in Korinsak et al. using an output file of Tassel (Bradbury et al Bioinformatics 2007) results. Using a list of MSU pseudogenes (http://rice.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/). You can alter the file manually, or alternatively, the command `python tassel2geneTable_msu7.v1.0.py -h` in a terminal will give a list of argument options to set file names, thresholds (false discovery rate vs. bootstrap), output files prefix, and an option to output a full (somewhat unwieldy) results table. 

There should be five columns in the Tassel results file names: ['Trait', 'Marker', 'Chr', 'Pos', 'p']. If not contact me and I will help to generalise the program to accomodate different Tassel output formats. Tassel file should tab separated.

MSU gene list is as downloaded format.


2. `tassleManhattanPlotter.v1.0.py`

Makes Manhattan plots from Tassel output data. Draws both false discovery rate and bootstrap threshold lines on the graphs. Python recognised colours [`colors`], figure DPI [`dpi`], figure dimensions [`figsize`], fontize [`fontsize`], and output directory [`outdir`] can be set in the script.


