# TASSELmanip
Python scripts for dealing with TASSEL outputs and other aspects of GWAS analyses

1. `tassel2geneTable_msu7.v1.0.py`

The primary program cited in Korinsak et al. plots Table S1 using an output file of Tassel (Bradbury et al Bioinformatics 2007) results and a list of MSU pseudogenes (http://rice.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/). You can alter the file names etc manually, or alternatively, the command `python tassel2geneTable_msu7.v1.0.py -h` in a terminal will give a list of argument options to set file names, thresholds (false discovery rate vs. bootstrap), output files prefix, and an option to output a full (somewhat unwieldy) results table. 

There should be five columns in the Tassel results file names: ['Trait', 'Marker', 'Chr', 'Pos', 'p']. If not contact me and I will help to generalise the program to accomodate different Tassel output formats. The Tassel file should be tab separated.

MSU gene list is as downloaded format (reduced here as 'genes.csv'). NB runtime may take several (i.e. 3-4) hours - may parallelise in future. The files 'genes.csv' & 'tassel.txt' are redacted/reduced test run files to be placed in the same directory as the program scripts. It should run out of the box.

The PDF files is the Table S1 formatted output; 'PREFIX_MLM_BLB_tassel_LDgenes_fdr_msu7.xlsx' is the same tabular output to XL; 'PREFIX_MLM_BLB_tassel_LDfuncGenes_fdr_msu7.xlsx' is a tabular version of the annotated regions; 'PREFIX_MLM_BLB_tassel_LDgenes_other_sign_locs_fdr_msu7.csv' is a list of SNPs not contained in MSU gene regions; 'PREFIXFDR_identifed_loci_msu7.csv' is a simpler format of identified MSU regions.  

2. `tassleManhattanPlotter.v1.0.py`

Makes Manhattan plots from Tassel output data. Draws both false discovery rate and bootstrap threshold ablines on the graphs. Python recognised colours [`colors`], figure DPI [`dpi`], figure dimensions [`figsize`], fontize [`fontsize`], and output directory [`outdir`] can be set in the script (under `#set params`).


3. `tasselQQ.v1.0.py`

Same as #2 - plus a line type option.


4. `LDheatmapMaker_cbar.v1.0.py` depends on R output from `mkHM.R`

Makes an annotated LD heatmap featuring genes, markers and/or regions of interest. In the script, under `#PARAMS` you can change:

Input file - `file`

Work folder - `workfldr`

File prefix - `prefix`

DPI - `dpi`

Chromosomal region - `region`: entered as a list of start-finish positions

Genes to highlight - `genes`: entered as a list

Interesting markers - `markers_to_highlight`: entered as list of positions

Colours of gene labels -`colours`: entered as a list from ['white', 'blue', 'red', 'black', 'green', 'dkred']

Regions of Interest (bold red triangles on heatmap) - `rois`: entered as list of lists of start-finish positions (can leave blank = [])

NB lists for `genes`, `markers_to_highlight`, `markers_to_highlight` must be of equal length


