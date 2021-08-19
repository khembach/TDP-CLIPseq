# TDP-CLIPseq

this repository contains the code for the analysis of TDP-43 CLIP-seq data. 
It includes a snakemake workflow (adjusted from [ARMOR](https://github.com/csoneson/ARMOR)) for read preprocessing, alignment and peak identification with CLIPper.
Reproducible peaks are identified with the [merge_peaks](https://github.com/YeoLab/merge_peaks) pipeline from the Yeo lab.

Directory: 

* Rmd/ contains the rmarkdown files for data analysis
* scripts/ contains code for the Snakemake workflow and other utility functions in R
* IDR_merge_peaks/ contains notes on how I installed/fixed bugs in merge_peaks and the configuration files for the samples
* envs/ contains the environment files for snakemake
* snoDB/ contains the [snoDB](http://scottgroup.med.usherbrooke.ca/snoDB/) annotations downloaded from http://scottgroup.med.usherbrooke.ca/snoDB/csv
