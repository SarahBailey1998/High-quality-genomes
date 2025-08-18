# Rātā Moehau population genomics and demographic modelling

This directory contains the scripts required for processing genomic resequencing data, aligning data to a reference genome, calling and filtering SNPs, 
and conducting PSMC demographic modelling analysis. 

The SNP filtering steps produces a large number of filtered variant options based on various filtering parameters (e.g., missing data, minimum coverage, 
linkage disequilibrium). We recommend that the outputs of these different filtering strategies be carefully considered, and parameters adjusted as appropriate
for your use case. A single variant file was selected for final downstream population genomic analyses. 

## Software

* [TrimGalore](https://github.com/FelixKrueger/TrimGalore) v0.6.4

