# Rātā Moehau population genomics and demographic modelling

This directory contains the scripts required for processing genomic resequencing data, aligning data to a reference genome, calling and filtering SNPs, 
and conducting PSMC demographic modelling analysis. 

The SNP filtering steps produces a large number of filtered variant options based on various filtering parameters (e.g., missing data, minimum coverage, 
linkage disequilibrium). We recommend that the outputs of these different filtering strategies be carefully considered, and parameters adjusted as appropriate
for your use case. A single variant file was selected for final downstream population genomic analyses. 

### Software

* [TrimGalore](https://github.com/FelixKrueger/TrimGalore) v0.6.4
* [BWA](https://github.com/lh3/bwa) v0.7.17
* [SamTools](https://github.com/samtools/) v1.10 and v1.15.1
* [BCFtools](https://github.com/samtools/bcftools) v1.15.1, v1.10.2 and v1.19
* [VCFtools](https://github.com/vcftools/vcftools) v0.1.15
* [Stacks](https://catchenlab.life.illinois.edu/stacks/) v.2.61
* [PLINK](https://www.cog-genomics.org/plink/) v1.09b6.16
* [PSMC](https://github.com/lh3/psmc) v0.6.5
* [R](https://cran.r-project.org/index.html) v4.4.1 



