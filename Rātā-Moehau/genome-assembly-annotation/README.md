# Genome assembly and annotation

This directory contains all scripts required for assembly and annotation of the rātā Moehau genome, from pre-processing raw short- and long-read data, initial assembly, scaffolding, annotation, QC, and post-processing steps. 

We recommend that various tools are used for assessing the quality of the assembly at various steps in the process, including gathering standard assembly metrics (assembly length, contig/scaffold length, contig/scaffold N50, distribution of contig/scaffold lengths, number/length of gaps), ortholog presence via BUSCO or Compleasm, alignment of short-/long-read data to assess coverage spikes or gaps, and assessing assembly completeness with Merqury. We also used TIDK to look for potential telomeres, and FCS-GX to check for contamination. 

## Software requirements

### Data processing

* [Guppy](https://nanoporetech.com/software/other/guppy) v6.2.1
* [nanoQC](https://github.com/wdecoster/nanoQC) v0.9.4
* [NanoPlot](https://github.com/wdecoster/NanoPlot)
* [TrimGalore](https://github.com/FelixKrueger/TrimGalore) v0.6.4
* [Jellyfish](https://github.com/gmarcais/Jellyfish) v2.3.0
* [Porechop](https://github.com/rrwick/Porechop) v0.2.4
* [nanofilt](https://github.com/wdecoster/nanofilt) v2.6.0

### Genome assembly 

* [shasta](https://github.com/paoloshasta/shasta) v0.10.0

### Assembly duplicate purging

Following the [purge-dups pipeline](https://github.com/dfguan/purge_dups)

* [minimap2](https://github.com/lh3/minimap2) v2.24

### Assembly polishing

* [Racon](https://github.com/lbcb-sci/racon) v1.5.0
* minimap2 v2.24
* [medaka](https://github.com/nanoporetech/medaka) v1.11.1
* [SAMtools](https://github.com/samtools/samtools) v1.13
* [BWA](https://github.com/lh3/bwa) v0.7.17

### Assembly scaffolding

Following the [Dovetail Omni-C scaffolding pipeline](https://omni-c.readthedocs.io/en/latest/index.html)

* BWA v0.7.17
* [picard](https://broadinstitute.github.io/picard/) v2.21.8
* SAMtools v1.13
* [samblaster](https://github.com/GregoryFaust/samblaster) v0.1.26
* [fastp](https://github.com/OpenGene/fastp) v0.23.2
* SAMtools v1.15.1
* [pairtools](https://pairtools.readthedocs.io/en/latest/) v1.0.2
* [YaHS](https://github.com/c-zhou/yahs)
* [juicer_tools](https://github.com/aidenlab/JuicerTools) v1.9.9
* [LASTZ](https://github.com/lastz/lastz) v1.04.03

### Repeat masking

* [RepeatModeler](https://github.com/Dfam-consortium/RepeatModeler) v2.0.3
* [RepeatMasker](https://www.repeatmasker.org/) v4.1.0
* [SeqKit](https://bioinf.shenwei.me/seqkit/) v2.4.0
* [bioawk](https://github.com/lh3/bioawk) v1.0
* [RMBlast](https://www.repeatmasker.org/rmblast/) v2.10.0
* [seqtk](https://github.com/lh3/seqtk) v1.4
* [BEDTools](https://bedtools.readthedocs.io/en/latest/) v2.30.0

### Gene annotation

* [Kraken2](https://github.com/DerrickWood/kraken2) v2.1.3
* [Bracken](https://github.com/jenniferlu717/Bracken) v2.7
* [SortMeRNA](https://github.com/sortmerna/sortmerna) v4.3.6
* [STAR](https://github.com/alexdobin/STAR) v2.7.10b
* SAMtools/1.19
* picard/2.26.10
* bioawk/1.0
* [BRAKER](https://github.com/Gaius-Augustus/BRAKER) v3.0.8
* [AUGUSTUS](https://github.com/Gaius-Augustus/Augustus) v3.5.0
* [compleasm](https://github.com/huangnengCSU/compleasm) v0.2.5
* [TSEBRA](https://github.com/Gaius-Augustus/TSEBRA) v1.1.2.5
* [AGAT](https://agat.readthedocs.io/en/latest/index.html) v1.0.0
* [gffread](https://github.com/gpertea/gffread) v0.12.7
* [genometools](https://github.com/genometools/genometools) v1.6.1
* [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) v2.16.0
* [FASTX-Toolkit](https://github.com/agordon/fastx_toolkit) v0.0.14
* [DIAMOND](https://github.com/bbuchfink/diamond) v2.1.10
* [InterProScan](https://www.ebi.ac.uk/interpro/search/sequence/) v5.66-98.0

### R packages 

R version 4.4.0

* dplyr
* kableExtra
* ggplot2
* tidyr
* readr
