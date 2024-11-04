# Master Notes for the RNASeq Analysis Pipeline
This file contains notest about the workflow for this RNASeq analysis of C. albicans. 

# Project Data
Data for this analysis was provided by Dr. Ronda Rolfes. The files we analyzed were RNA sequences present in wild type C. albicans grown on thiamine deficient media.

# FastQC
We first analyzed the raw sequence data with FastQC to determine an optimal trimming strategy.

```
module load fastqc

fastqc WTA2_1.fq
fastqc WTA2_2.fq

module unload fastqc
```

The FastQC analysis showed considerable noise in the first 15 bases of the reads, prompting trimming of these initial bases, as well as removal of Illumina adapter sequences.


# Trimming Sequences
We used our FastQC analysis to inform an optimal trimming strategy. The script used for running trimmomatic is contained in the RNA_trim file.


# Post-Trimming FastQC
We again used FastQC to analyze the quality of the reads, this time on the files generated from sequence trimming. The same code was used as above.
