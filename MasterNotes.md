# Master Notes for the RNASeq Analysis Pipeline
This file contains notes about the workflow for this RNASeq analysis of *C. albicans*. 

# Project Objectives
This project seeks to determine the differential gene expression of the human commensal yeast, *Candida albicans*, when grown in the presence and absence of thiamine. Six total samples were included in this analysis: WTA1, WTB1, and WTC1, which were grown in the presence of thiamine, and WTA2, WTB2, and WTC2, which were grown in the absence of thiamine. 

For the purposes of this workflow, we only performed our analyses on the WTA2 sample, with data from the other five samples provided from other groups working on the project.

# Project Data
Data for this analysis was provided by Dr. Ronda Rolfes. The files we analyzed were reads of *C. albicans* RNA, sequenced with Illumina technology. Our group utilized two data files, WTA2_1.fq.gz and WTA2_2.fq.gz, which represent the forward and reverse reads of the transcripts for the WTA2 sample.

# FastQC
We first analyzed the raw sequence data with FastQC to determine an optimal trimming strategy.

```
$module load fastqc

$fastqc WTA2_1.fq
$fastqc WTA2_2.fq

$module unload fastqc
```

The FastQC analysis showed considerable noise in the first 15 bases of the reads, prompting trimming of these initial bases, as well as removal of Illumina adapter sequences.

See the attached WTA2_1.fastqc.html file for reference.

# Trimming Sequences
We used our FastQC analysis to inform an optimal trimming strategy. The script used for running trimmomatic is contained in the RNA_trim file.

Notable trimming parameters include:

- ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 , which removed the Illumina adapter sequences
- HEADCROP:15 , which removed the first 15 bases of each read
- TRAILING:20 , which removed the last 20 bases of each read
- SLIDINGWINDOW:4:15 , which removed reads containing sequences of 4 bases with an average Phred quality score of 15 or less
- MINLEN:75 , which removed any reads of less than 75 bases


# Post-Trimming FastQC
We again used FastQC to analyze the quality of the reads, this time on the files generated from sequence trimming. The same code was used as above.

See the attached WTA2_1.trPE.fastqc.html file for reference.

This post-trimming analysis showed an improved quality per base, less noise in the beginning of the reads, and successful removal of the adapter sequences.

# Sequence Alignment
With the quality of our sequences improved, we began sequence alignment with the *C. albicans* reference genome, GCF_000182965.3_ASM18296v3_genomic.fna, obtained from NCBI Databases. To perform this alignment, we used bowtie2/2.5.3, first building our indices from the reference genome:

`
$module load bowtie2/2.5.3

$bowtie2-build GCF_000182965.3_ASM18296v3_genomic.fna CalbID

$module unload bowtie2/2.5.3
`

With these indices, we then performed sequence alignment into a .sam file, using the script provided in the attached run_bowtie2 file.

We then converted this .sam file into a .bam file and sorted it withsamtools, using the following code:

`
$module load samtools

$samtools view -S -b WTA2.sam > WTA2.bam
$samtools WTA2.bam -o WTA.srt

$module unload samtools
`


