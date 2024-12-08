# Master Notes for the RNASeq Analysis Pipeline
This file contains notes about the workflow for this RNASeq analysis of *C. albicans*. 

# Project Objectives
This project seeks to determine the differential gene expression of the human commensal yeast, *Candida albicans*, when grown in the presence and absence of thiamine. Six total samples were included in this analysis: WTA1, WTB1, and WTC1, which were grown in the presence of thiamine, and WTA2, WTB2, and WTC2, which were grown in the absence of thiamine. 

For the purposes of this workflow, we only performed our analyses on the WTA2 sample, with analysis of the other five samples provided by other groups working on the project.

# Project Data
Data for this analysis was provided by Dr. Ronda Rolfes. The files we analyzed were reads of *C. albicans* RNA, sequenced with Illumina technology. Our group utilized two data files, WTA2_1.fq.gz and WTA2_2.fq.gz, which represent the forward and reverse reads of the transcripts for the WTA2 sample.

The reference genome we used for alignment, GCF_000182965.3, was obtained from NCBI database at the following link:

https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/943/734/735/GCF_943734735.2_idAnoGambNW_F1_1/GCF_943734735.2_idAnoGambNW_F1_1_genomic.gff.gz.

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
We again used FastQC to analyze the quality of the reads, this time on the files generated from sequence trimming. The same code was used as above in the initial FastQC section.

See the attached WTA2_1.trPE.fastqc.html file for reference.

The results of the trimming strategy for all six samples are compiled in the following document:

https://docs.google.com/spreadsheets/d/1AOa-XaTzR_PKMIRQDmu8oDTmawXXnkIwEjKOQkNC7Vs/edit?gid=0#gid=0

Prior to trimming, there were 21972519 reads, with poor per-base sequence content, high sequence duplication, and adapter contamination.

After trimming, there were 21012870 reads, with good per-base sequence content, low sequence duplication, and no adapter contamination.

This represents a loss of 959649 reads, with 95.63% reads retained, in addition to a significant improvement in per-base read quality, reduction in duplicate sequences, and removal of adapter contamination.

# Sequence Alignment
With the quality of our sequences improved, we began sequence alignment with the *C. albicans* reference genome, GCF_000182965.3_ASM18296v3_genomic.fna, obtained from NCBI Databases. 

The aboslute path to this file is:

https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/943/734/735/GCF_943734735.2_idAnoGambNW_F1_1/GCF_943734735.2_idAnoGambNW_F1_1_genomic.gff.gz. 

To perform this alignment, we used bowtie2/2.5.3, first building our indices from the reference genome:

```
$module load bowtie2/2.5.3

$bowtie2-build GCF_000182965.3_ASM18296v3_genomic.fna CalbID

$module unload bowtie2/2.5.3
```

With these indices, we then performed sequence alignment into a .sam file, using the slurm script provided in the attached run_bowtie2 file.

We then converted this .sam file into a .bam file and sorted it with samtools, using the following code:

```
$module load samtools

$samtools view -S -b WTA2.sam > WTA2.bam
$samtools WTA2.bam -o WTA.srt

$module unload samtools
```

Sequence alignment results for all six samples are compiled in the following document:

https://docs.google.com/spreadsheets/d/1fa-FXVMlCXOZkbHSx_mMg0OXLMy9BeBJg8uWrEMpKGo/edit?gid=0#gid=0

For the WTA2 sample, all 21012870 reads left over from trimming were detected by bowtie, and the overall alignment rate for the sample was 98.45%. Across all six samples included in the analysis, overall alignment rates were approximately 98%, which represents high overall alignment for all included samples.

# High-throughput Sequence Analysis
After sequence alignment, we next counted the number of reads per gene using the software HTSeq. This required creating a conda environment in the local terminal with the following code:

```
$module load anaconda3
$conda create --name htseq
$conda init
$module load anaconda3
$conda activate htseq
$conda install -c bioconda htseq
```

The slurm script for running this HTSeq analysis is included in the attached file run_htseq. The output files are included in the attached htseq_counts file, and the resulting counts of reads per gene were stored in the output file htseq_counts.txt.

The following reads were not able to be assigned to a specific gene:

- __no_feature    9820912
    - reads that could not be assigned to any gene
- __ambiguous     18124
    - reads that aligned to multiple genes
- __too_low_aQual 1114728
    - reads with alignment quality lower than 10
- __not_aligned   156874
    - reads with no alignment in the SAM file

# Differential Expression Analysis
We next performed a differential expression analysis with the R Software DESeq2.The code used for this analysis is included in the attached file, calb_DESeq_script_FINAL.R.

The output files of this analysis are included in the attached file DESeq2_analysis, and are depicted below:


![PCA Plot](https://github.com/user-attachments/assets/fb9634ac-17cf-47c3-b14e-5692cc7859f0 "PCA Plot for Differential Gene Expression of Thiamine + and - *C. albicans*")

![R_volcano_plot_correct.pdf](https://github.com/user-attachments/files/18053611/R_volcano_plot_correct.pdf)







# Gene IDs

The 14 genes identified as having significantly altered expression in thiamine deficient conditions were stored in a table in the attached file signif_TH-vTH+.csv (within DESeq2_analysis).

The locus tags from this file were stored in another file called signif_geneIDs in the command line. These locus tags were used to parse the NCBI gtf file for the gene names and gene IDs of this differentially expressed genes, using the following code:

```
grep -wFf signif_geneIDs GCF_000182965.3_ASM18296v3_genomic.gtf|grep "protein_coding"|cut -f9|cut -d ";" -f1,3,5 > signif_gene_annot_info
```

The resultant gene names and IDs are included in the attached file signif_gene_annot_inf, and this information was appended to the signif_TH-vTH+.csv file.

We input these gene names and IDs into the Candida genome page, http://www.candidagenome.org/, and UniProt, https://www.uniprot.org/, to determine the biological functions of these differentially expressed genes. These results were again appended to the signif_TH-vTH+.csv file.

# Gene Ontology Enrichment

The gene IDs, parsed from the NCBI gtf file, were input into PANTHER classification system (https://www.pantherdb.org/tools/compareToRefList.jsp) to perform a gene ontology enrichment analysis.

This analysis demonstrated that, among the upregulated genes in the thiamine deficient condition, there was a statistically significant number of genes involved in the thiamine biosynthetic pathway as well as the pyridoxal phophate biosynthetic pathways (adjusted p = 1.11 x 10^-07 and 9.19 x 10^-03, respectively). 

This result suggests that when thiamine is absent in its environment *C. albicans* upregulates its own production of thiamine to maintain adequate levels for functioning. Thiamine is an important cofactor in enzymatic function, with roles in the central metabolic pathway. This makes maintaining thiamine levels, either through the environment or through biosynthesis, a priority for survival.

# Discussion

