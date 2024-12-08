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

## Figure 1: PCA Plot for Differential Gene Expression of Thiamine + and - *C. albicans*
![PCA Plot](https://github.com/user-attachments/assets/fb9634ac-17cf-47c3-b14e-5692cc7859f0 "PCA Plot for Differential Gene Expression of Thiamine + and - *C. albicans*")

## Figure 2: Volcano Plot for Differential Gene Expression of Thiamine + and - *C. albicans*
![Volcano Plot](https://github.com/user-attachments/assets/4ee3d004-af48-4ef2-9d6b-99b4d0ee8cee)


The PCA Plot shows that the first and second principle components explain about 97% of the variance in the data. The Thi+ and Thi- samples are clustered together around similar values of Principal Component 2, at about +10 and -10, respectively. For Principal Component 1, two of the Thi+ and two of the Thi- samples are clustered around a value of +20, whereas one of each of these samples appear to be outliers at negative values of Principal Component 1. These results indicate that, although there is not a great difference in the variance in expression between the two groups, the clustering of the Thi+ and Thi- samples around different values for Principal Component 2 may indicate a small difference in gene expression.

The Volcano Plot affirms these findings. Most of the genes, represented by points on the plot, demonstrate insignificant difference in expression between the two experimental conditions. This can be determined by observing that their p-values fall below the threshold for significance indicated on the plot. Just 14 genes appear to have significantly different expresssion in the thiamine deficient samples, and all of these genes have positive log2-fold changes, meaning that their expression is upregulated in absence of thiamine. 

Given the finding that 14 *C. albicans* genes are upregulated in the absence of thiamine, the next step in our analysis was to determine the identity of these genes to discern any patterns in their function.


# Gene IDs

The 14 genes identified as having significantly altered expression in thiamine deficient conditions were stored in a table in the attached file signif_TH-vTH+.csv (within DESeq2_analysis), using additional code in the attached calb_DESeq_script_FINAL.R file.

The locus tags from this file were stored in another file called signif_geneIDs in the command line. These locus tags were used to parse the NCBI gtf file, linked above, for the gene names and gene IDs of this differentially expressed genes. The following code was used for parsing:

```
grep -wFf signif_geneIDs GCF_000182965.3_ASM18296v3_genomic.gtf|grep "protein_coding"|cut -f9|cut -d ";" -f1,3,5 > signif_gene_annot_info
```

The resultant gene names and IDs are included in the attached file signif_gene_annot_inf, and this information was appended to the signif_TH-vTH+.csv file.

We input these gene names and IDs into the Candida genome page, http://www.candidagenome.org/, and UniProt, https://www.uniprot.org/, to determine the biological functions of these differentially expressed genes. These results were again appended to the signif_TH-vTH+.csv file under the header "Biological Function". This file is linked below for analysis:

![signif_TH-vTH+.csv](https://github.com/user-attachments/files/18053681/signif_TH-vTH%2B.csv).

From parsing the NCBI gtf file, we identified 9 of the 14 genes' associated gene codes. Identification of these 9 genes' functions in UniProt demonstrated that the majority of these genes have functions involved in the thiamine biosynthesis pathway. This result is consistent with out hypothesis that, in a thiamine deficient environment, *C. albicans* will upregulate expression of genes required for endogenous synthesis of thiamine. Thiamine, or vitamin B1, is a crucial cofactor involved in the central metabolic pathway, meaning that its absence could be deleterious to *C albicans* and its ability to produce energy for cellular processes. Therefore, in the absence of environmental thiamine, it makes sense that endogenous thiamine synthesis would be upregulated to compensate for environmental deficits.

It should also be noted that several of these enzymes remain putative, with further structural and functional analysis needed to confirm their identities. Additionally, two genes DUO1 and ERG20, were determined to be involved in mitosis and sterol biosynthesis, respectively. These two growth processes seem somewhat counterintuitve for an organism in nutrient deficient conditions, making their upregulation somewhat unexpected. Additional research should be done to identify their roles in cellular maintenance in thiamine deplete conditions.


# Gene Ontology Enrichment

Given the fact that the majority of the upregulated genes appear to be involved in thiamine biosynthesis, we next sought to determine whether this result is statistically significant, or merely a product of random chance. To do so, we performed a gene ontology enrichment analysis.

To do so, the gene IDs, included in the signif_TH-vTH+.csv. file (linked in previous section), were input into PANTHER classification system (https://www.pantherdb.org/tools/compareToRefList.jsp), using the statistical overrepresentation test. Parameters used include:

- Reference List: Candida albicans
- Annotation Data Set: GO biological process complete
- Test Type: Fisher's Exact
- Correction: Calculate False Discovery Rate

The following output file was generated:

https://www.pantherdb.org/tools/compareToRefList.jsp

This analysis demonstrated that, among the upregulated genes in the thiamine deficient condition, there was a statistically significant number of genes involved in the thiamine biosynthetic pathway as well as the pyridoxal phophate biosynthetic pathways (adjusted p = 1.11 x 10^-07 and 9.19 x 10^-03, respectively). 

This result suggests that when thiamine is absent in its environment *C. albicans* upregulates endogenous production of the vitamin to maintain adequate levels for functioning. Thiamine is an important coenzyme in enzymatic function, with roles in the central metabolic pathway. This makes maintaining thiamine levels, either through the environment or through biosynthesis, a priority for survival.

# Conclusion

This RNASeq analysis demonstrated that in the absence of environmental thiamine, *Candida albicans* upregulates endogenous thiamine biosynthesis pathways to maintain adequate levels of thiamine for metabolic functioning. These results suggest the ability of *C. albicans* to maintain homeostasis during periods of environmental stress and also affirm the crucial role of thiamine to cellular functioning. Future study should examine the mechanisms by which thiamine biosynthesis enzymes are upregulated in *C. albicans*, identifying transcription factors present in higher concentrations in Thi- conditions.

