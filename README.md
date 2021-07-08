# Coregulation-of-complexes-in-aneuploidtumors

This is repository for the codes to run main analyses and produce data for the study "Regulation of protein complex partners as a compensatory mechanism in aneuploid tumors"

Authors: Gökçe Senger, Stefano Santaguida, Martin H. Schaefer

## Summary 

Aneuploidy, defined as whole chromosomal or chromosome arm-level changes, is a hallmark of human tumors, but its role in cancer still remains to be fully elucidated. We studied transcriptomic and proteomic changes induced by whole chromosome-level amplifications in 9830 TCGA and corresponding 298 CPTAC tumor samples by integrating cancer aneuploidy, transcriptomic and proteomic data. We found that only a relatively small number of genes on the amplified chromosomes changed expression while comparably many changes happened on other chromosomes. Moreover, we found an association between those changes on other chromosomes and co-complex members of proteins from amplified chromosomes. Integrating co-abundance analyses between aneuploid proteins and their co-complex members on other chromosomes with protein features revealed that aggregation-prone aneuploid proteins and those involved in a smaller number of complexes are related with stronger correlations with their partners. These findings highlight the importance of compensation for stoichiometric imbalance in protein complexes and suggest that avoiding proteotoxicity of unpaired complex members might allow cancer cells to deal with the excess amount of expression changes induced by aneuploidy. 

## Files

Note: Publicy available data used in this study can be downloaded from the corresponding consortium (Please see "Data availability" statement).

AneuploidyScores.R - R code to calculate whole chromosome-level aneuploidy scores and to detect cancer type-specific aneuploidies. 

DataProcessing.R - Transcriptomic and proteomic data processing.

human_genes.txt - Ensembl gene IDs and gene symbols obtained from ensembl BioMart (Human genome version GRCh38.p13 - downloaded on May, 2019). 

DetectingChanges.R - Detecting transcriptomic and proteomic changes induced by amplifications/deletions

AneuploidyTable.RData - Whole chromosome-level aneuploidy scores for 10522 TCGA samples. This file can also be produced by running the AneuploidyScores.R. 

CasesSelected_alltumours_amp.xlsx - Detected 86 whole chromosome-level amplifications. 

AssociationTests.R - Chi-square tests between changes on aneuploid chromosomes and those on other chromosomes.

HumanPCgenes_Subunits_allchr.RData - CORUM human known protein complexes and subunits. 

Table_cutoff.xlsx - File for which cut-off were used to detect DEGs in each aneuploidy case. 

DE_Transcriptome_Sig.RData - Transcriptomic changes for 86 detected cancer type-specific amplifications with annotation for being differentially expressed or not. 

DE_Proteome_Sig.RData - Proteomic changes for 13 detected cancer type-specific amplifications with annotation for being differentially abundant or not. 

Correlations.R - Calculation of spearman correlations for all possible pairwise correlations. 

ComplexNumber.R - Finding number of complexes a protein involved in (for differentially abundant proteins on aneuploid chromosomes).

SpearmanCorr_mrna&protein.RData - Spearman correlations between up-regulated proteins of aneuploid chromosomes and their co-complex members on other chromosomes. 



