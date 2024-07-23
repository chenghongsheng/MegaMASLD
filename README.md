# MegaMASLD: An interactive platform for exploring stratified transcriptomic signatures in MASLD progression

This project contains the R scripts of a mega-analysis of metabolic dysfunction-associated steatotic liver disease (MASLD) liver transcriptomes (MegaMASLD). The analysis used publicly accessible raw RNAseq data of over 800 livers and analysed them in a standardized and integrative manner, aiming to unravel targetable molecular events in MASLD. 
The manuscript is currently under review.

## Manuscript 
Title: MegaMASLD: An interactive platform for exploring stratified transcriptomic signatures in MASLD progression

DOI: 10.1101/2024.07.21.603199

Link to preprint: https://www.biorxiv.org/content/10.1101/2024.07.21.603199v1

R shiny tool: https://bioanalytics-hs.shinyapps.io/MegaMASLD/

## Data
Hepatic RNAseq datasets used in this study are GSE105127, GSE107650, GSE115193, GSE126848, GSE130970, GSE147304, GSE160016, GSE162694, GSE167523, GSE173735, GSE174478, GSE185051, GSE192959, GSE193066, GSE193080, GSE207310, and GSE213621.

## R scripts
**Data.integration_Deseq2_Boruta.R** contains the R script for batch and gender correction to integrate the RNAseq data, differential expression analysis between MASLD stages and key feature selection using Boruta algorithm.

**ssGSEA_ROC.staging.R** contains the R script for quantitative analysis of disease severity based on human-mouse harmonized MASLD Human-mouse harmonized MASLD signature and ROC analysis to determine the accuracy of the transcriptomic-based staging.

**Gender_MASLD.interaction.R** describes the differential expression analysis to study the interaction between gender and MASLD progression.

**cNMF.R** contains the R script for an unsupervised clustering of the MASLD transcriptomes using consensus non-negative matrix factorization (cNMF). The analysis enables a histologic-independeint risk stratication of the patients based on the liver transcriptomes.

**Pseudotime.R** contains the R script to predict the pseudotemporal projection of the disease progression based on the integrated liver transcriptomes. It also computes a take-off value, defined as a pseudotime point when the expression of a gene begins to exhibit a significant continuous upward trend, of every gene. 

