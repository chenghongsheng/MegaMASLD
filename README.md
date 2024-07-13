# MegaMASLD: An interactive webtool revealing transcriptomic shifts and key molecular events in human MASLD progression

This project contains the R scripts of a mega-analysis of metabolic dysfunction-associated steatotic liver disease (MASLD) liver transcriptomes (MegaMASLD). The analysis used publicly accessible raw RNAseq data of over 800 livers and analysed them in a standardized and integrative manner, aiming to unravel targetable molecular events in MASLD. 
The manuscript is currently under review.

## Manuscript 
Title: MegaMASLD: An interactive webtool revealing transcriptomic shifts and key molecular events in human MASLD progression

DOI: pending

Link: pending

R shiny tool: https://bioanalytics-hs.shinyapps.io/MegaMASLD/

## Data
Hepatic RNAseq datasets used in this study are GSE105127, GSE107650, GSE115193, GSE126848, GSE130970, GSE147304, GSE160016, GSE162694, GSE167523, GSE173735, GSE174478, GSE185051, GSE192959, GSE193066, GSE193080, GSE207310, and GSE213621.

## R scripts
Data.integration_Deseq2_Boruta.R contains the R script for batch and gender correction to integrate the RNAseq data, differential expression analysis between MASLD stages and key feature selection using Boruta algorithm.

ssGSEA_ROC.staging.R contains the R script 

Gender_MASLD.interaction.R

cNMF.R

Pseudotime.R

