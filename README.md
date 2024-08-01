# RNA-Seq data analysis of Melanoma and Melanocytes

This repository contains the data analysis workflow for RNA-Seq data developed by Daniel Zimmermann and used during the master's thesis "Molecular and Cell Biological Analysis of the Oncogenesis of Melanoma" by Magdalena Hohlrieder.

## Required data

-   Number of reads per gene, from `feature_count`
-   Sample metadata table
-   .gtf file for gene annotation

## Differential Expression Analysis

In `DESeq_analysis.Rmd`, differentially expressed genes are identified using the `DESeq2` package. For each comparison, both the complete gene list and the list of significant genes are saved for further analysis.

## Functional Analysis

In `functional_analysis.Rmd`, significantly up- or downregulated pathways are identified by Gene Set Enrichment Analysis (GSEA) using the `clusterProfiler` package. Specifically, the databases used are - Gene Ontology - Biological Processes - Molecular Function - Cellular Component - KEGG Pathways - Disease Ontology
