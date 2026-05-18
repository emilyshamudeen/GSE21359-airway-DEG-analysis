#   Differential-Gene-Expression-Analysis-of-Human-Airway-Epithelium

## Introduction 
This R script shows the downloading of GSE data from the NCBI GEO database and how to obain differenially expressed genes from microarray data.

## GSE21359
The microarray dataset being analysed is GSE21359. It contains expression profiles of small airway epithelium obtained from three groups of subjects:
1. Healthy non-smokers (n = 53)
2. Healthy smokers (n = 59)
3. Smokers with COPD (n = 23)

The primary gene of interes in this study is CXCL14, a chemokine expressed in various epithelia. The study hypothesized that the airway epithelium responds to cigarette smoking with altered CXCL14 expression as a disease relevant molecular phenotype. 

Full dataset details on NCBI GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21359

## What this script does 
1. Import the data
The GEOquery package is used to download the dataset directly from GEO. The data is stored as an ExpressionSet object containing the expression matrix, sample metadata, and gene annotation.

2. Check normalisation and scale
The data processing method is inspected to confirm whether log2 transformation is needed. A boxplot is drawn to confirm that all sample distributions are similar, indicating the data has been normalised.

3. Inspect the clinical variables
The sample metadata is explored to identify the columns needed for analysis. A group variable is created using stringr to classify each sample as non_smoker, healthy_smoker, or copd_smoker based on the sample title. Sex and age are also extracted as covariates.

4. Sample clustering and PCA
Two quality control visualisations are produced. A correlation heatmap shows how similar samples are to each other — samples in the same biological group should cluster together. A PCA plot compresses the 54,675-gene expression matrix down to two dimensions to check whether the three groups separate as expected.

5. Differential expression analysis
The limma package is used to perform differential expression analysis. Lowly-expressed probes are filtered out to reduce false positives. Array weights are calculated to downweight any poor-quality samples. A linear model is fitted and three pairwise contrasts are tested:
- Healthy smoker vs non-smoker: Effect of smoking alone
- COPD vs non-smoker: Combined effect of smoking and disease
- COPD vs healthy smoker: Disease effect within smokers
The eBayes function applies empirical Bayes shrinkage to stabilise variance estimates, and decideTests summarises the number of differentially expressed genes per contrast.

6. Gene annotation and volcano plots
Gene annotation is retrieved using fData and merged into the results table. Volcano plots are produced for each contrast, with adjusted p-value on the Y-axis and log2 fold change on the X-axis. Significant genes (adj.P < 0.05, |logFC| > 1) are highlighted.

7. Gene of interest lookup
Results for specific genes of interest are extracted by GenBank accession ID. CXCL14 (NM_004887), the primary study gene, and NKX3-1 (NM_006167), a lung epithelial transcription factor, are checked across all three contrasts.

## PCA Plot 
<img width="1350" height="900" alt="image" src="https://github.com/user-attachments/assets/f81bbcfe-078c-4b00-949f-abf98598b820" />

## Correlation Heatmap
<img width="2700" height="2400" alt="image" src="https://github.com/user-attachments/assets/bc7bb063-52fb-425a-aa51-59ef5b6ef010" />

## Volcano Plot: Healthy Smoker vs Non-Smoker
<img width="1200" height="900" alt="image" src="https://github.com/user-attachments/assets/3ea2878d-d594-49d7-bb0a-a3118936facf" />

## Volcano Plot — COPD vs Non-Smoker
<img width="1200" height="900" alt="image" src="https://github.com/user-attachments/assets/22d8fa6d-cab5-406f-8fb8-025e20e1414a" />

## Volcano Plot — COPD vs Healthy Smoker 
<img width="1200" height="900" alt="image" src="https://github.com/user-attachments/assets/824ae5e3-88d3-4c8d-9739-852dd69ab678" />

## R packages used 
- GEOquery — download data from NCBI GEO
- limma — differential expression analysis
- ggplot2 — data visualisation
- ggrepel — label positioning in plots
- pheatmap — correlation heatmap
- dplyr — data manipulation
- stringr — string detection for group assignment
