# PIAMA_allergic_disease
Rationale: Epigenetic signatures in the nasal epithelium, which is a primary
interface with the environment and an accessible proxy for the bronchial
epithelium, might provide insights into mechanisms of allergic disease. We aimed
to identify and interpret methylation signatures in nasal epithelial brushes
associated with rhinitis and asthma.

Cancan Qi, email: tracyqican@gmail.com

## Content of the scripts used in this project

### 1. The scripts of main anlayses (folder "Main")
* S1. EWAS of asthma, rhinitis and allergic disease;
* S2. meta analysis of discovery and replication;
* S3. eQTM analysis of replicated CpGs, and pathway analysis;
* S4. prediction of allergic disease;
* S5. single cell analysis of 2 asthma patients and 2 controls;
* S6. correlation with environmental factors;
* S7. DMR analysis;
* S8. IgE stratified analysis.
* S9. Other test

### 2. The scriptes of main plots (folder "Plots")
* P1. Manhattan plot and QQ plot of 3 phenotypes;
* P2. Boxplot of important CpGs; Boxplot of different time windows;
* P3. Venn diagram of 3 phenotypes;
* P4. Ohter useful plots.

### 3. The scripts for Data management (folder "Phenotype")
* D1 phenotype of allergic disease
* D2 phenotype of environmental factors
* D3 phenotype of lung function


## General pipelines for DNA methylation 450K data

### 1. Quality control (450K_QC_pipeline)
* QC of samples
* QC of probes
* Data normalization

### 2. EWAS (logistic regression model) (EWAS_pipeline)
* Data trimming
* SVA
* Fit logistic regression model
* QQplot and Manhattan plot



