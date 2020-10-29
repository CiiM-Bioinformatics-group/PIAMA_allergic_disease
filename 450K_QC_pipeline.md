# Quality Control of 450K DNA methylation data

## Preparation

### 1. Install R packages
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("minfi")
BiocManager::install("IlluminaHumanMethylation450kmanifest") ## or "IlluminaHumanMethylationEPICmanifest" for EPIC array
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19") ## or "IlluminaHumanMethylationEPICanno.ilm10b2.hg19" for EPIC array
BiocManager::install("wateRmelon")
BiocManager::install("GO.db")
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")
install.packages(c("ggplots2", "ggfortify", "RColorBrewer", "gridExtra", "pheatmap", "sva", "MASS", "compare", "tableone", "matrixStats", "plyr", "qqman", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"))
```

### 2. Input Data
* raw idat data
* sample sheet 
