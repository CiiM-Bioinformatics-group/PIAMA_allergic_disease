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

### 3. Load data
```R
library(minfi)

# set the locations of important files
loc <- "/data/f114798/"
loc.idat <- "Data/ImageData"
loc.sheet <- "Samplesheets"
loc.rgdata <- "Data/RG_data"
loc.data <- "Data"
setwd(loc)

# conversion and storing of idat to rg files 
targets <- read.metharray.sheet(loc.sheet, pattern="*.csv", recursive=TRUE)
targets$Basename <- apply(targets, 1, function(target) {paste(target[6], target[7], sep="_")})
RG.set <- read.metharray.exp(base=file.path(loc.idat, basename=targets$Slide), targets=targets, recursive=TRUE)
save(RG.set, file="filename")

```

## Quality Control

### 1. Identify bad quality samples

```R
#get p values of positions
p.values.probes <- detectionP(RG.set)

#get failed positions
failed.probes <- p.values.probes > 0.01

#fraction of failed positions per sample 
failed.fraction <- colMeans(failed.probes)
failed.fraction <- data.frame(failed.fraction)
failed.fraction <- cbind(pheno.data$Sample.ID,failed.fraction)

#get bad quality samples (samples with more than 0ne percent failed positions
failed.fraction <- data.frame(failed.fraction[failed.fraction[,2] > 0.01,])

#save data
save(failed.fraction, file = "filename")

```





