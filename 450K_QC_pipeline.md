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
Pheno.data <- pData(RG.set)

```

## Quality Control

### 1. Get QC report
```R
## creat QC report
qcReport(RG.set, sampNames = pd$Sample_Name, 
         sampGroups = pd$Sample_Group, pdf = "minfi_qcReport.pdf", maxSamplesPerPage = 6, 
         controls = c())

## creat a methylset
MSet.raw <- preprocessRaw(RG.set) # raw preprocessing
MSet.raw

## qc plot
qc <- getQC(MSet.raw)
pdf("minfi_qcplot.pdf")
plotQC(qc)
dev.off()

## sex prediction
predicted.sex <- getSex(gm.set, cutoff = -2)
pdf("predicted_sex.pdf")
plotSex(predicted.sex)
dev.off()
         
```

### 2. Identify bad quality samples

```R
#get p values of positions
det.p <- detectionP(RG.set)

# get bad quality samples 
bad.samples <- colMeans(det.p > 0.01) > 0.05
bad.sample.names.detP <- colnames(det.p[,bad.samples])
bad.sample.detp.rs <- data.frame(CallRate=1-colMeans(det.p[,bad.samples] > 0.01))
rownames(bad.sample.detp.rs) <- rg.set$Sample_Name[rg.set$Basename %in% bad.sample.names.detP]
write.csv(bad.sample.detp.rs, file=paste(loc.data, "Bad_Samples_CallRate.csv", sep="/"))
save(rg.set, targets, bad.sample.names.detP, bad.samples, file=paste(loc.rgdata, "RG_Channel_Set.Rdata", sep = "/"))

## 
targets <- targets[!bad.samples,]
rg.set <- read.metharray.exp(base=file.path(loc.idat, basename=targets$Slide), targets=targets, recursive=TRUE)
save(RG.set, targets, bad.sample.names.detP, bad.samples, file=paste(loc.rgdata, "RG_Channel_Set_Clean.Rdata", sep = "/"))


```

