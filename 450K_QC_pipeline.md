# Quality Control of 450K DNA methylation data

## Preparation

### 1. Install and load R packages
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

## also load the R packages

```

### 2. Input Data
* raw idat data
* sample sheet 

### 3. Load data
```R
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
qcReport(RG.set, sampNames = Pheno.data$Sample_Name, 
         sampGroups = Pheno.data$Sample_Group, pdf = "minfi_qcReport.pdf", maxSamplesPerPage = 6, 
         controls = c())

## creat a methylset
MSet.raw <- preprocessRaw(RG.set) # raw preprocessing
MSet.raw

## qc plot
qc <- getQC(MSet.raw)
pdf("minfi_qcplot.pdf")
plotQC(qc)
dev.off()
         
```

### 2. Check Gender Concordance
```R
## get gset
ratioSet <- ratioConvert(MSet.raw, what = "both", keepCN = TRUE)
gset <- mapToGenome(ratioSet)

## sex prediction
predictedSex <- getSex(gset, cutoff = -2)$predictedSex
gset <- addSex(gset, sex=predictedSex)

## 
pheno.order <- data.frame(Pheno.data[,c(2,3)])
pheno.order.plus.sex <- cbind(pheno.order,predictedSex)
save(pheno.order.plus.sex, file = "../../Data/pheno_order_plus_sex.R")

documentedSex <-Pheno.data$gender
levels(documentedSex) <- c("F","M") ## change the sex markers if it is not F and M in your origin data

gcheck <- function(docsex,predsex,samnam) {
  docsex <- as.character(docsex)
  sex.match <- identical(docsex,predsex)
  if (sex.match == FALSE) {
    sex.mismatch <- samnam[which(docsex != predsex)]
    write.csv(sex.mismatch, file="sexmismatch.csv") }
}
gcheck(documentedSex,predictedSex, Pheno.data$Sample_Name)

pdf("sex_concordance_plot.pdf")
plotSex(getSex(gset, cutoff = -2))
dev.off()

```

### 3. Identify bad quality samples

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

