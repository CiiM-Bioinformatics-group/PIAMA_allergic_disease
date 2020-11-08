# Quality Control of 450K DNA methylation data

Reference: <br>
http://bioconductor.org/packages/release/bioc/vignettes/minfi/inst/doc/minfi.html <br>
https://www.bioconductor.org/packages/release/bioc/vignettes/wateRmelon/inst/doc/wateRmelon.pdf <br>
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3592906/ <br>

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

### 2. Filter samples

* 2.1. Check Gender Concordance
```R
## get gset
ratioSet <- ratioConvert(MSet.raw, what = "both", keepCN = TRUE)
gset <- mapToGenome(ratioSet)

## sex prediction
predictedSex <- getSex(gset, cutoff = -2)$predictedSex
gset <- addSex(gset, sex=predictedSex)

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

* 2.2. Identify bad quality samples

```R
#get p values of positions
det.p <- detectionP(RG.set)

# get bad quality samples 
failed.probes <- det.p > 0.01

#fraction of failed positions per sample 
failed.fraction <- colMeans(failed.probes)
failed.fraction <- data.frame(failed.fraction)
failed.fraction <- cbind(Pheno.data$Sample_Name,failed.fraction)

#get bad quality samples (samples with more than 0ne percent failed positions)
failed.fraction <- data.frame(failed.fraction[failed.fraction[,2] > 0.01,])

fail.detp<-failed.fraction$Pheno.data.Sample_Name
write.csv(fail.detp, file="failed_samples_by_probe_frequency.csv")

```

* 2.3. Remove bad quality samples and sex mismatch samples

```R
## remove all failed samples from RGset
fail.sample<-unique(c(sex.mismatch,fail.detp))
fail.sample.index<-which(Pheno.data$Sample_Name %in% fail.sample)
filtered_RG.set <- RG.set[,-bad.sample.index]
save(filtered_RG.set,file="filename")
```

### 3. Filter probes

* 3.1 Filter on detection P values

```R
m.set <- preprocessRaw(filtered_RG.set)
ratioSet <- ratioConvert(m.set, what = "both", keepCN = TRUE)
gm.set <- mapToGenome(ratioSet)

det.p <- detectionP(filtered_RG.set)
bad.probes <- rowMeans(det.p > 0.01) > 0.1
table(bad.probes)
bad.probe.names.detP <- rownames(det.p[bad.probes,])

```

* 3.2 Filter targets from multiple genomic sites <br>
Reference: From supplementary material Chen et al. Epigenetics 2013 <br>
http://www.sickkids.ca/Research/Weksberg-Lab/Publications/index.html. 

```R
nonspec <- read.csv("Chen_nonspecific_probes_450k.csv", header = TRUE, sep = ";")
bad.probe.names.cr <- as.character(nonspec$TargetID[which((nonspec$bm50 + nonspec$bm49 + nonspec$bm48 + nonspec$bm47) > 0)]) 

```

* 3.3 Drop probes with SNP at interrogation or extension site

```R
## get gset

gset.nosnp <- dropLociWithSnps(gm.set, snps = c("CpG", "SBE"), maf = 0.05)
all.probes <- rownames(gm.set)
nosnp.probes <- rownames(gset.nosnp)
bad.probe.names.pm <- all.probes[!(all.probes %in% nosnp.probes)]

failed.probes <- unique(c(bad.probe.names.detP,bad.probe.names.cr,bad.probe.names.pm))
save(failed.probes,file="failed_probes.Rdata")

```

* 3.4 Remove probes on X and Y chr

```R
print("get probes on X,Y")
keepX <- featureNames(gm.set) %in% ann450k$Name[ann450k$chr %in% c("chrX")]
keepY <- featureNames(gm.set) %in% ann450k$Name[ann450k$chr %in% c("chrY")]

probe.names.x <- rownames(gm.set[keepX,])
probe.names.y <- rownames(gm.set[keepY,])

print("save probes on X,Y")
save(probe.names.x,probe.names.y, file="probes_XY.Rdata")

remove.probe<-unique(c(failed.probes,probe.names.x,probe.names.y))
save(remove.probe, file="remove_probes.Rdata")
```

### 4. Normalization 

```R
## probes filtering

MSet.raw <- preprocessRaw(filtered_RG.set)
r.probes <- which(rownames(MSet.raw) %in% remove.probe)
MSet.rp <- MSet.raw[-r.probes,]

## dasen normalization (P pacakge wateRmelon)
bn.dasen <- dasen(MSet.rp)
save(bn.dasen,"dasen_betas.RData")

## can also choose other normalization methods, such as "Funnorm", "ssNoob", "Quantile" from minfi package
## for example Quantile normalization
m.set.sq <- preprocessQuantile(filtered_RG.set)
m.set.flt <- m.set.sq[!(rownames(m.set.sq) %in% remove.probe),]

save(m.set.flt, file=”xxx")
M.val <- getM(m.set.flt) ## get M value
B.val <- getBeta(m.set.flt) ## get beta value

## DMSplot
phenoData <- pData(filtered_RG.set)
mdsPlot(B.val, numPositions = 1000, sampGroups = phenoData$status,sampNames = phenoData$Sample.Name)

```

### 5. estimate cell counts

```R
counts <- estimateCellCounts(RG.set,compositeCellType="Blood", cellTypes = c("CD8T","CD4T", "NK","Bcell","Mono","Eos","Neu"), meanPlot = FALSE)

```


