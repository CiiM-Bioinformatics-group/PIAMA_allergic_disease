# General EWAS pipeline

## Remove outliers using the IQR*3 (Tukey) method

```R
## load DNA methylation beta values
load("filename")

# transfer beta value to M value
M.brush<-log2(beta.brush/(1-beta.brush))
save(M.brush,file="PIAMA_brush_M.Rdata")

# trimming, remove outliers using the IQR*3 (Tukey) method
removeOutliers<-function(probes){
  require(matrixStats)
  if(nrow(probes) < ncol(probes)) warning("expecting probes are rows (long dataset)")
  rowIQR <- rowIQRs(probes, na.rm = T)
  row2575 <- rowQuantiles(probes, probs = c(0.25, 0.75), na.rm = T)
  maskL <- probes < row2575[,1] - 3 * rowIQR 
  maskU <- probes > row2575[,2] + 3 * rowIQR 
  initial_NAs<-rowSums(is.na(probes))
  probes[maskL] <- NA
  removed_lower <- rowSums(is.na(probes))-initial_NAs
  probes[maskU] <- NA
  removed_upper <- rowSums(is.na(probes))-removed_lower-initial_NAs
  N_for_probe<-rowSums(!is.na(probes))
  Log<-data.frame(initial_NAs,removed_lower,removed_upper,N_for_probe)
  return(list(probes, Log))
}

system.time(OutlierResults<-removeOutliers(M.brush))
M.brush2<-OutlierResults[[1]]
Log<-OutlierResults[[2]]

# save M values after trimming
save(M.brush2,file="PIAMA_brush_trimmed_M.Rdata")
save(Log,file="Outlier_log_M.Rdata")
```

## Generate surrogate variables (for tissue other than blood)

```R
library(sva)
library(tableone)
library(matrixStats)

## load trimmed M values
load("M_matrix.Rdata")
PHENO<-read.csv("phenotype.csv")
## samples in M_matrix and PHENO should match in orders

k = which(is.na(M_matrix), arr.ind=TRUE)
M_matrix[k] = rowMedians(M_matrix, na.rm=TRUE)[k[,1]]

mod1<-model.matrix(~.,PHENO)
mod0<-model.matrix(~.,PHENO[,-1])
n.sv <- num.sv(dat =M_matrix, mod = mod1, method = "be", vfilter = 2000, B = 1000, seed =1)
svobj1= sva(M_matrix,mod1,mod0,n.sv=n.sv)
SVs = as.data.frame(svobj1$sv)
colnames(SVs) <-paste0("sv",1:ncol(SVs))
modSv1 = cbind(PHENO,SVs)
save(modSV1,file="pheno_sva.Rdata")

```

## Fit logistic regression model
```R
# load libraries
library(foreign)
library(data.table) # to process results
library(MASS) # rlm function for robust linear regression
library(sandwich) # estimation of the standard error
library(lmtest) # to use coeftest
library(parallel) # to use multicore approach - part of base R
library(R.utils)
library(matrixStats)
library(plyr) 

## load trimmed M_matrix and phenotype
load("M_matrix.Rdata")
PHENO<-read.csv("phenotype.csv")
## samples in M_matrix and PHENO should match in orders

## function to perform logistic regression model
## Y: disease/phenotype, X1...Xn: covariates

GLMtest = function(meth_matrix,methcol,Y, X1, X2, X3, X4) {
  mod = glm(Y~meth_matrix[, methcol] +X1+X2+X3+X4,family = "binomial")
  cf = summary(mod)$coefficients
  cf[2, c("Estimate", "Std. Error", "Pr(>|z|)")]  
}

system.time(ind.res <- mclapply(setNames(seq_len(ncol(M_matrix)), dimnames(M_matrix)[[2]]), GLMtest, meth_matrix=M_matrix, Y=PHENO[,1], X1=PHENO[,2], X2=PHENO[,3], X3=PHENO[,4],X4<- PHENO[,5]))

# summary results
all.results<-ldply(ind.res,rbind)
names(all.results)<-c("probeID","BETA","SE", "P_VAL")

# add sample size for each probe
M_matrix<-t(M_matrix)
all.results<-all.results[match(rownames(M_matrix),all.results$probeID),]
all.results$N<-rowSums(!is.na(M_matrix))

# calculate lambda value
median(qchisq(as.numeric(as.character(all.results$P_VAL)),df=1,lower.tail = F),na.rm = T)/qchisq(0.5,1)

# save results
write.table(all.results, filename,na="NA")
gzip(filename)
```

## QQ plot and Manhattan Plot 

```R
library(qqman)
library(ggplot2)
library(MASS)
library(GWASTools)
library(ggrepel)
library(dplyr)

## load EWAS results and illumina annotation files (can be downloaded from illumina website)
results<-read.table("PIAMA_16_allergy_brush_sva_model2_be.gz")
anno<-read.csv("~/Documents/Projects/PIAMA/HumanMethylation450.csv",stringsAsFactors = F)

# add fdr and find sig
results<-results[order(results$P_VAL),]
pvalue<-results$P_VAL
FDR<-p.adjust(results$P_VAL,"fdr")
results<-cbind(results,FDR)

# add annotation
index<-match(results$probeID,anno$IlmnID)
anno1<-anno[index,]
results<-cbind(results,anno1)

sig<-results[(which(results$FDR<0.05)),]
write.table(sig,"PIAMA_brush_allergy_sig.txt",row.names = F)

########################  QQ plot  ########################
pvalue<-results$P_VAL
tiff("qqplot_brush_allergy.tiff",width=1500, height=1500,res=300)
qq(pvalue)
t <-median(qchisq(as.numeric(as.character(pvalue)),df=1,lower.tail = F),na.rm = T)/qchisq(0.5,1)
lamda<- round(t,digit=3)
text (2,5, paste0("lambda=",lamda),cex = 1)
dev.off()

########################  manhattan plot  ########################
# rename the data frame
man<-results[,c("probeID","CHR","MAPINFO","P_VAL")]
colnames(man)<-c("SNP","CHR","BP","P")
man<-droplevels(man)
man$SNP<-as.character(man$SNP)
man$CHR<-as.integer(man$CHR)


manhattan(man,suggestiveline = fdr,genomewideline = F,col = c("blue4", "orange3"),highlight=snpsOfInterest)

```
