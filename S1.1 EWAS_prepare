## This script is to
## 1. transfer beta value to M value
## 2. trimming
## performed on peregrine 19/02/2018

rm(list=ls())

setwd("/home/p282473/projects/PIAMA/data")
load("PIAMA_brush_corrected.Rdata")
# remove duplo samples
beta.brush<-beta.brush[,-grep("duplo",colnames(beta.brush))]

# transfer beta value to M value
M.brush<-log2(beta.brush/(1-beta.brush))
save(M.brush,file="PIAMA_brush_M.Rdata")

# trimming
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


