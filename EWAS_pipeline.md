# General EWAS pipeline

## Remove outliers using the IQR*3 (Tukey) method

## Generate surrogate variables (for tissue other than blood)

```R
load("M_values.Rdata")
PHENO<-read.csv("phenotype.csv")

mod1<-model.matrix(~.,PHENO)
mod0<-model.matrix(~.,PHENO[,-1])
n.sv <- num.sv(dat =M_matrix, mod = mod1, method = "be", vfilter = 2000, B = 1000, seed =1)
svobj1= sva(M_matrix,mod1,mod0,n.sv=n.sv)
SVs = as.data.frame(svobj1$sv)
colnames(SVs) <-paste0("sv",1:ncol(SVs))
modSv1 = cbind(PHENO,SVs)
```

## Fit logistic regression model

## QQ plot and Manhattan Plot 
