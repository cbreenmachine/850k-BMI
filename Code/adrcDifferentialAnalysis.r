##### Differential analysis of ADRC microarray data

### Load packages anda data for differential analysis

rm(list=ls())
getwd()
set.seed(1234)
library(limma)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(RColorBrewer)
library(missMethyl)
library(Gviz)
library(stringr)
library(NHMMfdr)
library(sva)
library(data.table)

load("adrcPreprocessedData.RData")

### set variables for making models
# only groupwise model shown below


#--> For smoking
# smoking <- targets$pack_years
# ix <- !is.na(smoking)
# smoking <- smoking[ix]
# targets <- targets[ix, ]
# mValues <- mValues[, targets$patient_id]

bmi <- targets$bmi

gran <- as.numeric(targets$Gran)
mono <- as.numeric(targets$Mono)
nk <- as.numeric(targets$NK)
bcell <- as.numeric(targets$Bcell)
cd8 <- as.numeric(targets$CD8T)
cd4 <- as.numeric(targets$CD4T)
sex <- as.factor(targets$sex)
age <- as.numeric(as.character(targets$age))
group <- as.factor(targets$cohort)
slide <- as.factor(targets$chip_id)
### generate model and identify surrogate variables

mod <- model.matrix( ~ bmi + group + sex + age + slide + gran + nk + mono + bcell + cd8 + cd4)
mod0 <- model.matrix(~ group + sex + age + slide + gran + nk + mono + bcell + cd8 + cd4)

n.sv <- num.sv(mValues,mod,method='leek')
svobj <- sva(mValues,mod,mod0,n.sv=n.sv)
modSv <- cbind(mod,svobj$sv)

### basic linear fit using limma

# get suset of informatino of interest from annotation

annEPICSub <- annEPIC[match(rownames(mValues),annEPIC$Name),c(1:4,22,24,19)]

fit <- lmFit(mValues,modSv)
fit2 <- eBayes(fit)
DMPs <- topTable(fit2,coef=2,num=Inf,sort.by='none',genelist=annEPICSub)

### get lambda to check fit for inflation

# lambda <- median(qchisq(1-DMPs$P.Value,1))/qchisq(0.5,1)
# lambda

### generate z-scores for NHMMfdr

z <- qnorm(1-DMPs$P.Value)

### fit NHMMfdr

fitHMM <- fdr.nhmm(z,alttype="mixnormal",L=1)

# Z should be a matrix of covariates (try CpG island status)
# Also try 
# dist is a vector indicating distances?
#fitNHMM <- fdr.nhmm(z,Z = annEPICSub$Relation_to_Island, alttype="mixnormal",L=2)
# ix <- annEPICSub$chr == "chr1"
# fitNHMM <- fdr.nhmm(x = z[ix], dist = annEPICSub[ix, ]$pos)


LIS.adjust <- LIS.adjust(fitNHMM$LIS,fdr=0.05,adjust=TRUE)
DMPs <- cbind(DMPs,LIS.adjust$States,LIS.adjust$aLIS)
length(which(LIS.adjust$States==1))

fwrite(DMPs,file="adrcSmoking.DMPs.csv",quote=F,row.names=T)
save(c(DMPs, fit2),file="adrcSmoking.DMPs.RData")

### repeat analysis with other variables of interest