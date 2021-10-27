##### Differential analysis of ADRC microarray data
# conda environment 850k
# Writtn by Andy Madrid
# 

rm(list=ls())
getwd()
set.seed(1234)
library(limma)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(missMethyl)
library(stringr)
library(sva)
library(data.table)

load("adrcPreprocessedData.RData")
out_dir_name <- "./LOAD/"

### set variables for making models
# only groupwise model shown below

mValues <- mValues[, targets$patient_id]

bmi <- as.numeric(targets$bmi)
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

bmi <- as.numeric(targets$bmi)


### generate model and identify surrogate variables
mod <- model.matrix( ~ bmi + group + sex + age + slide + gran + nk + mono + bcell + cd8 + cd4)
mod0 <- model.matrix(~  group + sex + age + slide + gran + nk + mono + bcell + cd8 + cd4)

n.sv <- num.sv(mValues,mod,method='leek')
svobj <- sva(mValues,mod,mod0,n.sv=n.sv)
modSv <- cbind(mod, svobj$sv)

### basic linear fit using limma

# get suset of informatino of interest from annotation
# Grab the most important columns--chrom, pos, rs id, gene set, CpG island, etc.
annEPICSub <- annEPIC[match(rownames(mValues),annEPIC$Name),c(1:4,22,24,19)]

fit <- lmFit(mValues,modSv)
fit2 <- eBayes(fit)

tmp <- rep(0, ncol(modSv))
tmp[2] <- 1
contrast <- matrix(tmp, length(tmp))
contrasts.fit(fit, contrast)
fit.cont <- lmFit(mValues)

DMPs <- topTable(fit2,coef=2,num=Inf,sort.by='none',genelist=annEPICSub)

saveRDS(fit2, file = paste0(out_dir_name, "fit2.RDS"))
fwrite(DMPs, file = paste0(out_dir_name, "DMPs.csv"))
