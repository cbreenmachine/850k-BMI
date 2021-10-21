# Written by Andy Madrid

### Load packages for preprocessing and get annotation

rm(list=ls())
getwd()
set.seed(1234)
library(limma)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(missMethyl)
library(stringr)
library(data.table)

annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(annEPIC)

### Import raw idat file into R environment

targets <- read.csv("../Data/Samplesheets/master-WB.csv",header=T)
rel.names <- paste("../Data/Raw/", targets$patient_folder, targets$patient_id, sep = "/")

# Pull in red-green channel set
RGSet <- read.metharray(rel.names)
RGSet

### Basic QC based on detection levels of probes
detP <- detectionP(RGSet)

# non-processed	data for more QC of raw data
# preprocessRaw just "maps" Reg/Green channel into methylation 
# without any normalization
mSetRaw	<- preprocessRaw(RGSet)
qcRaw <- getQC(mSetRaw)
plotQC(qcRaw)

# save QC report

#qcReport(RGSet,sampNames=targets$adrcnum,sampGroups=targets$COHORT,pdf="qcReport_minfi_ADRC_Project.pdf")

# filter poor samples with high mean detection P-values (>0.05)

keep <- colMeans(detP) < 0.05
RGSet <- RGSet[,keep]
targets <- targets[keep,]
detP <- detP[,keep]


### Normalize data

# background and control normalize with "GenomeStudio" standards
mSetIllumina <- preprocessIllumina(RGSet,bg.correct=TRUE,normalize='controls')

# Memory management
rm(mSetRaw)
rm(qcRaw)

# within array normalization
mSetSWAN <- preprocessSWAN(RGSet,mSet=mSetIllumina,verbose=TRUE)
rm(mSetIllumina)
gc()

# check QC of normalized data

#qcSWAN <- getQC(mSetSWAN)
#plotQC(qcSWAN)

# check predicted sex of each sample from normalized data

mSetSWAN <- mapToGenome(mSetSWAN)
#pSex <- getSex(mSetSWAN)

# get estimated cell counts for blood tissue
# I've noticed some errors get thrown depending on which version of minfi is being used...

#cellCounts <- estimateCellCounts(RGSet,compositeCellType="Blood", referencePlatform = "IlluminaHumanMethylationEPI")

######################################

### Filter probes for low quality, sex chromosomes, cross-reactive, CH

# filter if at least one sample has detP > 0.01

detP <- detP[match(featureNames(mSetSWAN),rownames(detP)),]
keep <- rowSums(detP < 0.01) == ncol(mSetSWAN)
mSetFiltered <- mSetSWAN[keep,]

rm(detP)
gc()

# filter probes on X and Y chromosomes

keep <- !(featureNames(mSetFiltered) %in% annEPIC$Name[annEPIC$chr %in% c("chrX","chrY")])
mSetFiltered <- mSetFiltered[keep,]

# filter probes with SNPs and at CH sites

mSetFiltered <- dropLociWithSnps(mSetFiltered)
mSetFiltered <- dropMethylationLoci(mSetFiltered)
gc()


# filter probes of known cross-reactive sites

xRtvProbes <- read.csv("pub-2016McCartney-CrossHybridCpG.csv", header=F, stringsAsFactors = F)
keep <- !(featureNames(mSetFiltered) %in% xRtvProbes)
mSetFiltered <- mSetFiltered[keep,]

rm(mSetSWAN)
rm(RGSet)
gc()


### Get Beta (bValues) and logit M-values (mValues)

mValues <- getM(mSetFiltered) %>% as.data.frame()
bValues <- getBeta(mSetFiltered) %>% as.data.frame()

### Write files of mValues, bValues, estimated cell counts

fwrite(bValues,file="adrcPreprocessed.bValues.csv",quote=F,row.names=T)
fwrite(mValues,file="adrcPreprocessed.mValues.csv",quote=F,row.names=T)
#write.table(cellCounts,file="adrcCellCounts.csv",quote=F,row.names=T.sep=',')
save(bValues,mValues,annEPIC,targets,file="adrcPreprocessedData.RData")
