library()
### generate z-scores for NHMMfdr

z <- qnorm(1-DMPs$P.Value)

### fit NHMMfdr

fitHMM <- fdr.nhmm(z,alttype="mixnormal",L=1)

# Z should be a matrix of covariates (try CpG island status)
# Also try 
# dist is a vector indicating distances?
#fitNHMM <- fdr.nhmm(z,Z = annEPICSub$Relation_to_Island, alttype="mixnormal",L=2)
# ix <- annEPICSub$chr == "chr1"
fitNHMM <- fdr.nhmm(x = z[ix], dist = annEPICSub[ix, ]$pos)


#LIS.adjust <- LIS.adjust(fitNHMM$LIS,fdr=0.05,adjust=TRUE)
#DMPs <- cbind(DMPs,LIS.adjust$States,LIS.adjust$aLIS)
#length(which(LIS.adjust$States==1))

#fwrite(DMPs,file="adrcSmoking.DMPs.csv",quote=F,row.names=T)
#save(c(DMPs, fit2),file="adrcSmoking.DMPs.RData")

### repeat analysis with other variables of interest