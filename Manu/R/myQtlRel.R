library(QTLRel)
library(parallel)

rm(list = ls())

chr <- commandArgs(TRUE) 

load("~/Dropbox/GitHub/wgs2/Manu/R/mdata.rdt")
load("~/Dropbox/GitHub/wgs2/Manu/R/kin.rdt")
load(paste0("/data/xwang/adsp3/R/chr", chr, ".rdt"))

geno <- geno[, mdata$ADSP.Sample.ID] + 1 # QTLRel use 1,2,3
kin <- as.matrix(kin[mdata$ADSP.Sample.ID, mdata$ADSP.Sample.ID])

n.core <- 15 
u.core <- round(nrow(geno) / n.core)
idx1 <- ((1:n.core) - 1) * u.core + 1 
idx2 <- c((1:(n.core - 1)) * u.core, nrow(geno))
gList <- mclapply(1:15, function(x) geno[idx1[x]:idx2[x], ], mc.preschedule = F)

# VC
# Code 1
# Y <- rep(0, nrow(mdata)) 
# Y[mdata$AD1 == "Possible"] = 0.25
# Y[mdata$AD1 == "Probable"] = 0.5
# Y[mdata$AD1 == "Definite"] = 1

# Code 2
Y <- rep(0, nrow(mdata)) 
Y[mdata$AD1 == "Possible"] = 0.33
Y[mdata$AD1 == "Probable"] = 0.66
Y[mdata$AD1 == "Definite"] = 1

X <- mdata[c("Age", "Sex")]
o <- estVC(y=Y, x=X, v=list(AA=kin, DD=NULL, HH=NULL, AD=NULL, MH=NULL, EE=diag(570)))

fit <- mclapply(gList, mc.cores = n.core, FUN = function (gdat) {
         scanOne(Y, x = X, gdat = t(gdat), numGeno = TRUE, vc = o)})

# save(fit, file = paste0("/data/xwang/adsp3/QTLRel/chr", chr, ".rdt")) # Code 1
save(fit, file = paste0("/data/xwang/adsp3/QTLRel2/chr", chr, ".rdt")) # Code 2

