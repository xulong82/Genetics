library(GMMAT)

setwd("/data/xwang/wes2")
load("~/Dropbox/GitHub/wes/mdata.rdt")

mdata = mdata[c("AD", "Age", "Sex")]

eigenvec <- read.table("~/Dropbox/GitHub/wes/plink/wes_pca.eigenvec")[-c(1:2)]
names(eigenvec) <- paste0("PC", 1:20)

mdata = cbind(mdata, eigenvec[1:2])

# NULL model

load("./plink/kin.rdt")
# model0 = glmmkin(AD ~ Age + Sex, data = mdata, kins = as.matrix(kin), family = binomial(link = "logit"))
# model0 = glmmkin(AD ~ Age + Sex + PC1, data = mdata, kins = as.matrix(kin), family = binomial(link = "logit"))
model0 = glmmkin(AD ~ Age + Sex + PC1 + PC2, data = mdata, kins = as.matrix(kin), family = binomial(link = "logit"))

# save(model0, file = "./glmm/model0.rdt")

# Score test

# load("./glmm/model0.rdt")
glmm.score(model0, infile = "./plink/wes", outfile = "./glmm/gmmat.score.bed.pc1.pc2.txt")

# Wald test

# load("./plink/kin.rdt")
# snps = read.table("./plink/wes.bim", stringsAsFactors = F)

# gmmat.wald = glmm.wald(AD ~ Age + Sex, data = mdata, kins = as.matrix(kin), family = binomial(link = "logit"), infile = "./plink/wes", snps = snps$V2)
# save(gmmat.wald, file = "./glmm/gmmat.wald.rdt")

