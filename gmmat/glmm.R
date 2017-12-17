library(rstan)
library(parallel)
library(dplyr)
library(GMMAT)

rm(list = ls())
setwd("~/GitHub/Genetics/gmmat")

# setup

load("./mdata.rdt") 
load("./sampling.rdt")
load("./genotypes.rdt")
load("./kin.rdt")

mcmc = mcmc[order(mcmc$P), ]
mcmc = rownames(mcmc[1:100, ]) # top 100 variants

all(mcmc %in% rownames(geno))
all(mdata$ADSP.Sample.ID %in% colnames(geno))
geno = geno[mcmc, mdata$ADSP.Sample.ID]

kin <- kin[mdata$ADSP.Sample.ID, mdata$ADSP.Sample.ID]

# bayes-glmm estimation of binary phenotypes

cov <- mdata[c("Age", "Sex")]
dat <- list(N = 570, K = 4, D = 2, cov = cov)
dat <- within(dat, { Ad = as.numeric(mdata$AD1); L = t(chol(kin)) })
dat$Ad <- as.numeric(dat$Ad > 2)

glmm <- stan_model("./glmm.stan")
dat0 <- within(dat, {g = rep(0, 570)})

fit0 = optimizing(glmm, data = dat0)
init = list(a = fit0$par["a"], sigma = fit0$par["sigma"], beta = fit0$par[paste0("beta[", 1:2, "]")])
init$z = fit0$par[paste0("z[", 1:570, "]")]
init$p = 0

mymle <- lapply(1:nrow(geno), function(g1) {
  dat$g = geno[g1, ]
  mle1 = optimizing(glmm, data = dat, hessian = TRUE, algorithm = "LBFGS", init = init)
  mle1.se = tryCatch(sqrt(diag(solve(-mle1$hessian)))["p"], error=function(e) NULL)
  data.frame(effect = mle1$par["p"], effect.se = mle1.se)
}); names(mymle) = rownames(geno)

mymle = as.data.frame(do.call(rbind, mymle))
mymle$Pval <- pnorm(abs(mymle$effect), sd = mymle$effect.se, lower.tail = F) * 2

save(mymle, file = "glmm.rdt")

# GMMAT

gmmat.dt = mdata
gmmat.dt$Ad = as.numeric(mdata$AD1 %in% c("Probable", "Definite"))
gmmat.fit = GMMAT::glmmkin(Ad ~ Age + Sex, data = gmmat.dt, kins = as.matrix(kin), family = binomial(link = "logit"))

# score test
write.table(geno, file = "./gmmat.infile", sep = "\t")
glmm.score(gmmat.fit, infile = "./gmmat.infile", outfile = "./gmmat.score.txt", infile.nrow.skip = 1, infile.ncol.skip = 1)

# summarize

load("./glmm.rdt")
gmmat.score = read.csv("./gmmat.score.txt", sep = "\t", stringsAsFactors = F)

all(rownames(mymle) == gmmat.score$SNP)

plot(-log10(mymle$Pval), -log10(gmmat.score$PVAL), 
     xlab = "-log10(Pval) :: Bayes-GLMM", ylab = "-log10(Pval) :: GMMAT")
