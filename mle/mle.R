library(rstan)
library(dplyr)

# ? consistency of inferences by MLE and MC

rm(list = ls())
setwd("~/gitHub/Genetics/mle")

load("./mdata.rdt") 
load("./sampling.rdt")
load("./genotypes.rdt")

mcmc = mcmc[order(mcmc$P), ]
mcmc = rownames(mcmc[1:100, ]) # top 100 variants

all(mcmc %in% rownames(geno))
all(mdata$ADSP.Sample.ID %in% colnames(geno))
geno = geno[mcmc, mdata$ADSP.Sample.ID]

glm <- stan_model("./glm.stan")
data <- list(N = 570, K = 4, D = 2, x = mdata[c("Age", "Sex")], y = as.numeric(mdata$AD1), g = geno[1, ])

mymle <- lapply(1:nrow(geno), function(g1) {
  data$g = geno[g1, ]
  mle1 = optimizing(glm, data = data, hessian = TRUE, algorithm = "LBFGS")
  mle1.se = tryCatch(sqrt(diag(solve(-mle1$hessian)))["p"], error=function(e) NULL)
  data.frame(effect = mle1$par["p"], effect.se = mle1.se)
}); names(mymle) = rownames(geno)

mymle = as.data.frame(do.call(rbind, mymle))
mymle$Pval <- pnorm(abs(mymle$effect), sd = mymle$effect.se, lower.tail = F) * 2

mymcmc <- lapply(1:nrow(geno), function(g1) {
  data$g = geno[g1, ]
  mcmc1 = sampling(glm, data = data, chain = 1, iter = 1200, warmup = 200) 
  mcmc1 = summary(mcmc1, pars = c("p", "beta"))$summary
  mcmc1 = mcmc1["p", c("mean", "sd")]
}); names(mymcmc) = rownames(geno)

mymcmc = as.data.frame(do.call(rbind, mymcmc))
mymcmc$Pval <- pnorm(abs(mymcmc$mean), sd = mymcmc$sd, lower.tail = F) * 2

plot(mymle$effect, mymcmc$mean)
plot(mymle$effect.se, mymcmc$sd)
plot(-log10(mymle$Pval), -log10(mymcmc$Pval))
abline(a = 0, b = 1)

# MLE/Hessian method agrees with MCMC 
