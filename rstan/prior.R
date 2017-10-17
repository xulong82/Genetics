library(rstan)

rm(list = ls())
setwd("~/GitHub/Genetics/rstan")

# ?1: does prior affect MLE results?
# a1: yes, so in building objective function of mle, stan uses the full posterior distribution, instead of the likelihood

x = rnorm(100, 0, 1)
y = 3 * x + rnorm(100, 0, 1)
plot(x, y)

stan.lm <- stan_model("./prior.stan")

stan.dt <- list(N = 100, x = x, y = y, pmu = 3, psigma = .1)
stan.dt <- list(N = 100, x = x, y = y, pmu = 100, psigma = .1)
stan.dt <- list(N = 100, x = x, y = y, pmu = 3, psigma = 100)

fit.mle = optimizing(stan.lm, data = stan.dt)
fit.mle$par

fit.mc = sampling(stan.lm, data = stan.dt, chain = 1, iter = 1200, warmup = 200)
print(fit.mc)
hist(extract(fit.mc)$beta)

samples = as.data.frame(fit.mc)

var(samples$beta)
sqrt(var(samples$beta))
sd(samples$beta)
sd(samples$beta) / sqrt(nrow(samples))

# Compare flat prior vs normal prior

load("~/GitHub/wgs2/Manu/R/mdata.rdt")
load("~/GitHub/wgs2/Manu/R/sampling.rdt")
load("~/GitHub/wgs2/Manu/R/genotypes.rdt")

mcmc = mcmc[order(mcmc$P), ]
mcmc = rownames(mcmc[1:100, ]) # top 100 variants

all(mcmc %in% rownames(geno))
all(mdata$ADSP.Sample.ID %in% colnames(geno))
geno = geno[mcmc, mdata$ADSP.Sample.ID]

prior.norm = stan_model("~/GitHub/Genetics/rstan/prior2.stan") # N(0, 1) as prior of variant effect
prior.flat = stan_model("~/GitHub/Genetics/rstan/prior3.stan") # flat as prior of variant effect
data <- list(N = 570, K = 4, D = 2, x = mdata[c("Age", "Sex")], y = as.numeric(mdata$AD1), g = geno[1, ])

fit.norm = sampling(prior.norm, data = data, chain = 1, iter = 1200, warmup = 200)
fit.flat = sampling(prior.flat, data = data, chain = 1, iter = 1200, warmup = 200)

print(fit.norm)
print(fit.flat)

fit.norm = optimizing(prior.norm, data = data)
fit.flat = optimizing(prior.flat, data = data)

tops = lapply(1:nrow(geno), function(g1) {
  data$g = geno[g1, ]
  fit.norm = optimizing(prior.norm, data = data, hessian = TRUE, algorithm = "LBFGS")
  fit.norm.se = tryCatch(sqrt(diag(solve(-fit.norm$hessian)))["p"], error=function(e) NULL)
  fit.flat = optimizing(prior.flat, data = data, hessian = TRUE, algorithm = "LBFGS")
  fit.flat.se = tryCatch(sqrt(diag(solve(-fit.flat$hessian)))["p"], error=function(e) NULL)
  data.frame(norm = fit.norm$par["p"], norm.se = fit.norm.se, flat = fit.flat$par["p"], flat.se = fit.flat.se)
})

tops = as.data.frame(do.call(rbind, tops))

tops$P.norm <- pnorm(abs(tops$norm), sd = tops$norm.se, lower.tail = F) * 2
tops$P.flat <- pnorm(abs(tops$flat), sd = tops$flat.se, lower.tail = F) * 2

tops$P.norm <- -log10(tops$P.norm)
tops$P.flat <- -log10(tops$P.flat)

plot(tops$norm, tops$flat, xlab = "Standard Normal Prior", ylab = "Flat Prior", main = "Effect Size")
abline(a = 0, b = 1, col = "red")

plot(tops$P.norm, tops$P.flat, xlab = "Standard Normal Prior", ylab = "Flat Prior", main = "Tail Probability")
abline(a = 0, b = 1, col = "red")

# Standard normal distribution as priors does shrinkage large effect sizes, and hereby shrinkage P.values. 
# Standard normal distribution comes with the prior assumption that variant of large effect sizes are less likely.
# By shrinkage large effect size and P.value, standard normal as priors are effectively decreasing false positive rates. 
# In GWAS setup, false positive is more of a concern than false negative. 
