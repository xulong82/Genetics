library(rstan)

rm(list = ls())
setwd("~/GitHub/Genetics/rstan")

# ?1: does prior affect MLE results?
# a1: yes, so in building objective function of mle, stan uses the full posterior distribution, instead of the likelihood

x = rnorm(10, 0, 1)
y = 3 * x

stan.lm <- stan_model("./prior.stan")

stan.dt <- list(N = 10, x = x, y = y, prior = 3)
stan.dt <- list(N = 10, x = x, y = y, prior = 100)

fit.mle = optimizing(stan.lm, data = stan.dt)
fit.mle$par

fit.mc = sampling(stan.lm, data = stan.dt, chain = 1, iter = 1200, warmup = 200)
print(fit.mc)

samples = as.data.frame(fit.mc)

var(samples$beta)
sqrt(var(samples$beta))
sd(samples$beta)
sd(samples$beta) / sqrt(nrow(samples))
