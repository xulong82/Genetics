library(rstan)

rm(list = ls())
setwd("~/GitHub/Genetics/rstan")

# bayesian p-value: what is it?

n <- 20 # Number
x <- 1:n # Pedictor
a <- 3 # Intercept
b <- 1 # Slope
eps <- rnorm(n, mean = 0, sd = 1)
y <- a + b*x + eps # Assemble data

plot(x, y)
summary(lm(y ~ x))

stan.dt <- list(n = n, x = x, y = y)

stanfit <- stan(file = "bpvalue.stan", data = stan.dt, iter = 1200, warmup = 200, chain = 1) 
print(stanfit, dig = 3)

posterior <- extract(stanfit, inc_warmup = FALSE)

mu.mean = apply(posterior$mu, 2, mean) # mean point estimates
plot(mu.mean, y)

y.new = t(sapply(1:1000, function(iter) {
  rnorm(n, mean = posterior$mu[iter, ], sd = posterior$sigma[iter])
}))

sq.res <- t((y - t(posterior$mu))^2)
sq.res.new <- (y.new - posterior$mu)^2

fit = rowSums(sq.res)
fit.new = rowSums(sq.res.new)

bpvalue = mean(as.numeric(fit.new - fit > 0))

# bpvalue in this particular case is a metrics to quantify model fittness
# whether bpvalue is useful in model comparison and/or hypothesis test is still a question mark

# to reply the reviewer:
# we will report credible interval, to keep the Bayesian spirit
# we will make it clear that the p-value we reported are posterior tail probability
# we think this p-value is practically useful for communication purporses for its frequentist flavor
# we think implementation of the most rigorous model to identify the most probable candidate target is the most meaningful
# we think this p-value performs well in ranking the best candidates
