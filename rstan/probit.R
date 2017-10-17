# probit and logistic links

# Pr(Y=1|X) = exp(x)/(1+exp(x)) # cumulative logistic distribution 
# Pr(Y=1|X) = Pho(x) # cumulative normal distribution

# logit is more popular in health sciences since coefficients can be interpreted in terms of odds ratios
# x = ln(p/(1-p))

x = seq(-10, 10, 0.1)

cdf.logi = exp(x) / (1+exp(x))
cdf.prob = pnorm(x, mean = 0, sd = 1)

plot(x, cdf.logi, col = "blue", type = "l")
lines(x, cdf.prob, col = "red")

dens.logi = dlogis(x, location = 0, scale = 1)
dens.prob = dnorm(x, mean = 0, sd = 1)

plot(x, dens.logi, col = "blue", type = "l", xlim = c(-6, 6))
# logistic distribution is symmetric, and is pretty much a normal distribution with some standard deviation
# so it is straight forward to approximate a logistic by a normal, or the reverse
mean(dens.logi)
sd(dens.logi)

plot(x, dens.logi, col = "blue", type = "l", ylim = c(0, 0.5))
lines(x, dens.prob, col = "red")

# more comparisons on fit reasults, likelihood and parameters et al.

