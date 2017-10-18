# Ordered categorical models involve cut point parameters, as a vector of ever-growing real numbers 
# Stan implemented a function ordered() to generate it in its new update
# However, Bayes-GLMM took a different strategy to generate it in its original implementation
# This script compare performances of the two strategies

# Bayers-GLMM strategy:
# theta0 ~ dirichlet(1)
# theta = 10 * cumsum(theta0)

# dirichlet distribution
library(MCMCpack)
dirich = rdirichlet(1e3, alpha = rep(1, 3))
dirich = rdirichlet(1e3, alpha = rep(2, 3))
dirich = rdirichlet(1e3, alpha = c(1, 2, 3))
table(rowSums(dirich)) # each sample is a simplex vector
plot(density(c(dirich)))
summary(c(dirich))

library(rstan)
setwd("~/GitHub/Genetics/rstan/")

# Bayes-GLMM method accounts for the identificability issue by fixing the first cutpoint as 0, while the ordered() function does not
# identificability between intercept and cut points  
# we do not have the option to ignore fitting intercept in Stan

load("mdata.rdt") # ADSP
dat <- list(N = 570, K = 4, D = 2, x = mdata[c("Age", "Sex")], y = as.numeric(mdata$AD1))

# where the scale 10 came from?

logit = stan_model("logit.stan") # Stan manual
logit.mle = optimizing(logit, data = dat)

# a straight ordered() results suggests cut points ranges 0-10

logit3 = stan_model("logit3.stan") # modified version: Bayes-GLMM method for cutoff: 10*cumsum
logit4 = stan_model("logit4.stan") # modified version: Bayes-GLMM method for cutoff: 1*cumsum
logit5 = stan_model("logit5.stan") # modified version: Bayes-GLMM method for cutoff: 100*cumsum

logit3.mle = optimizing(logit3, data = dat) # 10 * cumsum
logit4.mle = optimizing(logit4, data = dat) # 1 * cumsum
logit5.mle = optimizing(logit5, data = dat) # 100 * cumsum

# logit4 fit was much worse, when cut values were restricted between [0-1] 
# logit3 and 5 both fit equally well, with similar cut point estimations 

# Are the 10 * cumsum() data-depedent?

# predictor variables?
dat2 = dat
dat2$x$Age = 2 * dat2$x$Age 
dat2$x$Sex = -20 + dat2$x$Sex

# does not depend on predictor range

# response variable?
# more than 10 levels?

dat2 = dat
dat2$y = sample(1:20, size = dat2$N, replace = T)
dat2$K = 20

logit.mle.dt1 = optimizing(logit, data = dat)
logit.mle.dt2 = optimizing(logit, data = dat2)

logit3.mle.dt2 = optimizing(logit3, data = dat2)
logit4.mle.dt2 = optimizing(logit4, data = dat2)
logit5.mle.dt2 = optimizing(logit5, data = dat2)

logit.mle.dt1
logit.mle.dt2

logit3.mle.dt2
logit4.mle.dt2
logit5.mle.dt2

# again, logit4 fit was much worse, while logit3 and logit5 fits are similar
# furhter, cut points estimation values ranges are still 0-10
# so the scale factor 10 does not depend on response variable neither
