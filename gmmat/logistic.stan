data {
  int<lower=1> N;  // Sample
  int<lower=1> D;  // Covariate
  row_vector[D] cov[N];  // Covariate
  vector[N] g;  // genotype
  int<lower=0,upper=1> Ad[N];  // phenotype
}

parameters {
  real a;  // intercept
  real p;  // variant effect
  vector[D] beta;  // covariate effect
} 

model {
  p ~ normal(0, 1);
  beta ~ normal(0, 1);

  for (n in 1:N)
    Ad[n] ~ bernoulli_logit(a + p * g[n] + cov[n] * beta);
}

generated quantities {
  vector[N] lp;

  for (n in 1:N)
    lp[n] <- bernoulli_logit_log(Ad[n], a + p * g[n] + cov[n] * beta);
}

