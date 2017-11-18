data {
  int<lower=1> N;  // Sample
  int<lower=1> K;  // Category
  int<lower=1> D;  // Covariate
  row_vector[D] x[N];  // Covariate
  vector[N] g;  // genotype
  int<lower=1,upper=K> y[N];  // phenotype
  real prior;
}

parameters {
  real p; // variant effect
  vector[D] beta; // covariant effect
  simplex[K-1] cut1;  // cut points

  real t; // prior
  real<lower=0.0001> sigma;
}

transformed parameters {
  ordered[K-1] cut2;
  real cutsum;

  real mu;
 
  mu = sigma * t;

  cut2[1] = 0;
  cutsum = 0;

  for (cutp in 2:(K-1)) {
    cutsum = cutsum + cut1[cutp-1];
    cut2[cutp] = 10 * cutsum;
  }
}

model {
  cut1 ~ dirichlet(rep_vector(1, K-1));
  beta ~ normal(0, 1);

  t ~ normal(prior, 1);
  sigma ~ inv_gamma(2, 1);
  p ~ normal(mu, sigma);
 
  for (n in 1:N)
    y[n] ~ ordered_logistic(p * g[n] + x[n] * beta, cut2);
}

