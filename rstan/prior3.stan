data {
  int<lower=1> N;  // Sample
  int<lower=1> K;  // Category
  int<lower=1> D;  // Covariate
  row_vector[D] x[N];  // Covariate
  vector[N] g;  // genotype
  int<lower=1,upper=K> y[N];  // phenotype
}

parameters {
  real p; // variant effect
  vector[D] beta; // covariant effect
  simplex[K-1] cut1;  // cut points
}

transformed parameters {
  ordered[K-1] cut2;
  real cutsum;
 
  cut2[1] = 0;
  cutsum = 0;

  for (cutp in 2:(K-1)) {
    cutsum = cutsum + cut1[cutp-1];
    cut2[cutp] = 10 * cutsum;
  }
}

model {
  cut1 ~ dirichlet(rep_vector(1, K-1));
# p ~ normal(0, 1); # flat
  beta ~ normal(0, 1);

  for (n in 1:N)
    y[n] ~ ordered_logistic(p * g[n] + x[n] * beta, cut2);
}
