data {
  int<lower=1> N; 
  real x[N]; // x
  real y[N]; // response
  real pmu;
  real psigma;
}

parameters {
  real beta;  // slope
  real sigma; // variance
} 

model {
  beta ~ normal(pmu, psigma);

  for (n in 1:N)
    y[n] ~ normal(beta * x[n], sigma);
}

