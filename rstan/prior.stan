data {
  int<lower=1> N; 
  real x[N]; // x
  real y[N]; // response
  real prior;
}

parameters {
  real beta;  // slope
} 

model {
  beta ~ normal(prior, 1);

  for (n in 1:N)
    y[n] ~ normal(beta * x[n], 1);
}

