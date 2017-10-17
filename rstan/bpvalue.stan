data{
  int n;
  real x[n];
  real y[n];
}

parameters{
  real alpha;
  real beta;
  real<lower=0, upper=100> sigma;
}

transformed parameters{
  real mu[n];
  for(i in 1:n){
    mu[i] = alpha + beta*x[i];
  }
}

model{
  alpha ~ normal(0, 100);
  beta ~ normal(0, 100);
  sigma ~ uniform(0, 100);

  for(i in 1:n){
    y[i] ~ normal(mu[i], sigma);
  }
}

