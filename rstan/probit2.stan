data {
  int<lower=2> K;
  int<lower=0> N;
  int<lower=1> D;
  row_vector[D] x[N];
  vector[N] g; // genotype
  int<lower=1,upper=K> y[N];
}
parameters {
  real p;
  vector[D] beta;
  ordered[K-1] c;
}
model {
  vector[K] theta;
  p ~ normal(0, 1);
  beta ~ normal(0, 1);
  for (n in 1:N) {
    real eta;
    eta = p * g[n] + x[n] * beta;
    theta[1] = Phi(c[1] - eta);
    for (k in 2:(K-1))
      theta[k] = Phi(c[k] - eta) - Phi(c[k-1] - eta);
    theta[K] = 1 - Phi(c[K-1] - eta);
    y[n] ~ categorical(theta);
  }
}
generated quantities {
  vector[2] test;
  test[1] = Phi(2);
  test[2] = Phi(-2);
}
