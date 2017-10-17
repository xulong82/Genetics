data {
  int<lower=2> K;
  int<lower=0> N;
  int<lower=1> D;
  int<lower=1,upper=K> y[N];
  row_vector[D] x[N];
}
parameters {
  vector[D] beta;
  ordered[K-1] c;
}
model {
  vector[K] theta;
  for (n in 1:N) {
    real eta;
    eta = x[n] * beta;
    theta[1] = inv_logit(c[1] - eta);
    for (k in 2:(K-1))
      theta[k] = inv_logit(c[k] - eta) - inv_logit(c[k-1] - eta);
    theta[K] = 1 - inv_logit(c[K-1] - eta);
    y[n] ~ categorical(theta);
  }
}
generated quantities {
  vector[2] test;
  test[1] = inv_logit(2);
  test[2] = inv_logit(-2);
}
