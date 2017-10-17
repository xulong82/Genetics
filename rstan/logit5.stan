data {
  int<lower=2> K;
  int<lower=0> N;
  int<lower=1> D;
  int<lower=1,upper=K> y[N];
  row_vector[D] x[N];
}
parameters {
  vector[D] beta;
# ordered[K-1] c;
  simplex[K-1] cut1;  // cut points
}
transformed parameters {
  ordered[K-1] cut2;
  real cutsum;
 
  cut2[1] = 0;
  cutsum = 0;

  for (cutp in 2:(K-1)) {
    cutsum = cutsum + cut1[cutp-1];
    cut2[cutp] = 100 * cutsum;
  }
}
model {
  cut1 ~ dirichlet(rep_vector(1, K-1));
  for (n in 1:N)
    y[n] ~ ordered_logistic(x[n] * beta, cut2);
}

