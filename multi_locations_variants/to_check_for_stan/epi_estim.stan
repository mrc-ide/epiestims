data {
  int<lower=0> nR;
  int<lower=0> tw;
  int I[nR, tw];
  matrix<lower=0>[nR, tw] O_I;
  row_vector[tw] U;
  real<lower=0> prior_alpha;
  real<lower=0> prior_beta;
}

parameters {
  vector<lower=0>[nR] Rt;
}

transformed parameters {
  matrix[nR, tw] R;
  matrix[nR, tw] lambda;
  
  R =  Rt *U;
  lambda = R .* O_I;
}

model{
  Rt ~ gamma(prior_alpha, prior_beta);
  for (i in 1:nR)
    I[i] ~ poisson(lambda[i]);
}
