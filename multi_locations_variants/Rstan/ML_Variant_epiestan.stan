data {
  int<lower=0> nt;
  int<lower=0> tw;
  int<lower=0> n_location;
  int<lower=1> n_variant;
  int I[n_location*n_variant, nt, tw];
  matrix<lower=0>[nt, tw] O_I[n_location*n_variant];
  row_vector[tw] U;
  matrix<lower=0>[nt, n_location] prior_shape;
  matrix<lower=0>[nt, n_location] prior_rate;
  vector<lower=0>[2] prior_beta;
}

parameters {
  matrix<lower=0>[nt, n_location] Rt;
  row_vector<lower=0>[n_variant-1] beta;
}

transformed parameters {
  row_vector<lower=0>[n_variant] Beta;
  matrix[nt, n_variant] Rl;
  matrix[nt, tw] R;
  matrix<lower=0>[nt, tw] lambda[n_location*n_variant];
  
  Beta = append_col(1, beta);
  for (i in 1:n_location){
    Rl = Rt[,i]*Beta;
    for (j in 1:n_variant){
      R =  Rl[,j] *U;
      lambda[(i-1)*n_location + j,,] = R .* O_I[(i-1)*n_location+j,,];
    }
  }
}

model{
  for(i in 1:(nt)){
    for(j in 1:n_location){
      Rt[i,j] ~ gamma(prior_shape[i,j], prior_rate[i,j]);
    }
  }
  beta ~ gamma(prior_beta[1],prior_beta[2]);
  
  for(j in 1:(n_variant*n_location)){
    for (i in 1:nt){
      // I[j,i,] ~ poisson(lambda[j,i,]);
      target += 1;
    }
  }
      
}
