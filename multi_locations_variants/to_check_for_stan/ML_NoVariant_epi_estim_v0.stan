data {
  int<lower=0> nt;
  int<lower=0> tw;
  int<lower=0> n_location;
  // int<lower=1> n_variant;
  int I[n_location, nt, tw]; //int I[n_location*n_variant, nt, tw];
  matrix<lower=0>[nt, tw] O_I[n_location]; //matrix<lower=0>[nt, tw] O_I[n_location*n_variant];
  row_vector[tw] U;
  real<lower=0> prior_shape;
  real<lower=0> prior_rate;
  // real<lower=0> prior_beta;
}

parameters {
  matrix<lower=0>[nt, n_location] Rt;
  // row_vector<lower=0>[n_variant-1] beta;
}

transformed parameters {
  // row_vector<lower=0>[n_variant] Beta;
  //matrix[nt, n_variant] Rl;
  matrix[nt, tw] R;
  matrix<lower=0>[nt, tw] lambda[n_location]; //matrix<lower=0>[nt, tw] lambda[n_location*n_variant];
  
  // Beta = append_col(1, beta);
  for (i in 1:n_location){
    // Rl = Rt[,i]*Beta;
    // for (j in 1:n_variant){
      R =  Rt[,i] *U; //R =  Rl[,j] *U;
      lambda[i,,] = R .* O_I[i,,];
      // lambda[(i-1)*n_location + j,,] = R .* O_I[(i-1)*n_location+j,,];
    //}
  }
}

model{
  to_vector(Rt) ~ gamma(prior_shape, prior_rate);
  for(j in 1:(n_location)){ //for(j in 1:(n_variant*n_location)){
    for (i in 1:nt){
      I[j,i,] ~ poisson(lambda[j,i,]);
    }
  }
      
}
