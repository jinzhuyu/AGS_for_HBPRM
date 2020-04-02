data {
  int<lower=1> I;            // number of points in each group
  int<lower=1> J;            // number of groups
  int<lower=1> K;            // number of coefficients 
  int<lower=1> N;
  int<lower=1> Y[N];         // Count outcome    
  real<lower=1> X[J,N];            // predictor matrix 
  int<lower=1> index_group[N]; //
}

parameters {
  real coeff[K, J];           
  real<lower=0> sigma2[K];
  real mu[K];
}

transformed parameters {
  real log_lambda[N]; //or named eta
  for (i in 1:N){
    log_lambda[i]= dot_product(coeff[, index_group[i]], X[, i]);
  }
  // log_lambda = rows_dot_product(coeff[, index_group[i]], X);
}

model {
  target +=normal_lpdf(mu|0,2);         // flat enough
  target +=inv_gamma_lpdf(sigma2|1,1);

  // coefficients
  for (k in 1:K){
    target +=normal_lpdf(coeff[k, ]|mu[k],sigma2[k]);
  }
  
  for (i in 1:N){
     target += poisson_log_lpmf(Y[i]|log_lambda[i]);
  }
}

//generated quantities {
// real <lower=1> y_hat[I_new*J_new];
//  for (ij in 1:(I_new*J_new)) {
//    yhat[ij] = ;
//  }
//}



