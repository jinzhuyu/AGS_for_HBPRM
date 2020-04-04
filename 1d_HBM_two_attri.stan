
  data {
    int<lower=1> I;            // number of points in each group
    int<lower=1> J;            // number of groups
    int<lower=1> N;
    int<lower=1> Y[N];         // Count outcome    
    real<lower=0> time[N];      
    real<lower=1> intensity[N];
    int<lower=1> index_group[N];
  }
  parameters {
    real beta[J];           
    real gamma[J];          
    real<lower=0> sigma2[2];
    real mu[2];
  }
  transformed parameters {
    real log_lambda[N]; //or named eta
    for (i in 1:N){
      log_lambda[i]=beta[index_group[i]]*time[i]+gamma[index_group[i]]*intensity[i];
    }
  }
  model {
    // hyperpriors
    target +=normal_lpdf(mu|0,2);         // flat enough
    target +=inv_gamma_lpdf(sigma2|1,1);
    // coefficients
    target +=normal_lpdf(beta|mu[1],sigma2[1]);
    target +=normal_lpdf(gamma|mu[2],sigma2[2]);
  
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



