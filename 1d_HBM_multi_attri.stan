
data {
  int<lower=1> I;            // number of points in each group
  int<lower=1> J;            // number of groups
  int<lower=1> K;            // number of coefficients 
  int<lower=1> N;
  int<lower=1> Y[N];         // count outcome    
  real<lower=1> X[K,N];      // attribute matrix 
  int<lower=1> group_id[N]; 
}


parameters {
  real coeff[K, J];           
  real mu[K];
  real<lower=0> sigma2[K];
}



transformed parameters {
  real log_lambda[N];
  for (i in 1:N) {
    log_lambda[i] = dot_product(coeff[, group_id[i]], X[, i]);
  }
}


// log_lambda = columns_dot_product(coeff_rep, X)

model {
  target += normal_lpdf(mu|0,2);
  target += inv_gamma_lpdf(sigma2|1,1);
  
  for (k in 1:K) {
   target += normal_lpdf(coeff[k,]|mu[k], sqrt(sigma2[k]));
  }
  
  target += poisson_log_lpmf(Y|log_lambda);
}

//generated quantities {
// real <lower=1> y_hat[I_new*J_new];
//  for (ij in 1:(I_new*J_new)) {
//    yhat[ij] = ;
//  }
//}

// functions {
//   vector rep_each (vector arr, int n) {
//     int len_arr;
//     len_arr = num_elements(arr);
//     
//     int len_new_arr;
//     len_new_arr = len_arr*n;
//     
//     for (i in 1:len_new_arr) {
//       
// 
//     return arr;
//   }
// }





