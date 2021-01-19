# library(dplyr)  # convert factor to numeric
library('tidyr')
library('rstan')
setwd('C:/GitHub/approx_gibbs_for_HBM')

# dependent source code
source('data_preprocessing.R')
source('gibbs_all_vectorized.R')
source('diagnostics_and_plotting.R')
options(digits = 5)


## generate simulated data
generate_data = function(n_point_group, n_group, n_cov, normalize=1, add_noise=1){
  
  n_total = n_point_group*n_group
  
  # generate covariates
  x_group = runif(n_group, 1e3, 1e5)  
  x_group = matrix(rep(x_group, each = n_point_group), n_total,1)
  
  x_indiv = matrix(NA, n_total, n_cov-1)
  lb = c(10,1,1,1,50,5,10,10,100)
  ub = c(100,20,10,5,1500,50,500,1000,1000)
  for (i in 1:ncol(x_indiv)){
    for (j in 1:n_group){
      row_index = seq(((j-1)*n_point_group+1), j*n_point_group, by=1)
      x_indiv[row_index,i] = sort(runif(n_point_group,lb[i],ub[i]))
    }
  }
  
  # generate coeff
  mu = 1e-3
  sd = mu
  alpha = abs(rnorm(n_group, mu, sd))
  alpha_rep = rep(alpha, each = n_point_group)
  # coeff_indiv = abs(rnorm(n_group, mu, 1e-4))
  # coeff_indiv_rep = rep(coeff_indiv, each = n_point_group)
  coeff_indiv = matrix(abs(rnorm((n_cov-1)*n_group, mu, sd)), n_group,n_cov-1)
  coeff_indiv_rep = matrix(NA, n_total, n_cov-1)
  coeff_indiv = matrix(abs(rnorm((n_cov-1)*n_group, mu, sd)), n_group,n_cov-1)
  coeff_indiv_rep = matrix(NA, n_total, n_cov-1)
  for (j in 1:n_group){
    row_index = seq(((j-1)*n_point_group+1),j*n_point_group,by=1)
    coeff_indiv_rep[row_index,] = rep_row(coeff_indiv[j,], n_point_group)
  }
  
  # calculate y. y=normalized rate*n
  rate = exp(alpha_rep+rowSums(coeff_indiv_rep*x_indiv))
  rate_norm = rate
  for (j in 1:n_group){
    row_index = seq(((j-1)*n_point_group+1),j*n_point_group,by=1)
    lb_norm = max(abs(rnorm(1,5e-3,5e-3)), 1e-3)
    rate_norm[row_index] = min_max(rate[row_index],lb_norm)
  }
  
  y = round(rate_norm*x_group, 0)
  
  x = as.matrix(cbind(x_indiv, x_group))
  # scale x_group. Offset term for poisson regression for rate data. offset = log(n)
  x[,n_cov] = log(x[,n_cov])
  
  # normalize x
  if (normalize==1){
    for (i in 1:n_cov){
      x[,i] = min_max(x[,i],lb_norm=1e-4)
    }
  }
  
  # add noise
  if (add_noise==1){
    cv = 0.50
    mu_y_noise = y/100
    sigma_y_noise = mu_y_noise*cv
    len_y = length(y)
    y = y + (-1)^sample(1:len_y, replace=T)*round(rnorm(len_y, mu_y_noise, sigma_y_noise), 0)
    
    # for (i in 1???:n_cov){
    #   x_curr = x[,i]
    #   mu_x_noise = x_curr/20
    #   sigma_x_noise = mu_x_noise*cv
    #   len_x = length(x_curr)
    #   x[,i] = x_curr + (-1)^sample(1:len_x, replace=T)*rnorm(len_x, mu_x_noise, sigma_x_noise)
    #   x[,i] = abs(x[,i])
    # }
  }
  
  # handle values of x and y that are too small
  y[y<10] = 10
  x[x==0] = 1e-5
  
  data_df = as.data.frame(cbind(x, y))
  
  return(data_df)
}



##################
# metric_all_mat = matrix(NA, n_case, n_metric*2)
# 
# # simulated data
# n_case = 3L
# n_metric = 4L
# n_cov_vec = c(2L,2L,4L)    #seq(2L, 8L, length=n_case)
# n_group_vec = seq(10L, 30L, length=n_case)
# n_point_group_vec = seq(20L, 60L, length=n_case)

# for (i in 1:n_case){

# n_cov = n_cov_vec[i]
# n_group = n_group_vec[i]
# n_point_group = n_point_group_vec[i]
# group_attr_id = n_cov_vec[i]

n_cov = 6
n_group = 10
n_point_group = 20
group_attr_id = n_cov

syn_data = generate_data(n_point_group, n_group, n_cov, normalize=1, add_noise=0)
data_all = prepare_data(syn_data, group_attr_id)

# # sampling parameters
n_keep = 15000L
n_warmup = 7500L
n_chain = 4L

fit_model(data_all, group_attr_id)


# y = data_all[[2]]$y
# 
# group_id = data_all[[2]]$group_id


# # check data using scatter plots
# par(mfrow=c(2,2))
# for (i in 1:max(group_id)){
# 
#   print(i)
# 
#   xx = seq(1,length(group_id[group_id==i]))
#   n_col = length(syn_data[1,])
#   yy = syn_data[,n_col][group_id==i]
#   yy = sort(yy)
#   yy_log = log(yy)
# 
#   # fit model
#   fit_lm = lm(yy_log~xx)
#   model_sum = summary(fit_lm)
#   r2 = model_sum$r.squared
# 
#   # plot
#   plot(xx, yy_log, col = i, pch=19, main=paste(i))
#   R2_label = bquote(italic(R)^2 == .(format(r2, digits = 3)))
#   x_pos = min(xx) + 0.85*(max(xx)-min(xx))
#   y_pos = min(yy_log) + 0.25*(max(yy_log)-min(yy_log))
#   text(x=x_pos, y=y_pos, labels=R2_label)
#   print(r2)
# }
# 
