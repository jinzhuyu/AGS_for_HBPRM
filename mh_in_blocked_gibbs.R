
################################################################################
####  metropolis-hastings within blocked gibbs sampler for HBMs
################################################################################

library('invgamma')
library('mvtnorm')
library('coda')
library('MCMCpack')
library('rstan')
# library('profvis') # profiling the code to enhance efficiency
library('abind')


set.seed(10)


################################################################################
#### 0 data
setwd('C:/GitHub/HBKM_grouped_count_data')

# load data
data_all_0 = read.csv("C:/GitHub/HBKM_grouped_count_data/hours_peopleRegainPower_grouped.csv", header = F)

data_all = as.matrix(data_all_0)

Y_mat = data_all[2:nrow(data_all), 6:10]
y = matrix(as.numeric(unlist(Y_mat)), nrow(Y_mat), ncol(Y_mat))

time_mat_char = data_all[2:nrow(data_all), 1:5]
intensity = data_all[1, 6:10]
n_attr = 2
n_group = ncol(Y_mat)
n_point_each = nrow(Y_mat)
n_time = n_point_each

## Hyperpriors
#### vague normal
mu_0 = 0
sigma_0_square = 0.5
#### vague inverse gamma_mat
a_0 = 2
b_0 = 2

# parameters for hyperpriors
tau2_all = 1

a_all = 2
b_all = 2


Y_vec=array(as.numeric(unlist(Y_mat[train,])))
time_vec=array(as.numeric(unlist(time_mat_char[train,]))) 
intensity_vec=rep(array(as.numeric(unlist(intensity))), each=length(train))
intensity_vec1=round(intensity_vec/1000, 0)
index_group_vec=rep(1:n_group, each=length(train))
len_vec=length(Y_vec)


# individual level attributes
time_vec=array(as.numeric(unlist(time_mat_char)))   # careful here
time_mat = matrix(time_vec, n_point_each, n_group)

X_indiv_3d_arr = array(NA, dim=c(n_point_each, n_group, length(indiv_coeff_id)))
X_indiv_3d_arr[,,1] = time_mat # values of one attribute in a matrix

intensity_numeric = array(as.numeric(unlist(intensity)))
X_group_mat = intensity_numeric  # values of one attribute in a column. If there is one, then it becomes a vector
X_group_mat_ext = cbind(rep(1,n_group), X_group_mat)

X_ext_3d = array(NA, dim = c(n_point_each, n_group, n_coeff))
X_ext_3d[,,1] = array(1, dim=c(n_point_each, n_group))

for (ii in group_coeff_id[-1]) {
  id_group_temp = which(group_coeff_id == ii)
  X_ext_3d[,,ii] = rep_row(t(X_group_mat_ext[,id_group_temp]), n_point_each)
}


for (ii in indiv_coeff_id) {
  id_indiv_temp = which(indiv_coeff_id==ii)
  X_ext_3d[,,ii] = X_indiv_3d_arr[,,id_indiv_temp]
}


#####################################################################
#### 1 functions to sample from conditional posterior distributions

# function for repeating rows
rep_row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

# proposal function

prop_norm = function(coeff_curr, coeff_curr_id){
  #  input: coeff - a matrix; 
  # e.g.: alpha1, gamma1
  # alpha2, gamma2
  # alpha3, gamma3
  # coeff_curr_id - a vector.
  #  output: a matrix
  coeff_curr_new = coeff_curr
  for (ii in 1:n_group){
    coeff_curr_new[coeff_curr_id,ii] = rmvnorm(1, mean = coeff_curr[coeff_curr_id, ii],
                                               sigma=diag(length(coeff_curr_id)))
    # would it lead to more efficient proposals to set sd = mean, when mean is very small, say 0.05
  }
  return(coeff_curr_new)
}

# metropolis draw

mh_draw <- function(coeff_curr, coeff_curr_id, hyper_para_curr, targ_dens, prop = prop_norm,...) {
  # input: coeff_curr - a 2d array of samples of coeffs;
  # coeff_curr_id - id of group-level or indiv-level coeff to be updated 
  # output: 
  # generate new x_ii, evaluate with x, return new x_ii or reject it
  coeff_curr_new <- prop(coeff_curr, coeff_curr_id) # draw a new sample from current sample
  # new group or indiv coeff values, but evaluating the dens requires group and indiv coeff and hyper paras.
  log_rate <- targ_dens(coeff_curr_new, coeff_curr_id, hyper_para_curr) -
    targ_dens(coeff_curr, coeff_curr_id, hyper_para_curr)
  is_accept = log(runif(n_group)) < log_rate
  
  coeff_curr[coeff_curr_id,is_accept] <- coeff_curr_new[coeff_curr_id,is_accept]
  
  return(coeff_curr)
}


# hyper parameters mu and sigma2
# conditional posterior for hyper parameters. Conjugate.

a_sigma2_x = a_all + n_group

draw_mu_x = function(x, sigma2_x){
  # input: x - an array, e.g. vector of alpha_j; sigma2_x: a positive 'number'double'
  # output: a 'double'
  mu_mu_x = tau2_all*sum(x)/(n_group*tau2_all+sigma2_x)
  sigma2_mu_x = tau2_all*sigma2_x/(n_group*tau2_all+sigma2_x)
  return(rnorm(1, mu_mu_x, sqrt(sigma2_mu_x)))
}

draw_sigma2_x = function(x, mu_x){
  b_sigma2_x = b_all+sum((x-mu_x)^2)/2
  # note that invgamma pdf is defined over the support x > 0
  return(rinvgamma(1, a_sigma2_x, b_sigma2_x))
}

# regression coefficients

# Stirling's formula to calculate factorial of n to avoid exceeding the limit of computer arithmetic
log_fac_approx = function(x){
  return(x*log(x)-x+0.5*log(x)+0.5*log(2*pi*x))
}

log_post_orig = function(coeff_curr, coeff_curr_id, hyper_para_curr, ...){
  # prior: p(beta_j|mu_beta, sigma_beta)
  # output: a vector of conditional posterior of each indiv coeff, e.g. p(beta_j|rest)
  coeff_curr_subset = coeff_curr[coeff_curr_id,] ## change the index if there are multiple individual level paras
  if (length(coeff_curr_id)==1){
    # individual level. Updated one at a time.
    log_prior = dnorm(coeff_curr_subset, hyper_para_curr[1], hyper_para_curr[2], log = T)
  }
  else{
    # group-level coeffs updated in the same block
    log_prior = colSums(dnorm(coeff_curr_subset, hyper_para_curr[,1], hyper_para_curr[,2], log = T))
  }
  
  # log of likeli: the prod over i of p(y_ij|log_lambda_ij). log_lambda_ij = alpha_j + beta_j*x1_ij + gamma_j*x2_j
  # coeff_curr = coeff #coeff[,ii,]
  coeff_curr_rep = apply(coeff_curr, 1, function(x) rep_row(x,n_point_each))
  coeff_curr_3d = array(as.vector(coeff_curr_rep), dim = c(n_point_each,n_group,n_coeff))
  
  log_lambda = rowSums(coeff_curr_3d*X_ext_3d, dim=2) # log_lambda_ij = alpha_j*1 + beta_j*t_ij + gamma_j*I_j
  log_likeli = colSums(-log_lambda + y*log_lambda - log_fac_approx(y))
  
  # log of posteriror
  log_post = log_prior + log_likeli
  
  return(log_post)
}

################################################################################
#### 2 - run the mh-in-gibbs
################################################################################

# metropolis check (metropolis step) for each proposed sample to make sure that 
# the stationary distribution is the target distribution


## Initilization
n_coeff = 3
n_group = ncol(Y_mat)
group_coeff_id = c(1, 3) # coeff[c(1,3)]: alpha_j and gamma_j
indiv_coeff_id = seq(1, n_coeff)[-group_coeff_id]
n_group_coeff = length(group_coeff_id)
n_hyper_para = 2

n_keep = 3000
n_chain = 1
n_warm_up = max(max(500, n_keep/5), 2000)
n_sample = n_warm_up + n_keep


hyper_para = array(NA, dim = c(n_coeff, n_sample, n_hyper_para))
hyper_para[,1,] = array(1, dim = c(n_coeff, n_hyper_para))
coeff = array(NA, dim = c(n_coeff, n_sample, n_group))
coeff[,1,] = array(1, dim = c(n_coeff, n_group))

t_start_gibbs = Sys.time()

for (ii in 2:n_sample){
  # print(ii)
  # update hyperpriors
  for (jj in 1:n_hyper_para){
    if (jj==1){
      for (kk in 1:n_coeff){
        hyper_para[kk,ii,jj] = draw_mu_x(coeff[kk,ii-1,], hyper_para[kk,ii-1,-jj])
      }
    }
    else{
      for (kk in 1:n_coeff){
        hyper_para[kk,ii,jj] = draw_sigma2_x(coeff[kk,ii-1,], hyper_para[kk,ii,-jj])
      }
    }
  }
  

  coeff_temp = coeff[,ii-1,]
  
  # update indiv-level coeffs
  for (i_indiv_id in indiv_coeff_id){
    # print(i_indiv_id)
    coeff_temp[i_indiv_id,] = mh_draw(coeff_temp, i_indiv_id, hyper_para[i_indiv_id,ii,],
                                      targ_dens = log_post_orig,  prop = prop_norm)[i_indiv_id,]
  }
  coeff[,ii,] = coeff_temp
  
  # update group-level coeffs in the same block
  coeff_temp = mh_draw(coeff_curr = coeff_temp, coeff_curr_id = group_coeff_id,
                       hyper_para_curr=hyper_para[group_coeff_id,ii,],
                       targ_dens = log_post_orig, prop = prop_norm)
  coeff[,ii,] = coeff_temp
  
}

t_end_gibbs = Sys.time()


t_gibbs = t_end_gibbs-t_start_gibbs

cat('Excecution time using Gibbs sampler: ', t_gibbs)


################################################################################
##### 3 - Summarize and visualize posterior distributions 
################################################################################

coeff_name = c('alpha','beta','gamma')
for (jj in 1:length(coeff_name)){
  par(mfrow=c(ceiling(n_group/2),2))
  for (ii in 1:n_group){
    # print(paste('Posterior of', coeff_name[jj]))
    plot(coeff[jj, (n_warm_up+1):n_sample, ii], ylab = paste(coeff_name[jj], ii))
  }
}

hyper_para_name = c('mu','sigma')
for (jj in 1:length(hyper_para_name)){
  par(mfrow=c(ceiling(n_coeff/2),2))
  for (ii in 1:n_coeff){
    plot(hyper_para[ii, (n_warm_up+1):n_sample, jj], ylab = paste(hyper_para_name[jj], ii))
  }
}



################################################################################
#### 34 - Assess bias and coverage
################################################################################







