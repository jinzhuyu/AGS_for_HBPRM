
################################################################################
####  metropolis-hastings within gibbs sampler for HBMs
################################################################################

library('invgamma')
library('mvtnorm')
library('coda')
library('MCMCpack')
library('rstan')
library('profvis') # profiling the code to enhance efficiency
library('abind')

set.seed(10)


################################################################################
#### 0 data
setwd('C:/GitHub/approx_gibbs_for_HBM')

# load data
data_all_0 = read.csv("C:/GitHub/approx_gibbs_for_HBM/power_recovery_data.csv", header = F)

data_all = as.matrix(data_all_0)

Y_mat = data_all[2:nrow(data_all), 6:10]
y = matrix(as.numeric(unlist(Y_mat)), nrow(Y_mat), ncol(Y_mat))

time_mat_char = data_all[2:nrow(data_all), 1:5]
intensity = data_all[1, 6:10]
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


set.seed(1)

train=1:(nrow(Y_mat)-2)

# I_new=2,J_new=n_group,
# x1_new=time_mat[-train,])#, x2_new=intensity)

Y_vec = array(as.numeric(unlist(Y_mat[train,])))
time_vec = array(as.numeric(unlist(time_mat_char[train,])))    # careful here
intensity_vec = rep(array(as.numeric(unlist(intensity))), each=length(train))
intensity_vec1 = round(intensity_vec/1000, 0)
index_group_vec = rep(1:n_group, each=length(train))
len_vec = length(Y_vec)


# individual level attributes
time_vec = array(as.numeric(unlist(time_mat_char)))   # careful here
time_mat = matrix(time_vec, n_point_each, n_group)

X_indiv_3d_arr = array(NA, dim=c(n_point_each, n_group, length(indiv_coeff_id)))
X_indiv_3d_arr[,,1] = time_mat # values of one attribute in a matrix

intensity_numeric = array(as.numeric(unlist(intensity)))
X_group_mat = intensity_numeric  # values of one attribute in a column. If there is one, then it becomes a vector


X_3d = array(NA, dim = c(n_point_each, n_group, n_coeff))

if (length(group_coeff_id)==1){
  X_3d[,,group_coeff_id] = rep_row(X_group_mat, n_point_each)
}else{
  for (ii in group_coeff_id) {
    id_group_temp = which(group_coeff_id==ii)
    X_3d[,,ii] = rep_row(t(X_group_mat[,id_group_temp]), n_point_each)
  }
}

for (ii in indiv_coeff_id) {
  id_indiv_temp = which(indiv_coeff_id==ii)
  X_3d[,,ii] = X_indiv_3d_arr[,,id_indiv_temp]
}
# test
# X_group = matrix(1:(n_group*length(group_coeff_id)),c(n_group,length(group_coeff_id)))
# X_indiv_1 = X_group*rnorm(1, 1, 1)
# X_indiv_2 = X_group*rnorm(1, 1, 1)
# X_indiv_3d_arr = abind(X_indiv_1,X_indiv_2, along=3)
# coeff_curr = matrix(rnorm(n_group*n_coeff, 1, 1), c(n_coeff,n_group))


#####################################################################
#### 1 functions to sample from conditional posterior distributions

# function for repeating rows
rep_row = function(x,n){
  matrix(rep(x,each=n),nrow=n)
}


# proposal function
prop_norm = function(coeff_curr, coeff_curr_id, hyper_param_curr, ...){
  # input: coeff_curr - a vector; 
  # e.g.: beta1, beta2
  # coeff_curr_id - an integer.
  # hyper_param_curr is redundant here, but is needed in the proposal using approximate likelihood
  # output: a vector
  coeff_curr_1d = coeff_curr[coeff_curr_id, ]
  coeff_curr_new_1d = coeff_curr_1d
  sd = c(0.02/32, 3e-6/13)
  
  for (jj in 1:n_group) {
    coeff_curr_new_1d[jj] = rnorm(1, mean=coeff_curr_1d[jj], sd=sd[coeff_curr_id])
  }
  
  return(coeff_curr_new_1d)
}


# the prop_normal does not depende on other coeff and hyper params
prop_approx_likeli = function(coeff_curr, coeff_curr_id, hyper_param_curr, ...){
  # hyper param curr: a vector
  
  denominator = hyper_param_curr[2]^2*colSums(y*X_3d[,,coeff_curr_id]^2) + 1
  
  # coeff_curr = coeff #coeff[,ii,]
  coeff_curr_rep = apply(coeff_curr, 1, function(x) rep_row(x, n_point_each))
  coeff_curr_3d = array(as.vector(coeff_curr_rep), dim = c(n_point_each, n_group, n_coeff))
  
  prod_part_temp = rowSums(coeff_curr_3d*X_3d, dim=2)-
    rep_row(coeff_curr[coeff_curr_id, ], n_point_each)*X_3d[coeff_curr_id,,coeff_curr_id]
  
  
  numerator = hyper_param_curr[1] + 
    hyper_param_curr[2]^2*colSums(y*X_3d[,,coeff_curr_id]*
      (log(y)-prod_part_temp))
  
  
  mu_hat = numerator/denominator
  
  simga_hat = sqrt(hyper_param_curr[2]^2/denominator)
  
  coeff_curr_new_1d = coeff_curr[coeff_curr_id, ]
  for (jj in 1:n_group) {
  coeff_curr_new_1d[jj] = rnorm(1, mean=mu_hat[jj], sd=simga_hat[jj])
  }
  
  return(coeff_curr_new_1d)
}


# metropolis draw

mh_draw = function(coeff_curr, coeff_curr_id, hyper_param_curr, targ_dens_mh, prop_mh) {
  # input: coeff_curr - a 2d array of samples of coeffs;
  # coeff_curr_id - id of group-level or indiv-level coeff to be updated 
  # output: 
  # generate new x_ii, evaluate with x, return new x_ii or reject it
  coeff_curr = as.matrix(coeff_curr)  # keep the dimension
  coeff_curr_new = coeff_curr
  coeff_curr_new[coeff_curr_id, ] = prop_mh(coeff_curr, coeff_curr_id, hyper_param_curr) # draw a new sample from current sample
  
  # new coeff values, but evaluating the dens requires the current values of other coeffs and hyper params.
  log_rate = targ_dens_mh(coeff_curr_new, coeff_curr_id, hyper_param_curr) - 
    targ_dens_mh(coeff_curr, coeff_curr_id, hyper_param_curr)
  is_accept = log(runif(n_group)) < log_rate
  
  coeff_curr[coeff_curr_id, is_accept] = coeff_curr_new[coeff_curr_id, is_accept]
  
  accept_count[coeff_curr_id, is_accept] <<- accept_count[coeff_curr_id, is_accept] + 1
  
  return(coeff_curr)
}


gibbs_with_mh = function(coeff_curr, hyper_param_curr, targ_dens, prop_fun){
  
  
  # update hyperpriors
  for (i_coeff in 1:n_coeff) {
    hyper_param_curr[i_coeff, 1] = draw_mu_x(
      coeff_curr[i_coeff, ], hyper_param_curr[i_coeff, 2])
    
    hyper_param_curr[i_coeff, 2] = draw_sigma2_x(
      coeff_curr[i_coeff, ], hyper_param_curr[i_coeff, 1])
  }
  
  
  for (i_coeff in 1:n_coeff){
    
    coeff_curr = mh_draw(coeff_curr, i_coeff, hyper_param_curr[i_coeff,],
                         targ_dens_mh = targ_dens,  prop_mh = prop_fun)
  }
  
  
  return(list(coeff_curr, hyper_param_curr))  # can R return a list of matrix?
}



gibbs_no_mh = function(coeff_curr, hyper_param_curr, approx_fun){
  
  
  # update hyperpriors
  for (i_coeff in 1:n_coeff) {
    hyper_param_curr[i_coeff, 1] = draw_mu_x(
      coeff_curr[i_coeff, ], hyper_param_curr[i_coeff, 2])
    
    hyper_param_curr[i_coeff, 2] = draw_sigma2_x(
      coeff_curr[i_coeff, ], hyper_param_curr[i_coeff, 1])
  }
  
  
  for (i_coeff in 1:n_coeff) {
    
    coeff_curr[i_coeff, ] = approx_fun(coeff_curr, i_coeff, hyper_param_curr[i_coeff,])
  }
  
  
  return(list(coeff_curr, hyper_param_curr))  # can R return a list of matrix?
}


# prop=prop_norm
# coeff_curr = coeff_temp
# coeff_curr_id = i_coeff
# hyper_param_curr = hyper_param[i_coeff,i_samp,]
# targ_dens = log_post_orig
# n_group number of log rate


# hyper parameters mu and sigma2
# conditional posterior for hyper parameters. Conjugate.

a_sigma2_x = a_all + n_group

draw_mu_x = function(x, sigma2_x) {
  # input: x - an array, e.g. vector of alpha_j; sigma2_x: a positive 'number'double'
  # output: a 'double'
  mu_mu_x = tau2_all*sum(x)/(n_group*tau2_all+sigma2_x)
  sigma2_mu_x = tau2_all*sigma2_x/(n_group*tau2_all+sigma2_x)
  return(rnorm(1, mu_mu_x, sqrt(sigma2_mu_x)))
}

draw_sigma2_x = function(x, mu_x) {
  b_sigma2_x = b_all+sum((x-mu_x)^2)/2
  # note that invgamma pdf is defined over the support x > 0
  return(rinvgamma(1, a_sigma2_x, b_sigma2_x))
}

# regression coefficients


# conditional posterior for invididual-level parameters, beta_Approximate normal
# beta_approx = function(alpha, gamma, mu_beta, sigma2_beta){
#   mu_beta_hat = (mu_beta+sigma2_beta*colSums(y_mat_gibbs*time_gibbs_mat*
#                                                (log(y_mat_gibbs)-rep_row(alpha,n_row_train)-rep_row(gamma*intensity_gibbs,n_row_train))))/
#     (sigma2_beta*sum_y_time2+1)
#   sigma2_beta_hat = 1/(sum_y_time2+1/sigma2_beta)
#   return(rnorm(n_group, mu_beta_hat, sqrt(sigma2_beta_hat)))
# }

# combine the attributes

# conditional joint posterior for alpha and gamma. 
# 1. Origninal complex distribution. Target distribution. mvnorm will be adopted as proposal dist.
# 2. Approximate the likelihood only with a normal distribution. Used as the proposal dist.
# 3. Approximate product of normal variables with a normal distribution. Used as the proposal dist.

log_post_orig = function(coeff_curr, coeff_curr_id, hyper_param_curr) {
  # prior: p(beta_j|mu_beta, sigma_beta)
  # output: a vector of conditional posterior of each indiv coeff, e.g. p(beta_j|rest)
  
  # log of prior
  coeff_curr_subset = coeff_curr[coeff_curr_id,]
  log_prior = dnorm(coeff_curr_subset, hyper_param_curr[1], hyper_param_curr[2], log=T)
  
  # log of likeli: the prod over i of p(y_ij|log_lambda_ij). log_lambda_ij = alpha_j + beta_j*x1_ij + gamma_j*x2_j
  # coeff_curr = coeff #coeff[,ii,]
  coeff_curr_rep = apply(coeff_curr, 1, function(x) rep_row(x, n_point_each))
  coeff_curr_3d = array(as.vector(coeff_curr_rep), dim = c(n_point_each, n_group, n_coeff))
  
  log_lambda = rowSums(coeff_curr_3d*X_3d, dim=2)   # log_lambda_ij = beta_j*t_ij + gamma_j*I_j
  
  # the constant term, - log_fac_approx(y), is removed,
  # because this term is cancelled out in log_post_new - log_post_current
  log_likeli = colSums(-exp(log_lambda) + y*log_lambda)
  
  # log of posteriror
  log_post = log_prior + log_likeli
  
  return(log_post)
}



# change the log dens group orig into a single function with coeff and hyper_param as inputs

# log_dens_group_approx_likeli =function(log_lambda_curr){
#   log_lambda_j = rowSums(log_lambda_curr)
#   log_joint_approx_likeli = sum(mapply(function(log_lambda,y) dnorm(log_lambda,log(y),sqrt(1/y)),
#                                        log_lambda_j, y_j))
#   log_joint_post_approx_likeli = log(prior_group_coeff_curr) + log_joint_approx_likeli
#   return(log_joint_post_approx_likeli)
# }
# 
# log_dens_group_approx_product_norm = function(coeff, mu, sigma,){
#   log_lambda_j = rowSums(log_lambda_curr) # alpha_j*1 + betaj*t_ij + gammaj*I_j
#   approx_product_norm = prod(coeff_group, log_lambda_j) # product of the group-level parameters alpha*gamma
#   mu_approx_product_norm = prod(mu,log(y_j))
#   var_approx_product_norm = prod(mu + sigma^2, log(y_j)^2 + 1/y_j) - mu_approx_product_norm^2
#   log_dens_approx_product_norm = dnorm(approx_product_norm, mu_approx_product_norm,
#                                        var_approx_product_norm, log = T)
# }



################################################################################
#### 2 - run the mh-in-gibbs
################################################################################

# tranform into a metropolis sample if there is no analytical solution due to non-conjugacy
# sample from normal distribution and check

# metropolis check (metropolis step) for each proposed sample to make sure that 
# the stationary distribution is the target distribution




## Initilization
n_coeff = 2
n_group = ncol(Y_mat)
group_coeff_id = 2 # coeff[c(1,3)]: alpha_j and gamma_j
indiv_coeff_id = seq(1, n_coeff)[-group_coeff_id]
n_group_coeff = length(group_coeff_id)
n_hyper_param = 2

n_keep = 30000
n_chain = 4
n_warm_up = n_keep/3
n_sample = n_warm_up + n_keep


# values of hyper params 
tau2_all = 1
a_all = 2
b_all = 2
a_sigma2_x = a_all + n_group


# for (i_chain in 1:n_chain){

meth_name = c('Original conditional','Conditional using approx. likelihood')

dens_fun_list = c(prop_norm) #, prop_approx_likeli)
n_meth = length(dens_fun_list)
  
hyper_param = array(NA, dim = c(n_coeff, n_sample, n_hyper_param, n_chain, n_meth))

coeff = array(NA, dim = c(n_coeff, n_sample, n_group, n_chain, n_meth))

accept_count_all = array(0, dim = c(n_coeff, n_group, n_chain, n_meth))

t_run = array(0, dim = c(n_meth))


for ( i_meth in 1:n_meth) {
  
  t_start  = Sys.time()

  for (i_chain in 1:n_chain){
    cat('\n\n\nChain: ', i_chain)
    
    # initialization
    hyper_param[, 1, , i_chain, i_meth] = array(abs(rnorm(1, 1, 1/2)), dim = c(n_coeff, n_hyper_param))
    coeff[, 1, , i_chain, i_meth] = array(rnorm(1, 0.001, 0.001/2), dim = c(n_coeff, n_group))
    accept_count = array(0, dim = c(n_coeff, n_group))
    
    for (i_samp in 2:n_sample){
      if (i_samp%%1000==0){
        cat('\n', i_samp)
      }
      
      coeff_hyper_param_return = gibbs_with_mh(
        coeff[, i_samp-1, , i_chain, i_meth], hyper_param[, i_samp-1, , i_chain, i_meth],
        log_post_orig, dens_fun_list[[i_meth]])
      
      coeff[, i_samp, , i_chain, i_meth] = as.matrix(coeff_hyper_param_return[[1]])
      
      hyper_param[, i_samp, , i_chain, i_meth] = as.matrix(coeff_hyper_param_return[[2]])
    }
    
    accept_count_all[,,i_chain, i_meth] = accept_count
    
  }

  t_end = Sys.time()
  
  
  t_run[i_meth] = t_end - t_start 
}


cat('Excecution time using Gibbs sampler: ', t_run)


################################################################################
##### 3 - Summarize and visualize posterior distributions 
################################################################################

plot_colors = c('blue','red', 'pink', "orange")
legend_names = sapply(1:n_chain, toString)
index_good = (20000):n_sample
# if(max(coeff)<=1e-4){
#   axis(2,at=marks,labels=format(marks,scientific=T))
# }
plot_sample = function(coeff, coeff_name) {
  for (i_meth in 1:n_meth) {
    for (i_coeff in 1:length(coeff_name)){
      par(mfrow=c(2, ceiling(n_group/2)))
      for (i_group in 1:n_group){
        for(i_chain in 1:n_chain){
          if (i_chain==1){
            plot(index_good, coeff[i_coeff, index_good, i_group,i_chain, i_meth], type='l',
                 xlab=paste(coeff_name[i_coeff],i_group), ylab='', col=plot_colors[i_chain])
          } else{
            lines(index_good, coeff[i_coeff, index_good, i_group,i_chain, i_meth],
                  xlab=paste(coeff_name[i_coeff],i_group), col=plot_colors[i_chain])
          }
        }
      }
      legend("bottom", legend=legend_names, col=plot_colors, lwd=1.0, 
             cex=1.0, bty='n', xpd =NA, horiz=T, inset=c(0,-0.5)) 
    }
  }
}

plot_density = function(coeff, coeff_name) {
  for (i_meth in 1:n_meth) {
    for (i_coeff in 1:length(coeff_name)){
      par(mfrow=c(2, ceiling(n_group/2)))
      for (i_group in 1:n_group){
        for(i_chain in 1:n_chain){
          if (i_chain==1){
            plot(density(coeff[i_coeff, index_good, i_group,i_chain, i_meth]), type='l', main = '',
                 xlab=paste(coeff_name[i_coeff],i_group), ylab='Density', col=plot_colors[i_chain])
          } else{
            lines(density(coeff[i_coeff, index_good, i_group,i_chain, i_meth]), main = '',
                  xlab=paste(coeff_name[i_coeff],i_group), ylab='Density', col=plot_colors[i_chain])
          }
        }
      }
      legend("bottom", legend=legend_names, col=plot_colors, lwd=1.0, 
             cex=1.0, bty='n', xpd = NA, horiz = T, inset = c(0,-0.5)) 
    }
  }
}

coeff_name = c('beta ','gamma ')
plot_sample(coeff, coeff_name)
plot_density(coeff, coeff_name)


for (i_meth in 1:n_meth){
  for (i_coeff in 1:n_coeff) {
    cat('\nAcceptance rate of', paste(coeff_name[i_coeff]), ':',
        rowSums(round(accept_count_all[i_coeff, , ,i_meth]/(n_chain*n_sample), 3)))
  }
}

for (i_meth in 1:n_meth){
  for (jj in 1:length(coeff_name)) {
    for (ii in 1:n_group){
      # print(paste('Posterior of', coeff_name[jj]))
      cat('\n', mean(coeff[jj, index_good, ii, , i_meth]))
    }
    # mtext(paste('Posterior of', coeff_name[jj]), outer=TRUE,  cex=1, line=-0.5)
  }
}

glm_fit = glm(as.vector(y) ~ as.vector(X_3d[,,1]) + as.vector(X_3d[,,2]),family = 'poisson')
cat('Poisson regression:\n', glm_fit$coefficients)



# hyper_param_name = c('mu','sigma')
# for (jj in 1:length(hyper_param_name)){
#   par(mfrow=c(ceiling(n_coeff/2),2))
#   for (ii in 1:n_coeff){
#     plot(hyper_param[ii, (n_warm_up+1):n_sample, jj], ylab = paste(hyper_param_name[jj], ii))
#   }
#   # mtext(0.5,0.55,paste('Posterior of', hyper_param_name[jj]), outer=TRUE,  cex=1, line=-0.5)
#   print(mean(hyper_param[, (n_warm_up+1):n_sample, jj]))
# }

# change sequence of coeff. Checked. Same results.
# check post of hyper parameters

# reduce the problem to a simpler problem, say, mh in gibbs for two parameters without hyper parameters


################################################################################
#### 34 - Assess bias and coverage
################################################################################

# Autocorrrelatoin vs lag
# Posterior summaries
# https://support.sas.com/rnd/app/stat/examples/BayesSalm/new_example/index.html


# posterior summaries
# effect size
# https://www.flutterbys.com.au/stats/tut/tut7.2b.html


# useful ref.
# https://github.com/stablemarkets/BayesianTutorials/blob/master/MultipleLinearReg/multiplelinearreg.r


# Hamiltonian MC
## Ref.: https://datascienceplus.com/bayesian-regression-with-stan-part-1-normal-regression/
##       https://jrnold.github.io/bayesian_notes/introduction-to-stan-and-linear-regression.html
##       https://rpubs.com/kaz_yos/stan-pois1


# # alpha = min(tar(x_prop)*prop(x_curr|x_prop)/tar(x_curr)prop(x_prop|x_curr), 1)
# # if the proposal distribution is normal, symmetrical, then prop(x|y) = prop(y|x), alpha = min(tar(x_prop)/tar(x_curr), 1)
# # e.g. N(x2|x1, sigma2) = N(x1|x2, sigma2), proportional to exp(-0.5*|x1-x2|^2/sigma2)
# 
# 
# sample_x = mh_draw(start_val, n_sample)
# n_burn = 5000
# acceptance = 1-mean(duplicated(chain[-(1:n_burn),]))

# cache in MH to boost the efficiency
# ref.: https://www.cra.com/Figaro_ScalaDoc/com/cra/figaro/library/cache/MHCache.html/
#     : https://stablemarkets.wordpress.com/2019/03/02/efficient-mcmc-with-caching/  
#     : https://theoreticalecology.wordpress.com/2010/09/17/metropolis-hastings-mcmc-in-r/
