# Inference and prediction in Hierarhical Bayesian model with application in predicting recovery from power outage

# rm(list = ls()) #clear variables

# cat("\014")  #clear console

library('rstan')
library('coda')
library('MCMCpack')
library('invgamma')

setwd('C:/GitHub/approx_gibbs_for_HBM')

# load data
data_all_0 = read.csv("C:/GitHub/approx_gibbs_for_HBM/power_recovery_data.csv", header = F)

data_all=as.matrix(data_all_0)

Y_mat=data_all[2:nrow(data_all),6:10]
time_mat=data_all[2:nrow(data_all),1:5]
intensity=data_all[1,6:10]

  
  
n_group=ncol(Y_mat)
n_time=nrow(time_mat)

# Hamiltonian MC
## Ref.: https://datascienceplus.com/bayesian-regression-with-stan-part-1-normal-regression/
##       https://jrnold.github.io/bayesian_notes/introduction-to-stan-and-linear-regression.html
##       https://rpubs.com/kaz_yos/stan-pois1



set.seed(1)

train=1:(nrow(Y_mat)-2)

# HB_data=list(Y=Y_mat[train,],I=(nrow(Y_mat)-2),J=n_group,x1=time_mat[train,],x2=intensity)
# ,
             # I_new=2,J_new=n_group,
             # x1_new=time_mat[-train,])#,x2_new=intensity)

Y_vec = array(as.numeric(unlist(Y_mat[train,])))
time_vec = array(as.numeric(unlist(time_mat[train,])))
intensity_vec = rep(array(as.numeric(unlist(intensity))),each=length(train))/1000

index_group_vec = rep(1:n_group,each=length(train))
len_vec = length(Y_vec)

X_vec = rbind(time_vec, intensity_vec)

HB_data_1d = list(Y = Y_vec, I=length(train), J=n_group, N=length(train)*n_group,
                  X = X_vec, group_id=index_group_vec, K=nrow(X_vec))
# time = time_vec, intensity = intensity_vec

# glm_fit = glm(Y_vec~time_vec+intensity_vec1, family = "poisson")
# summary(glm_fit)

n_sam = 12000
n_warmup = round(n_sam/3,0)
n_chain = 4

HBR_1d_fit=stan(file='1d_HBM_multi_attri.stan',data = HB_data_1d,iter =(n_sam+n_warmup), warmup=n_warmup,
                chains = n_chain)#, control = list(adapt_delta = 0.85, max_treedepth = 15))

rstan::traceplot(HBR_1d_fit, pars='coeff')
print(HBR_1d_fit, pars='coeff')

rstan::traceplot(HBR_1d_fit, pars=c('mu','sigma2'))
print(HBR_1d_fit, pars=c('mu','sigma2'))


# print(HBR_1d_fit, pars=c('beta','gamma'))
# rstan::traceplot(HBR_1d_fit, pars=c('beta','gamma'))


# trace_beta1=rstan::traceplot(HBR_1d_fit,pars=c('beta'))
# trace_beta1+ scale_color_discrete() + theme(legend.position = c(0.85, 0.17)) # shows up clearly when zoomed in
# trace_alpha1+ theme(legend.position = c(0.85, 0.2))
# rstan::traceplot(HBR_1d_fit,pars=c('alpha'))
# rstan::traceplot(HBR_1d_fit,pars=c('beta'))
# rstan::traceplot(HBR_1d_fit,pars=c('gamma'))

# }

# extract posteriors

HBR_fit=extract(HBR_1d_fit, permuted = TRUE)


# sampler diagnostics
# ref.: https://cran.r-project.org/web/packages/rstan/vignettes/stanfit-objects.html

sampler_params <- get_sampler_params(HBR_1d_fit, inc_warmup = FALSE)
mean_accept_stat_by_chain <- sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))
print(mean_accept_stat_by_chain)

# convergence diagonostics
# https://rdrr.io/cran/rstan/man/Rhat.html
#If chains have not mixed well (ie, the between- and within-chain estimates don't agree), 
# R-hat is larger than 1. 
# We recommend running at least four chains by default and 
# only using the sample if R-hat is less than 1.05. 

##  Misbehaving chains indicate bias of seemingly good chains.

# Rhat(HBR_1d_fit)

### excecution time
time_stan = get_elapsed_time(HBR_1d_fit)
time_stan_total = sum(time_stan[,2])
cat('Time for each chain of samples with Stan: ', time_stan[,2])
cat('Total time for 4 chain of samples with Stan: ', time_stan_total)

### effective sample size
# summary(HBR_1d_fit)$summary[, "n_eff"]

## Traceplot and posterior distribution
#### plotting the posterior distribution for the parameters
# post_beta=As.mcmc.list(m_norm,pars="beta")
# plot(post_beta[,1])


#### plot the correlation between the parameters
# pairs(HBR_1d_fit,pars="beta")

# plot credible intervals for the different betas
# plot(HBR_1d_fit,pars=c("alpha"))

trance_plot=plot(HBR_1d_fit, plotfun = "trace", pars = c('coeff'), inc_warmup = TRUE)
trance_plot + scale_color_discrete() + theme(legend.position = "top")

#### fitted curves
# #getting regression curves plus 95% credible intervals
# new_x<-data.frame(x1=new_X[,2],x2=rep(c("Min","Mean","Max"),each=20))
# new_y<-extract(m_norm,pars="y_pred")
# pred<-apply(new_y[[1]],2,quantile,probs=c(0.025,0.5,0.975)) #the median line with 95% credible intervals
# #plot
# plot(dat$x1,y_norm,pch=16)
# lines(new_x$x1[1:20],pred[2,1:20],col="red",lwd=3)
# lines(new_x$x1[1:20],pred[2,21:40],col="orange",lwd=3)
# lines(new_x$x1[1:20],pred[2,41:60],col="blue",lwd=3)
# lines(new_x$x1[1:20],pred[1,1:20],col="red",lwd=1,lty=2)
# lines(new_x$x1[1:20],pred[1,21:40],col="orange",lwd=1,lty=2)
# lines(new_x$x1[1:20],pred[1,41:60],col="blue",lwd=1,lty=2)
# lines(new_x$x1[1:20],pred[3,1:20],col="red",lwd=1,lty=2)
# lines(new_x$x1[1:20],pred[3,21:40],col="orange",lwd=1,lty=2)
# lines(new_x$x1[1:20],pred[3,41:60],col="blue",lwd=1,lty=2)
# legend("topright",legend=c("Min","Mean","Max"),lty=1,col=c("red","orange","blue"),bty = "n",title = "Effect of x2 value on\nthe regression")
















#-----------------------------

# Gibbs sampling with approximate closed-form conditionalf of priors

## Hyperpriors
#### vague normal
mu_0=0
sigma_0_square=0.5
#### vague inverse gamma_mat
a_0=2
b_0=2


## Complete conditional posterior distributions

# parameters for hyperpriors
tau2_all=1

a_all=2
b_all=2

rep_row<-function(x,n){
   matrix(rep(x,each=n),nrow=n)
}

a_sigma2_all=a_all+n_group

draw.mu_alpha=function(alpha,sigma2_alpha){
  mu_mu_alpha=tau2_all*sum(alpha)/(n_group*tau2_all+sigma2_alpha)
  sigma2_mu_alpha=tau2_all*sigma2_alpha/(n_group*tau2_all+sigma2_alpha)
  return(rnorm(1,mu_mu_alpha,sqrt(sigma2_mu_alpha)))
}

draw.mu_beta=function(beta,sigma2_beta){
  mu_mu_beta=tau2_all*sum(beta)/(n_group*tau2_all+sigma2_beta)
  sigma2_mu_beta=tau2_all*sigma2_beta/(n_group*tau2_all+sigma2_beta)
  return(rnorm(1,mu_mu_beta,sqrt(sigma2_mu_beta)))
}

draw.mu_gamma=function(gamma,sigma2_gamma){
  mu_mu_gamma=tau2_all*sum(gamma)/(n_group*tau2_all+sigma2_gamma)
  sigma2_mu_gamma=tau2_all*sigma2_gamma/(n_group*tau2_all+sigma2_gamma)
  return(rnorm(1,mu_mu_gamma,sqrt(sigma2_mu_gamma)))
}

draw.sigma2_alpha=function(alpha,mu_alpha){
  # a_sigma2_alpha=a_all+n_group
  b_sigma2_alpha=b_all+sum((alpha-mu_alpha)^2)/2
  return(rinvgamma(1,a_sigma2_all,b_sigma2_alpha))
}

draw.sigma2_beta=function(beta,mu_beta){
  # a_sigma2_beta=a_all+n_group
  b_sigma2_beta=b_all+sum((beta-mu_beta)^2)/2
  return(rinvgamma(1,a_sigma2_all,b_sigma2_beta))
}

draw.sigma2_gamma=function(gamma,mu_gamma){
  # a_sigma2_gamma=a_all+n_group
  b_sigma2_gamma=b_all+sum((gamma-mu_gamma)^2)/2
  return(rinvgamma(1,a_sigma2_all,b_sigma2_gamma))
}

# regression coefficients
# matrix(rep(vec,each=n),nrow=n)

draw.alpha=function(beta,gamma,mu_alpha,sigma2_alpha){
  mu_alpha_hat=(mu_alpha+
                sigma2_alpha*colSums(y_mat_gibbs*
                (log(y_mat_gibbs)-rep_row(beta,n_row_train)*time_gibbs_mat-
                   rep_row(gamma*intensity_gibbs,n_row_train))))/
                (sigma2_alpha*sum_y+1)
  sigma2_alpha_hat=1/(sum_y+1/sigma2_alpha)
  return(rep(0, n_group))
}

draw.beta=function(alpha,gamma,mu_beta,sigma2_beta){
  mu_beta_hat=(mu_beta+
               sigma2_beta*colSums(y_mat_gibbs*time_gibbs_mat*
                (log(y_mat_gibbs)-rep_row(alpha,n_row_train)-
                   rep_row(gamma*intensity_gibbs,n_row_train))))/
                (sigma2_beta*sum_y_time2+1)
  sigma2_beta_hat=1/(sum_y_time2+1/sigma2_beta)
  return(rnorm(n_group,mu_beta_hat,sqrt(sigma2_beta_hat)))
}

draw.gamma=function(alpha,beta,mu_gamma,sigma2_gamma){
  mu_gamma_hat=(mu_gamma+
                sigma2_gamma*colSums(y_mat_gibbs*rep_row(intensity_gibbs,n_row_train)*
                (log(y_mat_gibbs)-rep_row(alpha,n_row_train)-
                   rep_row(beta,n_row_train)*time_gibbs_mat)))/
                (sigma2_gamma*sum_y_inten2+1)
  sigma2_gamma_hat=1/(sum_y_inten2+1/sigma2_gamma)  
  return(rnorm(n_group,mu_gamma_hat,sqrt(sigma2_gamma_hat))) #!!!
}

# alpha=alpha_mat[i-1,]
# beta=beta_mat[i-1,]
# mu_gamma=mu[i,3]
# sigma2_gamma=sigma2[i,3]

# draw.gamma(alpha_mat[i-1,],beta_mat[i-1,],mu[i,3],sigma2[i,3])

## Sampling parameters
n_burn=n_warmup
n_good=n_sam
n_totl=n_burn+n_good
n_coeff=3

## Initialization of coefficients and hyperparameters

alpha_mat=array(NA,dim=c(n_totl,n_group,n_chain))
beta_mat=array(NA,dim=c(n_totl,n_group,n_chain))
gamma_mat=array(NA,dim=c(n_totl,n_group,n_chain))

mu=array(NA,dim=c(n_totl,n_coeff,n_chain))
sigma2=array(NA,dim=c(n_totl,n_coeff,n_chain))

y_mat_gibbs=matrix(as.numeric(unlist(Y_mat[train,])),length(train),n_group)
time_gibbs_mat=matrix(as.numeric(unlist(time_mat[train,])),length(train),n_group)
intensity_gibbs=c(as.numeric(unlist(intensity)))

n_row_train=nrow(y_mat_gibbs)
sum_y_inten2=colSums(y_mat_gibbs*rep_row(intensity_gibbs,n_row_train)^2)
sum_y_time2=colSums(y_mat_gibbs*time_gibbs_mat^2)
sum_y=colSums(y_mat_gibbs)

## value of first sample

for (i_chain in 1:n_chain){
  init_value = abs(rnorm(1, mean=1, sd=1))
  mu[1,,i_chain]=rep(init_value,n_coeff)
  sigma2[1,,i_chain]=rep(init_value,n_coeff)   
  gamma_mat[1,,i_chain]= rep(init_value,n_group)
  beta_mat[1,,i_chain]= rep(init_value,n_group)
  alpha_mat[1,,i_chain]=rep(0,n_group)
  
    
  ## Sampling
  # t_start_gibbs=Sys.time()
  for (i in 2:n_totl){
    # update hyperpriors
    # print('iter #')
    # print(i)
    
    mu[i,1,i_chain]=draw.mu_alpha(alpha_mat[i-1,,i_chain],sigma2[i-1,,i_chain])
    mu[i,2,i_chain]=draw.mu_beta(beta_mat[i-1,,i_chain],sigma2[i-1,,i_chain])
    mu[i,3,i_chain]=draw.mu_gamma(gamma_mat[i-1,,i_chain],sigma2[i-1,,i_chain])
    # print('mu')
    # print(mu[i,])
    
    sigma2[i,1,i_chain]=draw.sigma2_alpha(alpha_mat[i-1,,i_chain],mu[i,1,i_chain])
    sigma2[i,2,i_chain]=draw.sigma2_beta(beta_mat[i-1,,i_chain],mu[i,2,i_chain])
    sigma2[i,3,i_chain]=draw.sigma2_gamma(gamma_mat[i-1,,i_chain],mu[i,3,i_chain])
   
    # print('sigma2')
    # print(sigma2[i,])
    
    gamma_mat[i,,i_chain] = draw.gamma(alpha_mat[i-1,,i_chain],beta_mat[i-1,,i_chain],mu[i,3,i_chain],sigma2[i,3,i_chain])
    beta_mat[i,,i_chain] = draw.beta(alpha_mat[i-1,,i_chain],gamma_mat[i,,i_chain],mu[i,2,i_chain],sigma2[i,2,i_chain]) 
    alpha_mat[i,,i_chain] = draw.alpha(beta_mat[i,,i_chain],gamma_mat[i,,i_chain],mu[i,1,i_chain],sigma2[i,1,i_chain]) 
    
    # print('gamma')
    # print(gamma_mat[i,])
    # print('beta')
    # print(beta_mat[i,])
    # print('alpha')
    # print(alpha_mat[i,])
    # update priors
    
    # update regression coefficients
  }
}
# t_end_gibbs = Sys.time()


# t_gibbs=t_end_gibbs-t_start_gibbs

# cat('Excecution time using Gibbs sampler: ', t_gibbs)


### Abandon burn-in samples
index_good = (n_burn+1):n_totl
  
mu_good = mu[index_good,,]
sigma2_good = sigma2[index_good,,]
gamma_mat_good = gamma_mat[index_good,,]
beta_mat_good = beta_mat[index_good,,]
alpha_mat_good = alpha_mat[index_good,,]


### mean values

# rowMeans(colMeans(alpha_mat_good))
# 
rowMeans(colMeans(beta_mat_good))
# 
rowMeans(colMeans(gamma_mat_good))


print(HBR_1d_fit,pars=c('coeff'))

plot_colors = c('blue','red', 'pink', "orange")
legend_names = c('1','2', '3', '4')
plot_post_coeff = function(coeff_curr_3d, coeff_name){
  par(mfrow=c(2, ceiling(n_group/2)))
  for (ii in 1:n_group){
    for(i_chain in 1:n_chain){
      if (i_chain==1){
        plot(index_good, coeff_curr_3d[index_good,ii,i_chain], type='l', xlab='', ylab='',
             main=paste(coeff_name,ii), col=plot_colors[i_chain])
      } else{
        lines(index_good, coeff_curr_3d[index_good,ii,i_chain], col=plot_colors[i_chain])
      }
    }
  }
  legend("bottom", legend=legend_names, col=plot_colors, lwd=1.0, cex=1.0, bty='n', xpd = NA, horiz = T, inset = c(0,-0.5))
}
# improve legend: https://datascienceplus.com/mastering-r-plot-part-3-outer-margins/

# par(mfrow=c(2,3))
    # oma = c(4,3,0,0) + 2,# two rows of text at the outer left and bottom margin
    # mar = rep(2,4) + 2) # space for one row of text at ticks and to separate plots
### trace plot

plot_post_hyper = function(hyper_param_curr, hyper_param_name){
  par(mfrow=c(2, ceiling(n_coeff/2)))
  for (ii in 1:n_coeff){
    for(i_chain in 1:n_chain){
      if (i_chain==1){
        plot(index_good, hyper_param_curr[index_good,ii,i_chain], type='l', xlab='', ylab='',
             main=paste(hyper_param_name,ii), col=plot_colors[i_chain])
      } else{
        lines(index_good, hyper_param_curr[index_good,ii,i_chain], col=plot_colors[i_chain])
      }
    }
  }
  legend("bottom", legend=legend_names, col=plot_colors, lwd=1.0, cex=1.0, bty='n',
         xpd = NA, horiz = T, inset = c(0,-0.5))
}


# coeff_nam = c('beta ', 'gamma ')

plot_post_coeff(beta_mat, 'beta ')
plot_post_coeff(gamma_mat, 'gamma ')

plot_post_hyper(mu, 'mu')
plot_post_hyper(sigma2, 'sigma')


# # sigma2
# 
# par(mfrow=c(1,3))
# for (i in 1:3){
#   plot(index_good, sigma2[index_good,i],type='l',xlab='',ylab='',main=paste('sigma2 ',i),col=i)
# }


### effective sample size

# https://cran.r-project.org/web/packages/mcmcse/mcmcse.pdf mcmcse package

# lapply(alpha_mat,effectiveSize) 
# work with an mcmc or mcmc.list object. To get the size for each chain individually use lapply(x,effectiveSize)

## Approximate posterior predictive distribution
### negative binomial in r: The negative binomial distribution with size = n and prob = p has density
# ??(x+n)/(??(n) x!) p^n (1-p)^x, for x = 0, 1, 2, ., n > 0 and 0 < p ??? 1.
# The mean is ?? = n(1-p)/p and variance n(1-p)/p^2.

# a_hat=
# b_hat=
# # mu_mu_hat=
# # sigma_mu_hat=
# y_hat=rnbinom(n_good,size=a_hat,prob=1/(1+b_hat))



# R-INLA
# https://www.precision-analytics.ca/blog-1/inla
# http://www.maths.bath.ac.uk/~jjf23/brinla/
# https://stefansiegert.net/inla-project/inla-from-scratch  from scratch