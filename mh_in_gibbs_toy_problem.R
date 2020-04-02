# ref.: https://www.r-bloggers.com/quick-illustration-of-metropolis-and-metropolis-in-gibbs-sampling-in-r/
## Metropolis sampling
## x    - current value of Markov chain (numeric vector)
## log_targ - log of target density function
## prop - function with prototype function(x, ...) that generates 
##   a proposal value from a symmetric proposal distribution

# synthesize data
generate_one_dim_data = function(n_point_each){
  set.seed(2)
  x = abs(round(sort(runif(n_point_each,1,15)),2))
  beta = 0.5
  gamma = 1
  y = round(exp(x*beta+gamma), 0)
  return(list(y,x,beta,gamma))
}

# proposal function
prop_norm = function(x, x_id){
  # input: a single numeric value of current sample
  # output; a single numeric value of new sample
  if (x_id==1){
    sd = 0.0015
  }else{
    sd = 0.02
  }
  return(rnorm(1, mean=x, sd=sd))
}

# metropolis draw
mh_draw = function(x, x_id, log_targ, prop=prop_norm) {
  # input: x - a vector of current samples
  # output: a vector of new samples
  
  x_new = x
  # only update one of the params
  x_new[x_id] = prop(x[x_id], x_id)
  
  log_rate = log_targ(x_new, x_id) - log_targ(x, x_id)
  if (log(runif(1))<log_rate) {
    x = x_new
    counter[x_id] <<- counter[x_id] + 1
  }
  
  return(x)
}

## Metropolis-in-Gibbs sampling
## x - current vector of samples
## log_targ - target log density function
mh_gibbs_draw = function(x, log_targ){
  for (ii in 1:length(x)){
    x_new = mh_draw(x, ii, log_targ, prop=prop_norm)
    x = x_new
  }
  return(x)
}


################################################

log_targ_my = function(x, x_id){
  
  log_prior = dnorm(x[x_id], 0, 2, log=T)  # N(0,2), a weakly-informative prior
  
  log_lambda = x[1]*X_data + x[2]
  
  # the constant term,log_fac_approx(y_data), is removed
  log_likeli = sum(-exp(log_lambda) + y_data*log_lambda)
  
  # log posterior
  log_post = log_prior + log_likeli
  
  return(log_post)
  
}

## generate one dimensional data
n_point_each = 12
toy_data = generate_one_dim_data(n_point_each)
y_data = c(unlist(toy_data[1]))
cat('y data',y_data)
X_data = c(unlist(toy_data[2]))
cat('x data',X_data)


## initialize and sample using Metropolis-in-Gibbs
n_spl = 10000
n_good = 1000

# initialize
xcur = c(0.5,0.5)
x_spl = matrix(NA, n_spl, 2)
counter = c(0,0)

# plot(0,0, xlim = c(-0.5,2), ylim = c(-0.5,2), cex = 0.1,
#      xlab=expression(beta), ylab=expression(gamma), main="Metropolis-in-Gibbs")

# sample and plot the samples
for(j in 1:n_spl) {
  if ( j%%500 == 0) {
    print(j)
  }
  xcur = mh_gibbs_draw(xcur, log_targ_my)
  x_spl[j,] = xcur
  # if (j>=n_good){
  #   points(xcur[1], xcur[2], pch=20, col='#00000055')
  # }
}

par(mfrow=c(1,2))
plot(n_good:n_spl,x_spl[n_good:n_spl,1], xlab = 'Index', ylab = expression(beta))
plot(n_good:n_spl,x_spl[n_good:n_spl,2], xlab = 'Index', ylab = expression(gamma))


par(mfrow=c(1,2))
plot(density(x_spl[n_good:n_spl,1]), xlab = expression(beta), main = 'Posterior PDF')
plot(density(x_spl[n_good:n_spl,2]), xlab = expression(gamma), main = 'Posterior PDF')

# compare against poisson glm
cat('Metropolos in Gibbs:\n', colMeans(x_spl[n_good:n_spl,]))

glm_fit = glm(y_data~X_data,family = 'poisson')
cat('Poisson regression:\n', glm_fit$coefficients)

# accept rate
cat('Acceptance rate', counter/n_spl)
