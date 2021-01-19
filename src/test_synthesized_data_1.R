# import packages
library('tidyr')
library('rstan')

# set directory
setwd('C:/GitHub/approx_gibbs_for_HBM')

# dependent source code
source('data_preprocessing.R')
source('gibbs_all_vectorized.R')
source('diagnostics_and_plotting.R')

# set accuracy
options(digits = 5)

  # # sampling parameters
  n_keep = 8000L
  n_warmup = 1000L
  n_chain = 2L
  
  fit_model(data_all, group_attr_id)
  
  # metric_all_mat = matrix(NA,n_case,n_metric*2)
  # metric_all_mat[i,] = fit_model(syn_data, group_attr_id)
  # 
  # cat('\n current iter.: ', i)
  # 
  # write.csv(metric_all_mat, "metric_all_mat.csv")
  # save(metric_all_mat, file="metric_all_mat.RData")
  # load(metric_all_mat.RData)
  
# }




# # ## recovery data
# data_list = read.csv("C:/GitHub/approx_gibbs_for_HBM/power_recovery_data_vectorize.csv", header = T)
# group_attr_id = 2
# data_all_orig = prepare_data(data_list, group_attr_id)
# data_stan = data_all_orig[[1]]
# data_gibbs = data_all_orig[[2]]
# # normalization
# for (i in 1:group_attr_id) {
#   data_stan$X[,i] = min_max(data_stan$X[,i],lb_norm=0.005)
#   data_gibbs$X[,i] = min_max(data_gibbs$X[,i],lb_norm=0.005)
# }
# 
# data_all = list(data_stan, data_gibbs)


# stan
# stan_fit = fit_stan(file_stan='1d_HBM_multi_attri.stan', data_stan, n_keep, n_warmup, n_chain)
# stan_fit = stan(file='1d_HBM_multi_attri.stan', data = data_stan, iter = n_keep+n_warmup, warmup=n_warmup, chains = n_chain)

# post processing
# rstan::traceplot(stan_fit, pars='coeff')
# print(stan_fit, digits_summary = 3)
# 
 
 
# param_post_sampl = extract(stan_fit, pars=c('coeff','mu','sigma2'), permuted = TRUE)


# gibbs algorithms
y = data_gibbs$y
exact_mean = sapply(y, cal_exact_mean)
exact_var = sapply(y, cal_exact_var)
gibbs_fit = fit_gibbs(data_gibbs, n_keep, n_warmup, n_chain)
index_good = (n_warmup+1):(n_warmup+n_keep)
# coeff_gibbs = gibbs_fit[[1]]
# hyper_param = gibbs_fit[[2]]
# accept_count_all = gibbs_fit[[3]]

# index_good = (n_warmup+1):(n_warmup+n_keep)



# plot_colors <<- c('blue','red', 'purple', "orange2")
# legend_names = sapply(1:n_chain, toString)
# n_total_sample = n_warmup+n_keep
# index_good = (n_total_sample*0.5):n_total_sample
# if(max(coeff)<=1e-4){
#   axis(2,at=marks,labels=format(marks,scientific=T))
# }

# def_layout = function(n_group){
# 
#   if(n_group%%6==0){
#     m <- matrix(c(seq(1,6,by=1),rep(7,3)),nrow = 3,ncol = 3,byrow = TRUE)
#   }else if(n_group%%8==0){
#     m <- matrix(seq(1,9,by=1),nrow = 3,ncol = 3,byrow = TRUE)
#   }
#   layout(mat = m,heights = c(0.4,0.4,0.2))
# }

# plot_sample = function(coeff, coeff_name) {
#   for (i_meth in 1:n_meth) {
#     n_row = floor((length(coeff_name)*n_group)/6)    # 6/2)
#     m <- matrix(c(seq(1,n_row*6)), nrow = n_row, ncol = 6, byrow = TRUE)
#     layout(mat = m,heights = rep(1/n_row, n_row))
# 
#     for (i_coeff in 1:length(coeff_name)){
# 
#       for (i_group in 1:n_group){
#         par(mar = c(2,4,2,1))
#         for(i_chain in 1:n_chain){
#           if (i_chain==1){
#             plot(index_good, coeff[i_coeff, index_good, i_group,i_chain, i_meth], type='l',
#                  ylab=paste(coeff_name[i_coeff],i_group), xlab='', col=plot_colors[i_chain])
#           }else{
#             lines(index_good, coeff[i_coeff, index_good, i_group,i_chain, i_meth],
#                   ylab=paste(coeff_name[i_coeff],i_group), xlab= '', col=plot_colors[i_chain])
#           }
#         }
#         # if (i_group%%6==0){
#         #   par(mar = c(1,1,2,1))
#         #   plot(1, type = "n", axes=FALSE, xlab="", ylab="")
#         #   legend(legend=legend_names, title = 'Chain', col=plot_colors,
#         #          x='top', bty='n',inset=0, lwd=1.5,cex=1.1, horiz=T)
#         # }
#       }
#     }
#   }
# }
# coeff_name = paste('Cov_', 1:n_cov)
# plot_sample(coeff, coeff_name)

# # coeff dim = c(n_coeff, n_total_sample, n_group, n_chain, n_meth))
# n_coeff = dim(coeff)[1]
# n_group = dim(coeff)[3]
# n_meth = dim(coeff)[5]
# 
# coeff_mean = array(NA,dim=c(n_coeff,n_group, n_meth))
# for (i_coeff in 1:n_coeff) {
#   for (i_group in 1:n_group) {
#     for (i_meth in 1:n_meth) {
#       coeff_mean[i_coeff, i_group, i_meth] = mean(coeff[i_coeff,index_good, i_group, ,i_meth]) 
#     }
#   }
# }

# print(coeff_mean)

# coeff_name = paste('Cov_', 1:n_cov)
# plot_sample(coeff, coeff_name)

# for (i_meth in 1:n_meth){
#   for (i_coeff in 1:n_coeff) {
#     cat('\nAcceptance rate of', paste(coeff_name[i_coeff]), ':',
#         rowSums(round(accept_count_all[i_coeff, , ,i_meth]/(n_chain*n_total_sample), 3)),'\n')
#   }
# }

# print(coeff_mean)
# print(stan_fit, pars='coeff', digits_summary = 5)

metric_all = cal_metric(stan_fit, gibbs_fit, y, print_out=1)
  