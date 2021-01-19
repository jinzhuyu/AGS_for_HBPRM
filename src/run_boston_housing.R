################
# covid_19 test data

################
## import and adapt data
# import data
data_list = read.csv("C:/GitHub/approx_gibbs_for_HBM/test_data/boston_housing_data/boston_housing.csv", header = T)

# select useful columns
data_list = data_list[, c(2,3,4, 6:15)]

# add a group-level covariate~norm(0,1)
# find group id
group_attr_id_temp = 8
n_data = length(data_list[,1])
group_id = rep(1,n_data)
group_id_temp = 1
cov_group_rep = rep(abs(rnorm(1,0,0.5)), n_data)
for(i in 2:n_data){
  if(data_list[i,group_attr_id_temp]==data_list[i-1,group_attr_id_temp]){
    group_id[i] = group_id_temp
    cov_group_rep[i] = cov_group_rep[i-1]
    i = i+1
  }else{
    group_id_temp = group_id_temp+1
    group_id[i] = group_id_temp
    cov_group_rep[i] = abs(rnorm(1,0,0.5))
  }
}

# add covariates to data
data_list = cbind(data_list, cov_group_rep)

# reorder the columns
n_col = dim(data_list)[2]
selected_col = c(1,3:(group_attr_id_temp-1),(group_attr_id_temp+1):(n_col-2),n_col,n_col-1)
data_list_select = data_list[,selected_col]

# medv should be integer. Original price unit is 10^3
data_list_select$medv = data_list_select$medv*10


# data features
n_cov = length(data_list_select[1,])-1
group_attr_id = n_cov
n_group = max(group_id)
len_each_group = find_len_each_group(data_list_select, group_attr_id, n_data, n_group)


data_all = prepare_data(data_list_select, group_attr_id)

# # sampling parameters
n_keep = 8000L
n_warmup = 4000L
n_chain = 2L

# fit model
data_stan <<- data_all[[1]]
data_gibbs <<- data_all[[2]]

# stan
stan_fit <<- stan(file='1d_HBM_multi_attri.stan', data = data_stan, iter = n_keep+n_warmup, warmup=n_warmup, chains = n_chain)

# # gibbs
y <<- data_gibbs$y    # mean=225.33, var=8458.7
exact_mean <<- sapply(y, cal_exact_mean)
exact_var <<- sapply(y, cal_exact_var)
gibbs_fit <<- fit_gibbs(data_gibbs, n_keep, n_warmup, n_chain)
index_good <<- (n_warmup+1):(n_warmup+n_keep)

# value of performance metrics
metric_all = cal_metric(stan_fit, gibbs_fit, y, print_out=1)
metric_all