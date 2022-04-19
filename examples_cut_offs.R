library(truncnorm)
library(mvtnorm)
library(tmvtnorm)

vec = function(M){
  
  return(as.vector(M))
  
}

rmatnorm = function(M, U, V, tol = .Machine$double.eps^0.5, method = "chol"){
  n = nrow(M)
  p = ncol(M)
  Sigma = kronecker(U, V)
  vec.X = mvtnorm::rmvnorm(1, vec(M), Sigma, method = method)
  X = matrix(vec.X, nrow = n, ncol = p, dimnames = list(rownames(U),colnames(V)))
  return(X)
}


# Dataset

p = 5 # number of variables
K = 2 # number of levels of each variable (assumed to be equal for all)

n = 20 # number of observations

set.seed(1)

X = matrix(NA, n, p)

for(i in 1:n){
  X[i,] = sample(1:K, p, replace = TRUE)
}

X # a (n,p) dataset

set.seed(1)

# construct the matrix parameter Theta collecting the (true) cutoffs
# for each variable (column of Theta) we have K+1 cutoffs
# for each variable the first and last cutoffs are fixed beforehand at -Inf and Inf respectively
# finally, we have K-1 "free" cutoffs for each variable

Theta = matrix(NA, K+1, p)

for(j in 1:p){
  Theta[,j] = c(-Inf, sort(rtruncnorm(K-1, a = -Inf, b = Inf)), Inf)
}

Theta

###########################
## Dealing with cut-offs ##
###########################

source("functions_match_cut_offs.r")

## see functions_match_cut_offs.r for details on the various functions


match_x_Theta(X[1,], Theta, p) # cutoffs of individual 1 (each row corresponds to a variable)
match_x_Theta(X[2,], Theta, p)


cbind(match_x_theta_j(X[,1], Theta[,1]), X[,1]) # apply the function to variable 1
cbind(match_x_theta_j(X[,2], Theta[,2]), X[,2])
cbind(match_x_theta_j(X[,3], Theta[,3]), X[,3])


match_x_theta_j_lower(X[,1], Theta[,1]) # lower cut-offs for variable 1 (and all individuals)
match_x_theta_j_upper(X[,1], Theta[,1]) # upper cut-offs for variable 1 (and all individuals)


######################################
## Dealing with truncnorm functions ##
######################################

# To propose a candidate value for the cutoff we will need

x = rtruncnorm(1, mean = 0, sd = 1, a = 0, b = 2) # normal N(mean, sd) truncated at [a,b]
x

# To evaluate the proposal density of a cutoff

dtruncnorm(x, mean = 0, sd = 1, a = 0, b = 2) # normal N(mean, sd) truncated at [a,b]

# To evaluate the likelihood

ptruncnorm(x, mean = 0, sd = 1, a = 0, b = 2) # normal N(mean, sd) truncated at [a,b]



## Storage for the cutoffs

S = 10 # will be the number of MCMC iterations

Theta_chain = array(NA, c(S, K+1, p)) # each element is a (S, K+1) matrix collecting the S accepted cutoffs (rows) of variable j


# set initial values for the cutoffs

Theta_0 = matrix(NA, K+1, p)

for(j in 1:p){
  Theta_0[,j] = c(-Inf, sort(rtruncnorm(K-1, mean = 0, sd = 2, a = -Inf, b = Inf)), Inf)
}

Theta_0

Theta_chain[1,,] = Theta_0 # store the initial value Theta_0

Theta_chain[,1,]   = -Inf # set first and last cutoff equal to -Inf and Inf respectively
Theta_chain[,K+1,] = Inf

# (1) Propose the cutoffs

sigma_theta = 0.1 # sd of the proposal

# Recall that first (1) and last (K+1) cutoffs are fixed at -Inf and Inf
# Therefore there are (K-1) "free" cutoffs

s = 2 # iteration 2

j = 1 # first variable

theta_j_prop = c(-Inf,rep(NA,K-1),Inf)

for(k in 2:K){
  
  theta_j_prop[k] = rtruncnorm(1, mean = Theta_chain[s-1,k,j], sd = sigma_theta, a = theta_j_prop[k-1], b = Theta_chain[s-1,k+1,j])
  
}

theta_j_prop

# (2) Evaluate the proposal

theta_prop_eval = sum(log(dtruncnorm(theta_j_prop[2:K], mean = Theta_chain[s-1,2:K,j], sd = sigma_theta,
                                                        a    = Theta_chain[s-1,1:(K-1),j], b = theta_j_prop[3:(K+1)])))

theta_eval = sum(log(dtruncnorm(Theta_chain[s-1,2:K,j], mean = theta_j_prop[2:K], sd = sigma_theta,
                                                        a    = theta_j_prop[1:(K-1)], b = Theta_chain[s-1,3:(K+1),j])))

# log ratio for proposal evaluation

theta_eval - theta_prop_eval


# Likelihood evaluation

theta_x_j = match_x_theta_j(x = X[,j], theta = theta_j_prop)

match_x_theta_j_lower(x = X[,j], theta = theta_j_prop)

match_x_Theta_lower(x = X[i,], Theta = Theta_chain[s-1,,], p = p)


lower_X = t(apply(X = X, MARGIN = 1, FUN = match_x_Theta_lower, Theta = Theta_chain[s-1,,], p = p))
upper_X = t(apply(X = X, MARGIN = 1, FUN = match_x_Theta_upper, Theta = Theta_chain[s-1,,], p = p))

lower_X[1:10,]
upper_X[1:10,]


match_x_theta_j_lower(x = X[,j], theta = theta_j_prop)
match_x_theta_j_upper(x = X[,j], theta = theta_j_prop)

lower_X_tmp     = lower_X
lower_X_tmp[,j] = match_x_theta_j_lower(x = X[,j], theta = theta_j_prop)

upper_X_tmp     = upper_X
upper_X_tmp[,j] = match_x_theta_j_upper(x = X[,j], theta = theta_j_prop)


Sigma = diag(rep(1,p))

like_prop_eval = sum(log(sapply(1:n, function(i) pmvnorm(lower = lower_X_tmp[i,], upper = upper_X_tmp[i,], sigma = Sigma)))) # likelihood evaluated at the proposed cutoffs

like_eval = sum(log(sapply(1:n, function(i) pmvnorm(lower = lower_X[i,], upper = upper_X[i,], sigma = Sigma)))) # likelihood evaluated at the current cutoffs

like_prop_eval - like_eval
