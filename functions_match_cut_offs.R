## Functions

## Dealing with observations and cutoffs

match_x_Theta = function(x, Theta, p){
  
  # return the lower and upper cutoff associated to each observation of the individual vector x
  
  return(cbind(Theta[cbind(x,1:p)], Theta[cbind(x+1,1:p)]))
  
}

match_x_Theta_lower = function(x, Theta, p){
  
  # return the lower cutoff associated to each observation of the individual vector x
  
  return(cbind(Theta[cbind(x,1:p)]))
  
}

match_x_Theta_upper = function(x, Theta, p){
  
  # return the upper cutoff associated to each observation of the individual vector x
  
  return(cbind(Theta[cbind(x+1,1:p)]))
  
}

# Or for a given variable j return the lower and upper cutoff for the all individuals in a (n,2) matrix

match_x_theta_j = function(x, theta){
  
  # x     : (n,1) vector of observations from variable j
  # theta : set of cutoffs associated to variable j
  
  # return a (n,2) matrix collecting the lower and upper cutoffs associated to each individual (row) relative to variable j
  
  return(cbind(theta[x], theta[x+1]))
  
}

match_x_theta_j_lower = function(x, theta){
  
  # x     : (n,1) vector of observations from variable j
  # theta : set of cutoffs associated to variable j
  
  # return a (n,1) matrix collecting the lower cutoffs associated to each individual (row) relative to variable j
  
  return(cbind(theta[x]))
  
}

match_x_theta_j_upper = function(x, theta){
  
  # x     : (n,1) vector of observations from variable j
  # theta : set of cutoffs associated to variable j
  
  # return a (n,1) matrix collecting the upper cutoffs associated to each individual (row) relative to variable j
  
  return(cbind(theta[x+1]))
  
}

from_latent_to_discrete = function(z, theta){
  
  # z : (n,1) vector of observations from the latent
  # theta : vector with K+1 thresholds with theta[1] = -Inf and theta[K+1] = Inf
  
  # return a discrete random variable with K levels
  
  K = length(theta) - 1
  
  sapply(1:n, function(i) max((1:K)[z[i] > theta]))
  
}
