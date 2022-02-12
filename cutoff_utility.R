lower_bounds = function(tetha, x){
  lowb = matrix(rep(tetha, dim(x)[1]), nrow=dim(x)[1], byrow=TRUE)
  lowb[x==0] = -Inf
  return(lowb)
}

upper_bounds = function(tetha, x){
  lowb = matrix(rep(tetha, dim(x)[1]), nrow=dim(x)[1], byrow=TRUE)
  lowb[x==1] = Inf
  return(lowb)
}

generate_Z = function(Sigma, lower, upper, algorithm="gibbs"){
  #Sigma is the covariance matrix of the Gaussian
  #lower is a n*q matrix which stores in every row the lower bound of every z_i
  #upper is a n*q matrix which stores in every row the upper bound of every z_i
  #algorithm is the algorithm we use to sample every z_i 
  
  Z = matrix(NA, nrow=dim(lower)[1], ncol=dim(lower)[2])
  for(i in 1:dim(lower)[1]){
    Z[i,]=rtmvnorm(n=1, sigma=Sigma, lower=lower[i,], upper=upper[i,], algorithm=algorithm)
  }
  return(Z)
}

MH_tetha = function(Sigma, x, tau_prior, tetha, sd_proposal=2){
  #This function compute a step of metropolis harris to update the cutoffs vector
  
  for(j in 1:dim(x)[2]){
    prop_tetha_j = proposal_tetha_j(tetha[j], sd_proposal)
    temp_tetha = tetha
    temp_tetha[j] = prop_tetha_j
    
    acceptance = min(1, exp( log_density_tetha(temp_tetha, Sigma, x, tau_prior) - log_density_tetha(tetha, Sigma, x, tau_prior) 
                             #+ dproposal(tetha[j], prop_tetha_j, sd=sd_proposal, log=TRUE) 
                             #- dproposal(prop_tetha_j, tetha[j], sd=sd_proposal, log=TRUE)
                     ))
    
    if(rbern(1,acceptance)){
      tetha[j] = prop_tetha_j
    }
  }
  return(tetha)
}

proposal_tetha_j = function(tetha_j_old, sd=1){
  #This function sample from the proposal q(.|tetha_j_old)~N(tetha_j_old,sd)
  return(rnorm(n=1, mean=tetha_j_old, sd=sd))
}

dproposal = function(tetha_new_j, tetha_j, sd=1, log=FALSE){
  #the proposal q(tetha_new_j|tetha_j) ~ N(tetha_j, sd^2)
  #this function compute the density of the proposal
  
  return(dnorm(x=tetha_new_j, mean=tetha_j, sd=sd, log=log))
}

log_density_tetha = function(tetha, Sigma, x, tau_prior){
  # tetha is a q dimensional vector for which we want to evaluate the density
  # Sigma is the q*q covariance matrix
  # x is a n*q matrix with 0 and 1 corresponding to the categorical observations
  # tau_prior is the prior variance of the cutoffs
  
  lower_bounds_matrix = lower_bounds(tetha, x)
  upper_bounds_matrix = upper_bounds(tetha, x)
  
  log_dens = 0
  
  for(i in 1:dim(x)[1]){
    log_dens = log_dens + log(pmvnorm(lower=lower_bounds_matrix[i,], upper=upper_bounds_matrix[i,], sigma=Sigma)) 
  }
  
  for(j in 1:dim(x)[2]){
    log_dens = log_dens + dnorm(tetha[j], sd=tau_prior, log=TRUE)
  }
  
  return(log_dens)
}
