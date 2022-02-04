setwd("~/GitHub/GraphicalModels-BayesStat")
source("utilityFunctions.R")

# data is a dataframe of categorical variables, variables.names is a vector of strings
# corresponding to the names of the variables to consider and variables.values is a 
# vector whose elements corresponds to the values of the variables in variables.names.
computeCounts = function(data, variables.values, variables.names = NULL){
  if(is.null(variables.names)){
    variables.names = colnames(data)
  }
  if(length(variables.names) != length(variables.values)){
    stop("The lenghts of variables.names and variables.values do not match")
  }
  temp = data
  for(i in 1:length(variables.names)){
    temp = temp %>% filter(temp[variables.names[i]] == variables.values[i])
  }
  return(dim(temp)[1])
}

# Computes the marginal likelihood restricted to a specific subset of nodes of the graph.
# a is a real positive number such that the parameters of the Dirichlet distribution are given
# by a/dim(X_s) where X_s is the cartesian space of all the modalities of the variables in
# subset.
logMarginalLikelihoodSubset = function(adjacencyMatrix,data,subset,a){
  if(a <= 0){
    stop("a should be a positive number!")
  }
  # The matrix table contains all the possible combinations of the modalities of the
  # variables in subset.
  subset = as.vector(subset)
  values = list()
  for(i in 1:length(subset)){
    values[[i]] = unique(data[subset[i]])[[1]]
  }
  table = as.matrix(expand.grid(values))
  # For each combination of the modalities of the variables (i.e., for each row in table)
  # we compute the number of times such combination is present in the dataframe. Such values
  # are stored in the vector count.
  count = c()
  for(i in 1:dim(table)[1]){
    variables.values = table[i,]
    count = c(count,computeCounts(data,variables.values,subset))
  }
  # This is the vector of the parameters of the prior (Dirichlet) distribution. 
  aVec = rep(a/dim(table)[1],dim(table)[1])
  # Finally, we compute the value of the marginal likelihood as the ratio of the normalizing constants
  # of the prior and posterior distribution.
  ml = lgamma(sum(aVec)) - lgamma(sum(aVec + count)) + sum(lgamma(aVec + count) - lgamma(aVec))
  return(ml)
}

# Computes the (total) marginal likelihood via the factorization property of decomposable graphs.
logMarginalLikelihood = function(adjacencyMatrix,data,a){
  if(a <= 0){
    stop("a should be a positive number!")
  }
  # We compute the cliques and separators of the graph.
  decomposition = getCliquesAndSeparators(adjacencyMatrix)
  cliques = decomposition[[1]]
  separators = decomposition[[2]]
  # For each clique we compute the associated marginal likelihood and store it in the
  # vector mlC.
  mlC = c()
  for(i in 1:length(cliques)){
    clique = cliques[[i]]
    mlC = c(mlC,logMarginalLikelihoodSubset(adjacencyMatrix,data,clique,a))
  }
  # For each separator we compute the associated marginal likelihood and store it in the
  # vector mlS.
  mlS = c()
  if(length(separators) == 0){
    mlS = c(0)
  }
  else{
    for(i in 1:length(separators)){
      separator = separators[[i]]
      mlS = c(mlS,logMarginalLikelihoodSubset(adjacencyMatrix,data,separator,a))
    }
  }
  # We finally compute the total marginal likelihood via the factorization property.
  result = sum(mlC) - sum(mlS)
  return(result)
}

# Implemenation of the Metropolis-Hastings algorithm to make inference on the 
# graph generating the data collected in the dataframe "data". n.iter is the length
# of the chain, thin is the thinning and burnin is the burnin.
MetropolisHastingsCategorical = function(data,initialCandidate,n.iter,burnin = 0,thin = 1,prior,p = NULL,a = NULL,b = NULL){
  # We check that the passed parameters are correct
  if(!prior %in% c("Uniform","Binomial","Beta-Binomial")){
    stop("prior should be either 'Uniform', 'Binomial' or 'Beta-Binomial'!")
  }
  if(!isDecomposable(initialCandidate)){
    stop("Initial candidate graph should be decomposable!")
  }
  currentCandidate = initialCandidate
  # Run the burnin iterations
  message("BURN-IN")
  progressBarBI = txtProgressBar(min = 0, max = burnin, initial = 0, style = 3) 
  for(i in 1:burnin){
    setTxtProgressBar(progressBarBI,i)
    newCandidate = newGraphProposal(currentCandidate)
    num = logMarginalLikelihood(newCandidate,data,2)
    den = logMarginalLikelihood(currentCandidate,data,2)
    marginalRatio = exp(num - den)
    priorRatio = switch(prior, "Uniform" = 1, "Binomial" = binomialPrior(currentCandidate,newCandidate,p), "Beta-Binomial" = betaBinomialPrior(currentCandidate,newCandidate,a,b))
    acceptanceProbability = min(marginalRatio * priorRatio,1)
    accepted = rbern(1,acceptanceProbability)
    if(accepted == 1){
      currentCandidate = newCandidate
    }
  }
  close(progressBarBI)
  # Run the chain
  message("Metropolis-Hastings")
  progressBar = txtProgressBar(min = 0, max = n.iter, initial = 0, style = 3) 
  chain = list()
  c = 0
  for(i in 1:n.iter){
    setTxtProgressBar(progressBar,i)
    newCandidate = newGraphProposal(currentCandidate)
    num = logMarginalLikelihood(newCandidate,data,2)
    den = logMarginalLikelihood(currentCandidate,data,2)
    marginalRatio = exp(num - den)
    priorRatio = switch(prior, "Uniform" = 1, "Binomial" = binomialPrior(currentCandidate,newCandidate,p), "Beta-Binomial" = betaBinomialPrior(currentCandidate,newCandidate,a,b))
    acceptanceProbability = min(marginalRatio * priorRatio,1)
    accepted = rbern(1,acceptanceProbability)
    c = c + accepted
    if(accepted == 1){
      currentCandidate = newCandidate
    }
    if(i %% thin == 0){
      chain[[i/thin]] = currentCandidate
    }
  }
  close(progressBarBI)
  cat(paste0("\nThe average acceptance rate is: ", as.character(c / n.iter),"\n"))
  
  return(chain)
}

# Dirichlet Process Mixture
DPMixture = function(data,n.iter,burnin,a_alpha,b_alpha,a_pi,b_pi){
  n = nrow(data)
  q = ncol(data)
  baseline = sampleFromBaseline(S = n.iter, burn = burnin, q = q, a_pi, b_pi)
  
  Xi_chain = matrix(NA, n, n.iter)
  # n x n.iter matrix which collects the cluster indicators of each unit for every iteration 
  
  A_chain = vector(mode = "list", length = n.iter)
  # List which collect the adjacency matrices of the visited graphs. In particular, each element of the
  # list is an array of K(t) adjacency matrices, where K(t) is the (random) number of mixture components
  # at iteration t.
  
  alpha_0_chain = rep(NA, n.iter)
  # Vector which collects the values of the precision parameter (alpha_0) for every iteration
  
  ## Set the initial values
  K_inits = 2 # Number of clusters
  alpha_0 = 1 # Precision parameter
  A_0 = array(0, c(q, q, K_inits)) # Graphs
  colnames(A_0) = rownames(A_0) = 1:q
  A_chain[[1]] = A_0
  xi = sample(K_inits, n, replace = TRUE) # Cluster indicators
  while(length(table(xi)) < K_inits){ # This loop makes sure that both clusters have at least one unit
    xi = sample(K_inits, n, replace = TRUE)
  }
  Xi_chain[,1] = xi
  
  
  graphs = A_0
  r = table(xi)
  
  ## MCMC iterations
  for(t in 2:n.iter){
    ### Update cluster indicators
    maxCl = length(r) # Maximum number of clusters
    ind = which(r != 0) # Indices of non-empty clusters
    r = r[ind]
    
    #### Update the parameters u and v
    if(length(r) == 1){
      v = rbeta(length(r), 1 + r, alpha_0)
    }
    else{
      v = rbeta(length(r), 1 + r, alpha_0 + c(rev(cumsum(rev(r)))[-1]))
    }
    v = c(v, rbeta(1, 1, alpha_0))

    omega_tmp = v[1]
    
    for(k in 2:(length(r)+1)){
      omega_tmp[k] = v[k]*prod(1 - v[1:(k-1)])
    }
    
    R = omegatmp[length(omega_tmp)] # R is the weight for a potential new cluster
    
    omega = numeric(maxCl)
    
    omega[ind] = omega_tmp[-length(omega_tmp)]
    
    u = stats::runif(n)*omega[xi]
    u_star = min(u)
    
    h = 0 # h is the number of non-empty clusters
    
    while(R > u_star){
      
      h = h+1
      beta_temp = stats::rbeta(n = 1, shape1 = 1, shape2 = alpha_0)
      
      omega = c(omega, R*beta_temp) # weight of the new cluster
      
      # probabilità che un individuo venga assegnato ad un nuovo cluster
      
      R = R * (1 - beta_temp) # remaining weight
      
      Dag_star = out_baseline[,,sample(n_base - burn_base, 1)]
      
      Dags = abind(Dags, Dag_star)
      
    }
    
    ## [2] ## Update of indicator variables xi
    
    K_star = dim(Dags)[3]
    probs = matrix(0,n,K_star)
    
    ### NEW CODE
    for(k in 1:K_star){
      for(i in 1:n){
        boolvec = rep(FALSE,n)
        boolvec[i] = TRUE
        r_i = sum(!boolvec & xi == k)
        L_i = length(unique(xi[!boolvec]))
        data_num = data[boolvec | xi == k,]
        data_den = data[!boolvec & xi == k,]
        if(k <= L_i){
          prob = r_i * exp(logMarginalLikelihood(Dags[,,k],data_num, a = 1) - logMarginalLikelihood(Dags[,,k],data_den,a = 1))
        }
        else{
          prob = alpha_0 * exp(logMarginalLikelihood(Dags[,,k],data_num, a = 1) - logMarginalLikelihood(Dags[,,k],data_den,a = 1))
        }
        probs[i,k] = prob
      }
    }
    
    probs = probs / rowSums(probs)
    
    xi_star = sapply(1:n, function(i) sample(1:(K_star), size = 1, prob = probs[i,]))
    
    labs = as.integer(names(table(xi_star)))
    
    K_star = length(labs)
    
    Dags  = array(Dags[,,labs], c(q, q, K_star))
    
    # Riassegno le etichette ai cluster in modo che vadano sempre da 1 a K
    
    xi_star = as.factor(xi_star); levels(xi_star) = 1:K_star   # update labels
    
    r = table(xi_star)
    
    omega = omega[labs]
    
    xi = c(xi_star)
    
    K = dim(Dags)[3]
    
    
    ###############################
    ## Update of alpha_0 given K ##
    ###############################
    
    eta = rbeta(1, alpha_0 + 1, n)
    
    alpha_0 = c(rgamma(1, shape = a_alpha + K, rate = b_alpha - log(eta)), rgamma(1, shape = a_alpha + K - 1, rate = b_alpha - log(eta)))[sample(c(1,2), 1, prob = c(a_alpha + K - 1, n*(b_alpha - log(eta))))]
    
    alpha_0_chain[t] = alpha_0
    
    
    ###############################
    ## Update DAGs D_1, ..., D_K ##
    ###############################
    
    set = 1:K
    
    for(k in set){ # per ciascun cluster
      currentCandidate = Dags[,,k]
      newCandidate = newGraphProposal(currentCandidate)
      num = logMarginalLikelihood(newCandidate,data,2)
      den = logMarginalLikelihood(currentCandidate,data,2)
      marginalRatio = exp(num - den)
      priorRatio = 1
      acceptanceProbability = min(marginalRatio * priorRatio,1)
      accepted = rbern(1,acceptanceProbability)
      if(accepted == 1){
        Dags[,,k] = newCandidate
      }
    }
    
    A_chain[[t]]     = Dags
    Xi_chain[,t]     = xi
    
    if(t%%100 == 0) print(paste0("Iteration ", t))
    
    
  }
  
  
  return(list())
  
}
