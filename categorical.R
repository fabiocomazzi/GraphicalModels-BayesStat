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

# Computes the contingency table of the observations in data.
countDuplicates = function(data){
  if(is.null(dim(data))){
    return(table(data))
  }
  
  x = do.call('paste', c(data, sep = '\r'))
  ordered_x = order(x)
  rl = rle(x[ordered_x])
  
  return(cbind(data[ordered_x[cumsum(rl$lengths)],,drop=FALSE], count = rl$lengths))
}

# Computes the marginal likelihood restricted to a specific subset of nodes of the graph.
# a is a real positive number such that the parameters of the Dirichlet distribution are given
# by a/dim(X_s) where X_s is the cartesian space of all the modalities of the variables in
# subset.
logMarginalLikelihoodSubset = function(data,subset,a){
  if(a <= 0){
    stop("a should be a positive number!")
  }
  n = dim(data)[1]
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
  ml = lgamma(a) - lgamma(a + n) + sum(lgamma(aVec + count) - lgamma(aVec))
  return(ml)
}

# Computes the (total) marginal likelihood via the factorization property of decomposable graphs.
logMarginalLikelihood = function(adjacencyMatrix,data,a){
  if(a <= 0){
    stop("a should be a positive number!")
  }
  if(dim(data)[1] == 0){
    return(0)
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
    mlC = c(mlC,logMarginalLikelihoodSubset(data,clique,a))
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
      mlS = c(mlS,logMarginalLikelihoodSubset(data,separator,a))
    }
  }
  # We finally compute the total marginal likelihood via the factorization property.
  result = sum(mlC) - sum(mlS)
  return(result)
}

# Computes the (total) marginal likelihood starting from the clique-separator decomposition of the graph.
logMarginalLikelihoodFromDecomposition = function(data,cliques,separators,a){
  if(a <= 0){
    stop("a should be a positive number!")
  }
  if(dim(data)[1] == 0){
    return(0)
  }
  # For each clique we compute the associated marginal likelihood and store it in the
  # vector mlC.
  mlC = c()
  for(i in 1:length(cliques)){
    clique = cliques[[i]]
    mlC = c(mlC,logMarginalLikelihoodSubset(data,clique,a))
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
      mlS = c(mlS,logMarginalLikelihoodSubset(data,separator,a))
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
  close(progressBar)
  cat(paste0("\nThe average acceptance rate is: ", as.character(c / n.iter),"\n"))
  
  return(chain)
}

# Computes the ratio of the marginal likelihoods restricted to the variables in subset
# and computed on two datasets which differ in just one unit.
logMarginalLikelihoodRatio = function(countsTable,n,modalities,index,subset,a){
  if(a <= 0){
    stop("a should be a positive number!")
  }
  subset = as.vector(subset)
  # We compute the following quantities:
  cardinality = prod(modalities[subset]) # The number of possible configurations of the variables in subset
  aX = a/cardinality # The prior parameter associated to the modality of the i-th observation
  count = 0
  if(index != 0){ # By construction if index == 0 then count = 0
    aggregatedCounts = aggregate(countsTable$count, by = as.list(countsTable[subset]), FUN = sum) # Compute the counts of the modalities restricted to the variables in subset
    modality = countsTable[index,-dim(countsTable)[2]] # Retrieve the modality of observed for the i-th unit
    modality = modality[subset] # Restrict the modality to the variables in subset
    m = which(apply(as.data.frame(aggregatedCounts[,-dim(aggregatedCounts)[2]]), 1, function(x) return(all(x == modality)))) # Get the row index of the count corresponding to the modality of the i-th unit
    count = aggregatedCounts[m,]$x # Get the count corresponding to the modality of the i-th unit
  }
  # Finally, we compute the value of the ratio of the marginal likelihoods
  ml = log((aX + count) / (a + n - 1))
  return(ml)
}

# Computes the (log) predictive distribution of unit in position given by iMask belonging in the cluster
# made by the units given by clusterMask.
logPredictiveDistribution = function(countsTable,n,modalities,index,cliques,separators,a){
  if(a <= 0){
    stop("a should be a positive number!")
  }
  # For each clique we compute the associated marginal likelihood and store it in the
  # vector mlC.
  mlC = c()
  for(i in 1:length(cliques)){
    clique = cliques[[i]]
    mlC = c(mlC,logMarginalLikelihoodRatio(countsTable,n,modalities,index,clique,a))
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
      mlS = c(mlS,logMarginalLikelihoodRatio(countsTable,n,modalities,index,separator,a))
    }
  }
  # We finally compute the total marginal likelihood via the factorization property.
  result = sum(mlC) - sum(mlS)
  return(result)
}

# Dirichlet Process Mixture
DPMixture = function(data,n.iter,burnin,a_alpha,b_alpha,a_pi,b_pi,a = 1){
  n = nrow(data)
  q = ncol(data)
  
  ## Draw a sample from the baseline measure
  baseline = sampleFromBaseline(S = n.iter, burn = burnin, q = q, a_pi, b_pi)
  
  ## Initialize the chains
  Xi_chain = matrix(NA, n, n.iter)
  # n x n.iter matrix which collects the cluster indicators of each unit for every iteration 
  A_chain = vector(mode = "list", length = n.iter)
  # List which collect the adjacency matrices of the visited graphs. In particular, each element of the
  # list is an array of K(t) adjacency matrices, where K(t) is the (random) number of mixture components
  # at iteration t.
  alpha_0_chain = rep(NA, n.iter)
  # Vector which collects the values of the precision parameter (alpha_0) for every iteration
  
  ## Set the initial values
  K = 2 # Number of clusters
  alpha_0 = 1 # Precision parameter
  A = array(0, c(q, q, K)) # Group-specific graphs
  colnames(A) = rownames(A) = 1:q
  A_chain[[1]] = A
  xi = sample(K, n, replace = TRUE) # Cluster indicators
  while(length(table(xi)) < K){ # This loop makes sure that both clusters have at least one unit
    xi = sample(K, n, replace = TRUE)
  }
  Xi_chain[,1] = xi
  graphs = A
  r = table(xi)
  
  ## MCMC iterations
  message("RUNNING THE CHAIN")
  progressBar = txtProgressBar(min = 2, max = n.iter, initial = 2, style = 3)
  for(t in 2:n.iter){
    setTxtProgressBar(progressBar,t)
    
    ## Update of indicator variables xi
    maxCl = length(r) + 1 # The new (potential) number of cluster is given by the old number of cluster + 1
    graphs = abind(graphs, baseline[,,sample(n.iter - burnin, 1)]) # Sample a new graph from the baseline
    probs = matrix(0,n,maxCl) # The unit i is assigned to cluster j with probability = probs[i,j]
    for(k in 1:(maxCl-1)){ # For every cluster
      graph = graphs[,,k] # Group-specific graph
      # Compute the cliques and separators of the group-specific graph
      decomposition = getCliquesAndSeparators(graph)
      cliques = decomposition[[1]]
      separators = decomposition[[2]]
      # Compute the marginal likelihood with respect to all observations in cluster k
      data_clust = data[xi == k,]
      ml_clust = logMarginalLikelihoodFromDecomposition(data_clust,cliques,separators,a)
      for(i in 1:n){ # For every unit
        L_i = length(unique(xi[-i])) # Number of clusters in the sample excluding observation x_i
        if(L_i == (maxCl - 1) | xi[i] != k){ # If the unit is the only one in the cluster k then it cannot be reassigned to k (the corresponding probability is not updated and remains 0)
          iMask = rep(FALSE,n)
          iMask[i] = TRUE # iMask is TRUE in the position corresponding to the i-th unit, FALSE otherwise
          r_ki = sum(!iMask & xi == k) # Number of observations included in cluster k excluding observation x_i
          if(xi[i] == k){
            ml_num = ml_clust # Marginal likelihood of the numerator
            data_den = data[!iMask & xi == k,] # Dataset containing all the observation in cluster k excluding x_i
            ml_den = logMarginalLikelihoodFromDecomposition(data_den,cliques,separators,a) # Marginal likelihood of the denominator
          }
          else{
            ml_den = ml_clust # Marginal likelihood of the denominator
            data_num = data[iMask | xi == k,] # Dataset containing observation x_i plus all the observation in cluster k
            ml_num = logMarginalLikelihoodFromDecomposition(data_num,cliques,separators,a) # Marginal likelihood of the numerator
          }
          prob = r_ki * exp(ml_num - ml_den)
          probs[i,k] = prob
        }
      }
    }
    # Compute the probability that the unit is assigned to the new cluster (with index maxCl)
    graph = graphs[,,maxCl]
    # Compute the cliques and separators of the group-specific graph
    decomposition = getCliquesAndSeparators(graph)
    cliques = decomposition[[1]]
    separators = decomposition[[2]]
    for(i in 1:n){
      ml = logMarginalLikelihoodFromDecomposition(data[i,],cliques,separators,a = 1)
      prob = alpha_0 * exp(ml)
      probs[i,maxCl] = prob
    }
    probs = probs / rowSums(probs) # Normalize the matrix of probabilities
    xiNew = sapply(1:n, function(i) sample(1:(maxCl), size = 1, prob = probs[i,])) # New cluster allocation indices
    labels = as.integer(names(table(xiNew)))
    K = length(labels) # New number of clusters
    graphs = array(graphs[,,labels], c(q, q, K)) # Cluster-specific graphs
    ### Reassign the labels so that they span from 1 to K
    xiNew = as.factor(xiNew)
    levels(xiNew) = 1:K
    r = table(xiNew)
    xi = c(xiNew)
    
    ## Update of alpha_0
    eta = rbeta(1, alpha_0 + 1, n)
    alpha_0 = c(rgamma(1, shape = a_alpha + K, rate = b_alpha - log(eta)),
                rgamma(1, shape = a_alpha + K - 1, rate = b_alpha - log(eta)))[sample(c(1,2), 1, prob = c(a_alpha + K - 1, n*(b_alpha - log(eta))))]
    alpha_0_chain[t] = alpha_0
    
    ## Update of cluster-specific graphs
    for(k in 1:K){
      currentCandidate = graphs[,,k]
      newCandidate = newGraphProposal(currentCandidate)
      num = logMarginalLikelihood(newCandidate,data,2)
      den = logMarginalLikelihood(currentCandidate,data,2)
      marginalRatio = exp(num - den)
      priorRatio = 1
      acceptanceProbability = min(marginalRatio * priorRatio,1)
      accepted = rbern(1,acceptanceProbability)
      if(accepted == 1){
        graphs[,,k] = newCandidate
      }
    }
    
    ## Update the chain status
    A_chain[[t]] = graphs
    Xi_chain[,t] = xi
  }
  close(progressBar)
  
  ## Construct the similarity matrix
  similarityMatrix = matrix(0, nrow = n, ncol = n)
  for(t in (burnin + 1):n.iter){
    similarityMatrix = similarityMatrix + as.integer(matrix(Xi_chain[,t], nrow = n, ncol = n) == t(matrix(Xi_chain[,t], nrow = n, ncol = n)))
  }
  similarityMatrix = similarityMatrix / (n.iter - burnin)
  
  return(list(similarityMatrix = similarityMatrix))
}

# Dirichlet Process Mixture (efficient implementation)
DPMixture_Efficient = function(data,n.iter,burnin,a_alpha,b_alpha,a_pi,b_pi,a = 1){
  n = nrow(data)
  q = ncol(data)
  # Compute the total number of possible configurations of the variables in the dataset
  modalities = c()
  for(j in 1:dim(data)[2]){
    modalities = c(modalities, length(unique(data[,j])))
  }
  
  ## Draw a sample from the baseline measure
  baseline = sampleFromBaseline(S = n.iter, burn = burnin, q = q, a_pi, b_pi)
  
  ## Initialize the chains
  Xi_chain = matrix(NA, n, n.iter)
  # n x n.iter matrix which collects the cluster indicators of each unit for every iteration 
  A_chain = vector(mode = "list", length = n.iter)
  # List which collect the adjacency matrices of the visited graphs. In particular, each element of the
  # list is an array of K(t) adjacency matrices, where K(t) is the (random) number of mixture components
  # at iteration t.
  alpha_0_chain = rep(NA, n.iter)
  # Vector which collects the values of the precision parameter (alpha_0) for every iteration
  
  ## Set the initial values
  K = 2 # Number of clusters
  alpha_0 = 1 # Precision parameter
  A = array(0, c(q, q, K)) # Group-specific graphs
  colnames(A) = rownames(A) = 1:q
  A_chain[[1]] = A
  xi = sample(K, n, replace = TRUE) # Cluster indicators
  while(length(table(xi)) < K){ # This loop makes sure that both clusters have at least one unit
    xi = sample(K, n, replace = TRUE)
  }
  Xi_chain[,1] = xi
  graphs = A
  r = table(xi)
  
  ## MCMC iterations
  message("RUNNING THE CHAIN")
  progressBar = txtProgressBar(min = 2, max = n.iter, initial = 2, style = 3)
  for(t in 2:n.iter){
    setTxtProgressBar(progressBar,t)
    
    ## Update of indicator variables xi
    maxCl = length(r) + 1 # The new (potential) number of cluster is given by the old number of cluster + 1
    graphs = abind(graphs, baseline[,,sample(n.iter - burnin, 1)]) # Sample a new graph from the baseline
    probs = matrix(0,n,maxCl) # The unit i is assigned to cluster j with probability = probs[i,j]
    for(k in 1:(maxCl-1)){ # For every cluster
      clusterMask = xi == k
      countsTableTotal = countDuplicates(data[clusterMask,]) # Table of counts for the units in cluster k
      graph = graphs[,,k] # Group-specific graph
      # Compute the cliques and separators of the group-specific graph
      decomposition = getCliquesAndSeparators(graph)
      cliques = decomposition[[1]]
      separators = decomposition[[2]]
      for(i in 1:n){ # For every unit
        L_i = length(unique(xi[-i])) # Number of clusters in the sample excluding observation x_i
        countsTable = countsTableTotal
        index = which(apply(countsTable[,-dim(countsTable)[2]], 1, function(x) return(all(x == data[i,])))) # Row of countsTable corresponding to the configuration of the i-th unit
        if(length(index) == 0){ # Set index = 0 if the configuration of the i-th unit does not appear in the cluster k
          index = 0
        }
        if(xi[i] == k){ # If unit i belongs to cluster k than we need to reduce the count of the configuration corresponding to i by 1
          countsTable[index,]$count = countsTable[index,]$count - 1
        }
        if(L_i == (maxCl - 1) | xi[i] != k){ # If the unit is the only one in the cluster k then it cannot be reassigned to k (the corresponding probability is not updated and remains 0)
          iMask = rep(FALSE,n)
          iMask[i] = TRUE # iMask is TRUE in the position corresponding to the i-th unit, FALSE otherwise
          r_ki = sum(!iMask & clusterMask) # Number of observations included in cluster k excluding observation x_i
          n_obs = sum(clusterMask | iMask) # The number of observation in the considered cluster (plus eventually the i-th observation)
          prob = r_ki * exp(logPredictiveDistribution(countsTable,n = n_obs,modalities,index,cliques,separators,a))
          probs[i,k] = prob
        }
      }
    }
    # Compute the probability that the unit is assigned to the new cluster (with index maxCl)
    graph = graphs[,,maxCl]
    # Compute the cliques and separators of the group-specific graph and the predictive of the unit belonging to a new cluster
    decomposition = getCliquesAndSeparators(graph)
    cliques = decomposition[[1]]
    separators = decomposition[[2]]
    num = 1
    den = 1
    for(clique in cliques){
      num = num / prod(modalities[clique])
    }
    for(separator in separators){
      den = den / prod(modalities[separator])
    }
    predictive = num / den
    prob = alpha_0 * predictive
    for(i in 1:n){
      probs[i,maxCl] = prob
    }
    probs = probs / rowSums(probs) # Normalize the matrix of probabilities
    xiNew = sapply(1:n, function(i) sample(1:(maxCl), size = 1, prob = probs[i,])) # New cluster allocation indices
    labels = as.integer(names(table(xiNew)))
    K = length(labels) # New number of clusters
    graphs = array(graphs[,,labels], c(q, q, K)) # Cluster-specific graphs
    ### Reassign the labels so that they span from 1 to K
    xiNew = as.factor(xiNew)
    levels(xiNew) = 1:K
    r = table(xiNew)
    xi = c(xiNew)
    
    ## Update of alpha_0
    eta = rbeta(1, alpha_0 + 1, n)
    alpha_0 = c(rgamma(1, shape = a_alpha + K, rate = b_alpha - log(eta)),
                rgamma(1, shape = a_alpha + K - 1, rate = b_alpha - log(eta)))[sample(c(1,2), 1, prob = c(a_alpha + K - 1, n*(b_alpha - log(eta))))]
    alpha_0_chain[t] = alpha_0
    
    ## Update of cluster-specific graphs
    for(k in 1:K){
      currentCandidate = graphs[,,k]
      newCandidate = newGraphProposal(currentCandidate)
      num = logMarginalLikelihood(newCandidate,data,2)
      den = logMarginalLikelihood(currentCandidate,data,2)
      marginalRatio = exp(num - den)
      priorRatio = 1
      acceptanceProbability = min(marginalRatio * priorRatio,1)
      accepted = rbern(1,acceptanceProbability)
      if(accepted == 1){
        graphs[,,k] = newCandidate
      }
    }
    
    ## Update the chain status
    A_chain[[t]] = graphs
    Xi_chain[,t] = xi
  }
  close(progressBar)
  
  ## Construct the similarity matrix
  similarityMatrix = matrix(0, nrow = n, ncol = n)
  for(t in (burnin + 1):n.iter){
    similarityMatrix = similarityMatrix + as.integer(matrix(Xi_chain[,t], nrow = n, ncol = n) == t(matrix(Xi_chain[,t], nrow = n, ncol = n)))
  }
  similarityMatrix = similarityMatrix / (n.iter - burnin)
  
  return(list(similarityMatrix = similarityMatrix))
}

