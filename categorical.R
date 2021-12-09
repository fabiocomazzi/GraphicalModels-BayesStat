setwd("~/GitHub/GraphicalMdels-BayesStat")
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
MetropolisHastingsCategorical = function(data,initialCandidate,n.iter,burnin = 0,thin = 1,prior,p = NULL){
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
    priorRatio = switch(prior, "Uniform" = 1, "Binomial" = binomialPrior(currentCandidate,newCandidate,p), "Beta-Binomial" = 1)
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
    priorRatio = switch(prior, "Uniform" = 1, "Binomial" = binomialPrior(currentCandidate,newCandidate,p), "Beta-Binomial" = 1)
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