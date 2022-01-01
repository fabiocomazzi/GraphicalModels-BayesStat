setwd("~/GitHub/GraphicalModels-BayesStat")
source("utilityFunctions.R")


generateGaussianDataFromGraph = function(adjacencyMatrix=NULL, n.obs, n.variables, p=NULL, covariance=NULL){
  #This function generate a dataset of n.obs number of observation issued from a multivariate gaussian with dimension n.variables
  #The data are being generated accoring to the graph given by adjacencyMatrix or by p which is the vector of probability that an edge will be in the graph
  #If the covariance Sigma is null, the function will sample one from the Hyper Inverse Wishart distribution.
  #The function returns a list with the adjacency matrix, the data and the covariance which generated the data.
  if(is.null(adjacencyMatrix)){
    while(TRUE){
      graph = erdos.renyi.game(n.variables,p,type="gnp",directed = FALSE)
      adjacencyMatrix = as_adjacency_matrix(graph, sparse = 0)
      if(isDecomposable(adjacencyMatrix)){
        break
      }
    }
  }
  if(!isDecomposable(adjacencyMatrix)){
    stop("Graph should be decomposable.")
  }
  
  if(is.null(covariance)){
    inv.covariance = rgwish(1, adj = adjacencyMatrix, D = 10 * diag(1,n.variables))
    covariance = solve(inv.covariance)
  }
  
  mu = c(rep(0, n.variables))
  data = dataCopy = data.frame(rmvnorm(n.obs, mu, covariance))
  return (list(adjacencyMatrix, data, covariance))
}

logMarginalLikelihoodGaussian = function(adjacencyMatrix, data, b, D){
  #this function compute the log marginal likelihood for data and for the graph represented by adjacencyMatrix
  n = dim(data)[1]
  q = dim(data)[2]
  b_star = b + n
  D_star = D + t(data)*data
  return (-log(2*pi)*n*q/2 + logh(adjacencyMatrix, b, D) - logh(adjacencyMatrix, b_star, D_star))
}

logh = function(adjacencyMatrix, b, D){
  #this function compute log(h) for the graph represented by adjacencyMatrix and parameter b and D
  decomposition = getCliquesAndSeparators(adjacencyMatrix)
  cliques = decomposition[[1]]
  separators = decomposition[[2]]
  # Loop for each clique 
  num = 0
  for(i in 1:length(cliques)){
    clique = cliques[[i]]
    Dc = D[clique,clique]
    cardC = len(clique)
    num = num + ((b + cardC - 1)/2)*log(det(Dc)/2) - lmvgamma((b + cardC - 1)/2, cradC) 
  }
  # Loop for each separator
  den = 0
  if(length(separators) == 0){
    message("Encauntered null set of separator")
    den = 1
  }
  else{
    for(i in 1:length(separators)){
      separator = separators[[i]]
      Ds = D[separator, separator]
      cardS = len(separator)
      den = den + ((b + cardS - 1)/2)*log(det(Ds)/2) - lmvgamma((b + cardS - 1)/2, cardS)
    }
  }
  # We finally compute the and return the results.
  result = num - den
  return(result)
}

MetropolisHastingsGaussian = function(data, initialCandidate, n.iter, burnin = 0, thin = 1, prior, p = NULL, a = NULL, b = NULL){
  # We check that the passed parameters are correct
  if(!prior %in% c("Uniform","Binomial","Beta-Binomial")){
    stop("prior should be either 'Uniform', 'Binomial' or 'Beta-Binomial'!")
  }
  if(!isDecomposable(initialCandidate)){
    stop("Initial candidate graph should be decomposable!")
  }
  currentCandidate = initialCandidate
  b = 1 #degress of freedom of the Hyperinverse Whishart
  D = 10 * diag(1, dim(data)[2]) #parameter D of the Hyperinverse Whishart
  # Run the burnin iterations
  message("BURN-IN")
  progressBarBI = txtProgressBar(min = 0, max = burnin, initial = 0, style = 3) 
  for(i in 1:burnin){
    setTxtProgressBar(progressBarBI,i)
    newCandidate = newGraphProposal(currentCandidate)
    num = logMarginalLikelihoodGaussian(newCandidate,data,b, D)
    den = logMarginalLikelihoodGaussian(currentCandidate,data,b, D)
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
    num = logMarginalLikelihoodGaussian(newCandidate,data,b, D)
    den = logMarginalLikelihoodGaussian(currentCandidate,data,b, D)
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


