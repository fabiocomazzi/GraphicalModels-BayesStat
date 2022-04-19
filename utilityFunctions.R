library(extraDistr)
library(dplyr)
library(igraph)
library(matrixcalc)
library(Rlab)
library(BDgraph)
library(mvtnorm)
library(plyr)
library(pcalg)
library(gRbase)
library(abind)
library(fields)
library(network)

# n.obs is the number of observations to simulate (int), n.variables is the number of 
# random variables to generate and (optional) variables.names is a vector of 
# strings correpsonding to the names of the variables to simulate. The function takes as input
# either the adjacency matrix of an undirected decomposable graph or double p representing the
# probability of connecting two nodes of the graph in order to generate a random decomposable graph.
generateCategoricalDataFromGraph = function(adjacencyMatrix = NULL, n.obs, n.variables, p = NULL, variables.names = NULL){
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
  inv.covariance = rgwish(1, adj = adjacencyMatrix, D = 10 * diag(1,n.variables))
  covariance = solve(inv.covariance)
  mu = c(rep(0, n.variables))
  data = dataCopy = data.frame(rmvnorm(n.obs, mu, covariance))
  for(i in 1:n.variables){
    while(TRUE){
      gamma = runif(1, quantile(dataCopy[,i], 0.2), quantile(dataCopy[,i], 0.8))
      data[,i][dataCopy[,i] >= gamma] = 1
      data[,i][dataCopy[,i] < gamma] = 0
      if(length(unique(data[,i])) > 1){
        break
      }
    }
  }
  
  if(!is.null(variables.names)){
    if(length(variables.names) != dim(data)[2]){
      stop("Dimension of variables.names does not match.")
    }
    colnames(adjacencyMatrix) = rownames(adjacencyMatrix) = variables.names
    colnames(data) =  variables.names
  }
  return(list(adjacencyMatrix, data))
}

# Checks if adjacencyMatrix is the adjacency matrix of an undirected graph. 
isUndirectedGraph = function(adjacencyMatrix){
  if(!is.square.matrix(adjacencyMatrix)){
    message("Matrix is not the adjacency matrix of a graph as it is not square.")
    return(FALSE)
  }
  if(!isSymmetric(adjacencyMatrix)){
    message("Graph is not undirected as the adjacency matrix is not symmetric.")
    return(FALSE)
  }
  return(TRUE)
}

# Plots the graph given its adjacency matrix. Optional: variables.names is a vector
# of strings representing the names of the nodes of the graph.
plotGraph = function(adjacencyMatrix, variables.names = NULL, main = NULL){
  if(!isUndirectedGraph(adjacencyMatrix)){
    stop("Adjacency matrix does not represent an undirected graph.")
  }
  if(!is.null(variables.names)){
    if(length(variables.names) != dim(adjacencyMatrix)[1]){
      stop("Length of variables.names is not correct.")
    }
    colnames(adjacencyMatrix) = rownames(adjacencyMatrix) = variables.names
  }
  graph = graph_from_adjacency_matrix(adjacencyMatrix, mode = "undirected")
  plot(graph, main = main)
}

# Plots the graph in a circular shape given its adjacency matrix.
plotGraphCircular = function(adjacencyMatrix){
  labs = as.character(c(1:dim(adjacencyMatrix)[1]))
  graph = network(adjacencyMatrix, label = labs)
  vertex_col = "gray90"
  plot.network(graph, displaylabels = TRUE, vertex.col = vertex_col,
               mode = "circle",
               label.pos = 5,
               usecurve = TRUE, edge.curve = 0, vertex.cex = 2.5,
               label.cex = 0.8, edge.lwd = 0.1, arrowhead.cex = 0)
}

# Checks wheter the graph is decomposable or not. adjacencyMatrix
# is a matrix representing the adjacency matrix of an undirected graph.
# Optional: variables.names is a vector of strings representing
# the names of the nodes of the graph.
isDecomposable = function(adjacencyMatrix,variables.names = NULL){
    if(!isUndirectedGraph(adjacencyMatrix)){
      stop("Adjacency matrix does not represent an undirected graph.")
    }
    if(!is.null(variables.names)){
      if(length(variables.names) != dim(adjacencyMatrix)[1]){
        stop("Length of variables.names is not correct.")
      }
      colnames(adjacencyMatrix) = rownames(adjacencyMatrix) = variables.names
    }
    graph = graph_from_adjacency_matrix(adjacencyMatrix, mode = "undirected")
    return(is_chordal(graph)$chordal)
}

# Computes the (maximal) cliques and (minimal) separators of the graph. adjacencyMatrix
# is a matrix representing the adjacency matrix of an undirected graph. value is a list
# of two lists: the first one is a list of the cliques and the second one is a list of
# the separators. Optional: variables.names is a vector
# of strings representing the names of the nodes of the graph.
getCliquesAndSeparators = function(adjacencyMatrix,variables.names = NULL){
  if(!isUndirectedGraph(adjacencyMatrix)){
    stop("Adjacency matrix does not represent an undirected graph.")
  }
  if(!is.null(variables.names)){
    if(length(variables.names) != dim(adjacencyMatrix)[1]){
      stop("Length of variables.names is not correct.")
    }
    colnames(adjacencyMatrix) = rownames(adjacencyMatrix) = variables.names
  }
  decomposition = mpd(adjacencyMatrix)
  cliques = as.list(lapply(decomposition$cliques, as.numeric))
  cliques = unname(cliques[sapply(cliques, length) > 0])
  separators = as.list(lapply(decomposition$separators, as.numeric))
  separators = unname(separators[sapply(separators, length) > 0])
  # graph = graph_from_adjacency_matrix(adjacencyMatrix, mode = "undirected")
  # cliques = maximal.cliques(graph)
  # separators = min_separators(graph)
  return(list(cliques,separators))
}

# Computes the adjacency matrix of an undirected graph given a vector representing
# its edges. In particualar, edges is a vector such that the first edge is between
# the first and second element of edges, the second edge is between the third and the
# fourth element of edges and so on. Optional: variables.names is a vector
# of strings representing the names of the nodes of the graph.
getAdjacencyMatrixFromEdges = function(edges, variables.names = NULL){
  graph = make_graph(edges, directed = FALSE)
  adjacencyMatrix = as_adjacency_matrix(graph, sparse = 0)
  if(!is.null(variables.names)){
    if(length(variables.names) != dim(adjacencyMatrix)[1]){
      stop("Length of variables.names is not correct.")
    }
    colnames(adjacencyMatrix) = rownames(adjacencyMatrix) = variables.names
  }
  return(adjacencyMatrix)
}

# Given the adjacency matrix of a decomposable graph, computes the adjacency matrix
# of a new decomposable graph obtained by either adding or removing a single edge
# from the original graph. Observe that the proposal density is choosen so that
# q(G|G')/q(G'|G) = 1.
newGraphProposal = function(adjacencyMatrix){
  if(!isDecomposable(adjacencyMatrix)){
    stop("The graph must be decomposable!")
  }
  proposals = list()
  count = 1
  for(i in 1:(dim(adjacencyMatrix)[1] - 1)){
    for(j in 1:(dim(adjacencyMatrix)[2] - i)){
      newAdjacencyMatrix = adjacencyMatrix
      value = newAdjacencyMatrix[i,i+j]
      newAdjacencyMatrix[i,i+j] = newAdjacencyMatrix[i+j,i] = as.integer(!value)
      proposals[[count]] = newAdjacencyMatrix
      count = count + 1
    }
  }
  while(TRUE){
    proposalIndex = rdunif(1,1,count-1)
    if(isDecomposable(proposals[[proposalIndex]])){
      return(proposals[[proposalIndex]])
    }
  }
}

# Computes the prior ratio for the MH algorithm in the case of a Binomial prior
binomialPrior = function(currentProposal,newProposal,p){
  if(is.null(p)){
    stop("p should be between 0 and 1!")
  }
  if(p < 0 | p > 1){
    stop("p should be between 0 and 1!")
  }
  currentEdges = sum(currentProposal) / 2
  newEdges = sum(newProposal) / 2
  q = dim(currentProposal)[1]
  ratio = ((p ** newEdges) * ((1 - p) ** (q*(q-1)/2 - newEdges))) / ((p ** currentEdges) * ((1 - p) ** (q*(q-1)/2 - currentEdges)))
  return(ratio)
}

betaBinomialPrior = function(currentProposal,newProposal,a,b){
  num = logIntegralBetaBinomial(newProposal,a,b)
  den = logIntegralBetaBinomial(currentProposal,a,b)
  return(exp(num-den))
}

logIntegralBetaBinomial = function(graph,a,b){
  if(is.null(a) | is.null(b)){
    stop("a and b should be positive real numbers!")
  }
  if(a<=0 | b<=0){
    stop("a and b should be positive real numbers!")
  }
  edges = sum(graph) / 2
  q = dim(graph)[1]
  integral = lgamma(a + edges) + lgamma(b + q*(q-1)/2 - edges) - lgamma(q*(q-1)/2 + a + b) + lgamma(a + b) - lgamma(a) - lgamma(b)
  return(integral)
}

# Computes the graph that includes all the edges with posterior inclusion probability >= 0.5
medianProbabilityGraph = function(chain){
  mpg = chain[[1]]
  for(i in 2:length(chain)){
    mpg = mpg + chain[[i]]
  }
  mpg = mpg / length(chain)
  mpg = replace(mpg, mpg < 0.5, 0)
  mpg = replace(mpg, mpg >= 0.5, 1)
  
  return(mpg)
}

# Computes the most frequent graph
maximumPosterioriGraph = function(chain){
  uniqueValues = count(unlist(lapply(chain, toString)))
  map.index = which(uniqueValues$freq == max(uniqueValues$freq))
  map = uniqueValues$x[map.index]
  map = as.integer(unlist(strsplit(map,", ")))
  map = matrix(map,nrow = sqrt(length(map)))
  
  return(map)
}

# Computes the Structural Hamming Distance between two graphs given their adjacency matrices
computeSHD = function(adjacencyMatrix1, adjacencyMatrix2){
  graph1 = as_graphnel(graph_from_adjacency_matrix(adjacencyMatrix1,mode = "undirected"))
  graph2 = as_graphnel(graph_from_adjacency_matrix(adjacencyMatrix2,mode = "undirected"))
  return(shd(graph1,graph2))
}

# Encondes a graph with a unique integer value given its adjacency matrix
encodeGraph = function(adjacencyMatrix){
  result = ""
  for(i in 1:(dim(adjacencyMatrix)[1] - 1)){
    for(j in 1:(dim(adjacencyMatrix)[2] - i)){
      result = paste0(result,as.character(adjacencyMatrix[i,i+j]))
    }
  }
  result = strtoi(result,base = 2)
  return(result)
}

# Creates a sample drawn from the baseline measure over graphs. In particular, S is the number
# of desired draws, burn is the burn-in, q is the number of nodes in the graph and a_pi, b_pi are
# the hyper-parameters of the Beta prior on probability of edge inclusion pi.
sampleFromBaseline = function(S, burn, q, a_pi, b_pi){
  chain = array(NA, c(q, q, S)) # Initialize the chain with S adjacency matrices of NA's
  # The initial graph has no edges
  graph = matrix(0, q, q)
  chain[,,1] = graph
  # MCMC sampler to draw the remaining graphs from the baseline measure
  for(s in 2:S){
    newGraph = newGraphProposal(graph) # Draw a decomposable graph from the neighborhood of the current graph
    # Compute the multiplicity correction (log)prior.
    # Note that the number of edges of the graph is given by sum(adjacencyMatrix) / 2 since the graphs are decomposable
    logPriorNew = lgamma((sum(newGraph) / 2) + a_pi) + lgamma(q*(q-1)/2 - (sum(newGraph) / 2) + b_pi - 1) # New candidate
    logPriorOld = lgamma((sum(graph) / 2) + a_pi) + lgamma(q*(q-1)/2 - (sum(graph) / 2) + b_pi - 1) # Current candidate
    logPrior = logPriorNew - logPriorOld
    # Acceptance ratio
    acceptanceRatio = min(0, logPrior)
    # Check if the new candidate is accepted
    if(log(runif(1)) < acceptanceRatio){
      graph = newGraph
    }
    chain[,,s] = graph
  }
  chain = chain[,,(burn + 1):S] # Discard burn-in samples
  return(chain)
}