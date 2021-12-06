library(extraDistr)
library(dplyr)
library(igraph)
library(matrixcalc)
library(Rlab)
library(BDgraph)
library(mvtnorm)

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
  inv.covariance = rgwish(1, adj = adjacencyMatrix)
  covariance = solve(inv.covariance)
  mu = c(rep(0, n.variables))
  data = data.frame(rmvnorm(n.obs, mu, covariance))
  for(i in 1:n.variables){
    gamma = runif(1, quantile(data[,i], 0.05), quantile(data[,i], 0.95))
    data[,i][data[,i] >= gamma] = 1
    data[,i][data[,i] < gamma] = 0
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
  x11()
  plot(graph, main = main)
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
  graph = graph_from_adjacency_matrix(adjacencyMatrix, mode = "undirected")
  cliques = maximal.cliques(graph)
  separators = min_separators(graph)
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