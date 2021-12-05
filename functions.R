library(extraDistr)
library(dplyr)
library(igraph)
library(matrixcalc)

# n.obs is the number of observations to simulate (int), variables.names is a vector of 
# strings correpsonding to the names of the variables to simulate and variables.values
# is a list of vectors corresponding to the possible values of the variables (therefore
# the length of variables.values has to be the same of the length of variables.names).
generateCategoricalData = function(n.obs, variables.names, variables.values){
  data = data.frame(matrix(ncol = length(variables.names), nrow = n.obs))
  for(i in 1:length(variables.names)){
    values = variables.values[[i]]
    column = values[rdunif(n.obs,1,length(values))]
    data[,i] = column
  }
  colnames(data) = variables.names
  return(data)
}

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
  for(i in 1:length(variables.names)){
    data = data %>% filter(data[variables.names[i]] == variables.values[i])
  }
  return(dim(data)[1])
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
plotGraph = function(adjacencyMatrix, variables.names = NULL){
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
  plot(graph)
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
  cliques.list = list()
  for(i in 1:length(cliques)){
    clique = cliques[[i]]
    clique = colnames(adjacencyMatrix)[as.vector(clique)]
    cliques.list[[i]] = clique
  }
  separators.list = list()
  for(i in 1:length(separators)){
    separator = separators[[i]]
    separator = colnames(adjacencyMatrix)[as.vector(separator)]
    separators.list[[i]] = separator
  }
  return(list(cliques.list,separators.list))
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

