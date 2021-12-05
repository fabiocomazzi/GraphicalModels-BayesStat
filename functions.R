library(extraDistr)
library(dplyr)
library(igraph)
library(matrixcalc)
library(Rlab)

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
  temp = data
  for(i in 1:length(variables.names)){
    temp = temp %>% filter(temp[variables.names[i]] == variables.values[i])
  }
  return(dim(temp)[1])
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

# Checks wheter the graph is decomposable or not. adjacencyMatrix
# is a matrix representing the adjacency matrix of an undirected graph.
# Optional: variables.names is a vector of strings representing
# the names of the nodes of the graph.
# isDecomposable = function(adjacencyMatrix,variables.names = NULL){
#   if(!isUndirectedGraph(adjacencyMatrix)){
#     stop("Adjacency matrix does not represent an undirected graph.")
#   }
#   if(!is.null(variables.names)){
#     if(length(variables.names) != dim(adjacencyMatrix)[1]){
#       stop("Length of variables.names is not correct.")
#     }
#     colnames(adjacencyMatrix) = rownames(adjacencyMatrix) = variables.names
#   }
#   graph = graph_from_adjacency_matrix(adjacencyMatrix, mode = "undirected")
#   # Find all the cycles of the graph
#   cycles = NULL
#   for(vertex1 in V(graph)) {
#     for(vertex2 in neighbors(graph, vertex1, mode="out")) {
#       cycles = c(cycles, lapply(all_simple_paths(graph, vertex2,vertex1, mode="out"), function(p) c(vertex1,p)))
#     }
#   }
#   # Find the cycles of length at least 4
#   longCycles = length(cycles[which(sapply(cycles, length) >= 5)])
#   if(longCycles > 0){
#     return(FALSE)
#   }
#   else{
#     return(TRUE)
#   }
# }
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

# Computes the marginal likelihood restricted to a specific subset of nodes of the graph.
# a is a real positive number such that the parameters of the Dirichlet distribution are given
# by a/dim(X_s) where X_s is the cartesian space of all the modalities of the variables in
# subset.
marginalLikelihoodSubset = function(adjacencyMatrix,data,subset,a){
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
  ml = (gamma(sum(aVec)) / gamma(sum(aVec + count))) * (prod(gamma(aVec + count) / gamma(aVec)))
  return(ml)
}

# Computes the (total) marginal likelihood via the factorization property of decomposable graphs.
marginalLikelihood = function(adjacencyMatrix,data,a){
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
    mlC = c(mlC,marginalLikelihoodSubset(adjacencyMatrix,data,clique,a))
  }
  # For each separator we compute the associated marginal likelihood and store it in the
  # vector mlS.
  mlS = c()
  if(length(separators) == 0){
    mlS = c(1)
  }
  else{
    for(i in 1:length(separators)){
      separator = separators[[i]]
      mlS = c(mlS,marginalLikelihoodSubset(adjacencyMatrix,data,separator,a))
    }
  }
  # We finally compute the total marginal likelihood via the factorization property.
  result = prod(mlC) / prod(mlS)
  return(result)
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

learnGraph = function(data,n.iter,thin,burnin){
  variables.names = colnames(data)
  # Generate a first candidate proposal. Observe that we make sure the generated
  # graph is decomposable
  while(TRUE){
    graph = erdos.renyi.game(length(variables.names),0.5,type="gnp",directed = FALSE)
    currentCandidate = as_adjacency_matrix(graph, sparse = 0)
    if(isDecomposable(currentCandidate)){
      break
    }
  }
  colnames(currentCandidate) = rownames(currentCandidate) = variables.names
  # Run the burnin iterations
  message("BURN-IN")
  progressBarBI = txtProgressBar(min = 0, max = burnin, initial = 0) 
  for(i in 1:burnin){
    setTxtProgressBar(progressBarBI,i)
    newCandidate = newGraphProposal(currentCandidate)
    num = marginalLikelihood(newCandidate,data,40)
    den = marginalLikelihood(currentCandidate,data,40)
    acceptanceProbability = min(num/den,1)
    accepted = rbern(1,acceptanceProbability)
    if(accepted == 1){
      currentCandidate = newCandidate
    }
  }
  close(progressBarBI)
  # Run the chain
  message("Metropolis-Hastings")
  progressBar = txtProgressBar(min = 0, max = n.iter, initial = 0) 
  chain = list()
  for(i in 1:n.iter){
    setTxtProgressBar(progressBar,i)
    newCandidate = newGraphProposal(currentCandidate)
    num = marginalLikelihood(newCandidate,data,40)
    den = marginalLikelihood(currentCandidate,data,40)
    acceptanceProbability = min(num/den,1)
    accepted = rbern(1,acceptanceProbability)
    if(accepted == 1){
      currentCandidate = newCandidate
    }
    if(i %% thin == 0){
      chain[[i/thin]] = currentCandidate
    }
  }
  close(progressBarBI)
  
  return(chain)
}


data = generateCategoricalData(100,c("Var1","Var2","Var3","Var4","Var5"),list(c(0,1),c("a","b","c"),c(1,2,3,4),c(TRUE,FALSE),c("Low","Medium","High")))
# adj = getAdjacencyMatrixFromEdges(c(1,2,2,3,3,4,4,1,4,2,1,5))
# colnames(adj) = rownames(adj) = c("Var1","Var2","Var3","Var4","Var5")
# clique = as.vector(getCliquesAndSeparators(adj)[[1]][[1]])
# marginalLikelihoodSubset(adj,data,clique,71.61)
# plotGraph(adj)
# cliques = getCliquesAndSeparators(adj)[[1]]
# clique = cliques[[1]]
# marginalLikelihoodSubset(adj,data,clique,2)
# marginalLikelihood(adj,data,1)
learnGraph(data,100,10,5)
