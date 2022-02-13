rm(list = ls())
setwd("~/GitHub/GraphicalModels-BayesStat")
source("categorical.R")
library(mcclust)

generateDecomposableGraph = function(q,p){
  while(TRUE){
    graph = erdos.renyi.game(q,p,type="gnp",directed = FALSE)
    graph = as_adjacency_matrix(graph, sparse = 0)
    if(isDecomposable(graph)){
      return(graph)
    }
  }
}

similarityMatrixtoClustering = function(similarityMatrix){
  similarityMatrix = round(similarityMatrix)
  indices = c()
  for(i in n:1){
    indices[similarityMatrix[i,] == 1] = i
  }
  indices = as.factor(indices)
  levels(indices) = 1:(length(levels(indices)))
  return(indices)
}

# Generate two different graphs
q = 10 # Number of variables
n_1 = 100 # Sample size for group 1
n_2 = 100 # Sample size for group 2
p = 2 / q # Edge-inclusion probability
graph1 = generateDecomposableGraph(q,p)
graph2 = generateDecomposableGraph(q,p)
x11()
plotGraph(graph1)
x11()
plotGraph(graph2)

# Generate the data
X_1 = generateCategoricalDataFromGraph(adjacencyMatrix = graph1, n.obs = n_1, n.variables = q)[[2]]
X_2 = generateCategoricalDataFromGraph(adjacencyMatrix = graph2, n.obs = n_2, n.variables = q)[[2]]
Xi = c(rep(1, n_1), rep(2, n_2)) # True clusters indicators
X = rbind(X_1, X_2) # Dataset
n = n_1 + n_2 # Total sample size

# Run the MCMC
n.iter = 2500
burnin = 500
a_pi = 1
b_pi = 2*(q-2)/3
a_alpha = 1
b_alpha = 3
# mcmc = DPMixture(data = X, n.iter, burnin, a_alpha, b_alpha, a_pi, b_pi)
mcmc = DPMixture_Efficient(data = X, n.iter, burnin, a_alpha, b_alpha, a_pi, b_pi)


# Cluster estimates
similarityMatrix = mcmc$similarityMatrix
clustering = similarityMatrixtoClustering(similarityMatrix)
# Variation of Information (VI) between true and estimated clustering
vi.dist(Xi,clustering)