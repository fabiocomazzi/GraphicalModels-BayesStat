rm(list = ls())
setwd("~/GitHub/GraphicalModels-BayesStat")
source("gaussian_categorical.R")

# Generate 20 decomposable graph that will be used as the true graph to generate
# 20 different datasets
number_of_trial = 5
number_of_node = 6
trueGraphs = list()
encodedList = c()
for(i in 1:number_of_trial){
  while(TRUE){
    graph = erdos.renyi.game(number_of_node,0.3,type="gnp",directed = FALSE)
    newGraph = as_adjacency_matrix(graph, sparse = 0)
    encoded = encodeGraph(newGraph)
    if(isDecomposable(newGraph) & !encoded %in% encodedList){
      trueGraphs[[i]] = newGraph
      encodedList = c(encodedList,encoded)
      break
    }
  }
}

mpgs = list()
maps = list()
mpg_distances = c()
map_distances = c()
count = 1
for(trueGraph in trueGraphs){
  data = generateCategoricalDataFromGraph(adjacencyMatrix = trueGraph, n.obs = 1000, n.variables = number_of_node, p = 0.3)
  initialCandidate = matrix(0,number_of_node,number_of_node)
  print(paste("graph ", count))
  chain = MetropolisHastingsGaussianCategorical(data[[2]], initialCandidate, 100, 300, 1, prior = "Binomial", p=0.3)
  
  # Median Probability Graph
  mpg = medianProbabilityGraph(chain)
  mpgs[[count]] = mpg
  mpg_distances = c(mpg_distances,computeSHD(trueGraph,mpg))
  # Maximum a Posteriori Graph
  map = maximumPosterioriGraph(chain)
  maps[[count]] = map
  map_distances = c(map_distances,computeSHD(trueGraph,map))
  print(paste("map error: ", computeSHD(trueGraph, map)))
  print(paste("mpg error: ", computeSHD(trueGraph, mpg)))
  
  # Increase count
  count = count + 1
}

x11()
par(mfrow = c(2,1))
barplot(table(mpg_distances),main = "SHD (Median Probability Graphs)", col = rainbow(length(unique(mpg_distances))))
barplot(table(map_distances),main = "SHD (Maximum a Posteriori)", col = rainbow(length(unique(map_distances))))
print(mean(mpg_distances))
print(mean(map_distances))
