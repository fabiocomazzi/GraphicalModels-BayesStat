rm(list = ls())
setwd("~/GitHub/GraphicalMdels-BayesStat")
source("categorical.R")

data = generateCategoricalDataFromGraph(n.obs = 1000, n.variables = 6, p = 0.3)
plotGraph(data[[1]], main = "True Graph")

# Generate a first candidate proposal. Observe that we make sure the generated
# graph is decomposable
while(TRUE){
  graph = erdos.renyi.game(dim(data[[2]])[2],0.5,type="gnp",directed = FALSE)
  initialCandidate = as_adjacency_matrix(graph, sparse = 0)
  if(isDecomposable(initialCandidate)){
    break
  }
}
plotGraph(initialCandidate, main = "Initial Candidate")

# Run the chain
chain = MetropolisHastingsCategorical(data[[2]],initialCandidate,1000,1000,10,prior = "Binomial", p = 0.3)

# Median Probability Graph
mpg = medianProbabilityGraph(chain)
plotGraph(mpg,main = "Median Probability")

# Maximum a Posteriori Graph
map = maximumPosterioriGraph(chain)
plotGraph(map, main = "Maximum a Posteriori")
