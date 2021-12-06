setwd("~/GitHub/GraphicalMdels-BayesStat")
source("categorical.R")

data = generateCategoricalDataFromGraph(n.obs = 100, n.variables = 4, p = 0.6, variables.names = c("First","Second","Third","Fourth"))
plotGraph(data[[1]], main = "True Graph")
chain = MetropolisHastingsCategorical(data[[2]],100,1,100)
total = chain[[1]]
for(i in 2:length(chain)){
  total = total + chain[[i]]
}
total = total / length(chain)
total = replace(total, total < 0.5, 0)
total = replace(total, total >= 0.5, 1)
plotGraph(total)