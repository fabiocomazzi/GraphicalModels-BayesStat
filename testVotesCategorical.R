rm(list = ls())
setwd("~/GitHub/GraphicalModels-BayesStat")
source("categorical.R")
data = read.table("votes.txt",sep=",",header=TRUE)
data[data == "y"] = 1
data[data == "n"] = 0
data[data == "?"] = NA
data = data[complete.cases(data),]
dataRepublican = data %>% filter(Party == "republican") %>% dplyr::select(-Party)
dataDemocrat = data %>% filter(Party == "democrat") %>% dplyr::select(-Party)

initialCandidate = matrix(0,16,16)

chain = MetropolisHastingsCategorical(dataRepublican,initialCandidate,5000,1000,1,prior = "Binomial",p=0.3)
mpgRepublican = medianProbabilityGraph(chain)
mapRepublican = maximumPosterioriGraph(chain)

chain = MetropolisHastingsCategorical(dataDemocrat,initialCandidate,5000,1000,1,prior = "Binomial",p=0.3)
mpgDemocrat = medianProbabilityGraph(chain)
mapDemocrat = maximumPosterioriGraph(chain)

x11()
par(mfrow = c(1,2))
plotGraph(mpgRepublican,main = "Median Probability Graph")
plotGraph(mapRepublican,main = "Maximum a Posteriori")

x11()
par(mfrow = c(1,2))
plotGraph(mpgDemocrat,main = "Median Probability Graph")
plotGraph(mapDemocrat,main = "Maximum a Posteriori")