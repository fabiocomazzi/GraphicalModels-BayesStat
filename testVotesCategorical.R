rm(list = ls())
setwd("~/GitHub/GraphicalModels-BayesStat")
source("categorical.R")
library(fields)
library(network)

# Load and split the dataset
data = read.table("votes.txt",sep=",",header=TRUE)
data[data == "y"] = 1
data[data == "n"] = 0
data[data == "?"] = NA
data = data[complete.cases(data),]
dataRepublican = data %>% filter(Party == "republican") %>% dplyr::select(-Party)
dataDemocrat = data %>% filter(Party == "democrat") %>% dplyr::select(-Party)

# Run the MCMC chains
initialCandidate = matrix(0,16,16)

## Republicans
chain = MetropolisHastingsCategorical(dataRepublican,initialCandidate,5000,1000,1,prior = "Binomial",p=0.3)
mpgRepublican = medianProbabilityGraph(chain)
mapRepublican = maximumPosterioriGraph(chain)
probabilitiesRepublican = chain[[1]]
for(i in 2:length(chain)){
  probabilitiesRepublican = probabilitiesRepublican + chain[[i]]
}
probabilitiesRepublican = probabilitiesRepublican / length(chain)

## Democrats
chain = MetropolisHastingsCategorical(dataDemocrat,initialCandidate,5000,1000,1,prior = "Binomial",p=0.3)
mpgDemocrat = medianProbabilityGraph(chain)
mapDemocrat = maximumPosterioriGraph(chain)
probabilitiesDemocrat = chain[[1]]
for(i in 2:length(chain)){
  probabilitiesDemocrat = probabilitiesDemocrat + chain[[i]]
}
probabilitiesDemocrat = probabilitiesDemocrat / length(chain)

# Plot the estimated graphs for the two groups
x11()
par(mfrow = c(1,2))
plotGraph(mpgRepublican,main = "Median Probability Graph")
plotGraph(mapRepublican,main = "Maximum a Posteriori")

x11()
par(mfrow = c(1,2))
plotGraph(mpgDemocrat,main = "Median Probability Graph")
plotGraph(mapDemocrat,main = "Maximum a Posteriori")

# Plot the heatmap of edge-inclusion probabilities
pdf("Heatmap_Republican.pdf", width = 6.5, height = 5.5)
colours = colorRampPalette(c('white','black'))
par(mar = c(4,4,1,2), oma = c(0.5,0.5,0.5,0.5), cex = 1, mgp = c(3,1,0))
image.plot(probabilitiesRepublican, col = colours(100), zlim = c(0,1), cex.sub = 1, axes = F, horizontal = F, legend.shrink = 1)
axis(1, at = seq(0, 1, l = 16), lab = 1:16, las = 2)
axis(2, at = seq(0, 1, l = 16), lab = 1:16, las = 1)
dev.off()

pdf("Heatmap_Democrat.pdf", width = 6.5, height = 5.5)
colours = colorRampPalette(c('white','black'))
par(mar = c(4,4,1,2), oma = c(0.5,0.5,0.5,0.5), cex = 1, mgp = c(3,1,0))
image.plot(probabilitiesDemocrat, col = colours(100), zlim = c(0,1), cex.sub = 1, axes = F, horizontal = F, legend.shrink = 1)
axis(1, at = seq(0, 1, l = 16), lab = 1:16, las = 2)
axis(2, at = seq(0, 1, l = 16), lab = 1:16, las = 1)
dev.off()

# Plot the MPG estimates in a circular shape for easier comparison
labs = as.character(c(1:16))
graphRepublican = network(mpgRepublican, label = labs)
graphDemocrat = network(mpgDemocrat, label = labs)
vertex_col = "gray90"

pdf("Votes_Graph_Republican.pdf", width = 6, height = 6)
par(mar = c(0,0,0,0), oma = c(0,0,0,0))
plot.network(graphRepublican, displaylabels = TRUE, vertex.col = vertex_col,
             mode = "circle",
             label.pos = 5,
             usecurve = TRUE, edge.curve = 0, vertex.cex = 2.5,
             label.cex = 0.8, edge.lwd = 0.1, arrowhead.cex = 0)
dev.off()

pdf("Votes_Graph_Democrat.pdf", width = 6, height = 6)
par(mar = c(0,0,0,0), oma = c(0,0,0,0))
plot.network(graphDemocrat, displaylabels = TRUE, vertex.col = vertex_col,
             mode = "circle",
             label.pos = 5,
             usecurve = TRUE, edge.curve = 0, vertex.cex = 2.5,
             label.cex = 0.8, edge.lwd = 0.1, arrowhead.cex = 0)
dev.off()