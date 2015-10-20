
## BRAF Network exploration

### Load in datasets
mut_data <- read.csv("/Users/aolow/Box\ Sync/Documents/VantveerLab/MS1_SNPs/120214_mut_in_pepts_cosmic_corrected.csv")
nodes <- read.csv("/Users/aolow/Dropbox/MS1_Rcode/NEW_NODES3.csv", as.is=T)
#n_sub <- read.csv("/Users/aolow/Dropbox/MS1_Rcode/100214_nodes_substrates.csv")
nodes <- nodes[,1:3]

nodes <- nodes[!is.na(nodes$Substrate),]

## interesting stuff about proteins that can be weighter into the graph
# number peptides
# number of CDS mutations
# number of amino acid mutations
# gene_CDS length (protein size)

proteins <- unique(c(nodes$Kinase, nodes$Substrate))
proteins <- data.frame(Substrate=proteins)
CDS_muts <- aggregate(Mutation_CDS~Substrate, data=mut_data, function(x) length(unique(x)))
AA_muts <- aggregate(Mutation_AA~Substrate, data=mut_data, function(x) length(unique(x)))
gene_CDS <- unique(mut_data[,c(1,30)])

# kin sub or both

kins <- unique(nodes$Kinase)
subs <- unique(nodes$Substrate)
for(i in 1:1472){
  if(proteins$Substrate[i] %in% kins){
    proteins$kin[i] <- 1}
  else
  {proteins$kin[i] <- 0}}

for(i in 1:1472){
  if(proteins$Substrate[i] %in% subs){
    proteins$sub[i] <- 2}
  else
  {proteins$sub[i] <- 0}}
  
proteins <- within(proteins, mode <- kin + sub)

library(dplyr)
prots <- inner_join(CDS_muts, AA_muts, by="Substrate")
prots <- inner_join(prots, gene_CDS, by="Substrate")
prots <- merge(proteins[,c(1,4)], prots, by="Substrate", all=T)

network <- aggregate(Peptide~Kinase*Substrate, data=nodes, length)
names(network) <-c("from", "to", "Peptide")


## SUB FOR BRAF
net <- network[network$from == "BRAF" | network$to =="BRAF",]
x <- unique(c(net$from, net$to))
net <- network[network$from %in% x | network$to %in% x,]
x <- unique(c(net$from, net$to))
braf <- merge(prots, data.frame(Substrate=x), by="Substrate", all.y=T)
braf[is.na(braf)] <- 0

BRAF <- graph.data.frame(net, directed=T, vertices = braf)

# library(intergraph)
# BRAFn <- asNetwork(BRAF)
list.vertex.attributes(BRAF)
list.edge.attributes(BRAF)

######
network2 <- as.data.frame(network[,c(1,2)])
prots2 <- as.data.frame(prots)
prots3 <- unique(prots2)
prots3$Substrate <- as.character(prots3$Substrate)
prots3$mode <- as.numeric(prots3$mode)
prots3$Mutation_AA <- as.numeric(prots3$Mutation_AA)
prots3$Mutation_CDS <- as.numeric(prots3$Mutation_CDS)
prots3$Gene_CDS_length <- as.numeric(prots3$Gene_CDS_length)
prots3[is.na(prots3)] <- 0

g <- graph.data.frame(network, directed=T, vertices = prots3)

########### LOAD GRAPH FILES

load("g.RData")
load("BRAF.RData")

list.vertex.attributes(g)
list.edge.attributes(g)


## all close to BRAF proteins
# BRAF <- c("PTEN", "PDK1", "AKT1", "TSC1", "MTOR", "RPS6", 
#               "EIF4EBP1", "NRAS", "KRAS", "BRAF", "MAP2K1", "MAP2K2", "MAPK3", "MAPK1")
# 
# 
# BRAFg <- induced.subgraph(g,BRAF)



is.connected(BRAF)
is.simple(BRAF)
clusters(BRAF)
diameter(BRAF)
max(degree(BRAF))


colbar <- c("red", "dodgerblue", "gold")

l <- layout.kamada.kawai(BRAF)
e.widths <- sqrt(E(BRAF)$Peptide)
#v.shapes <- c("circle", "square", "none")[as.factor(V(BRAF)$mode)]
v.shapes <- "circle"
v.sizes <- 12*sqrt(graph.strength(BRAF))  
v.sizes <- sqrt(V(BRAF)$Gene_CDS_length+1)/10
v.colors <- colbar[as.factor(V(BRAF)$mode)]
V(BRAF)$label.cex <- 0.6


plot(BRAF, layout=l,
     vertex.color=v.colors,
     vertex.shape = v.shapes,
     vertex.color = v.colors,
     edge.width = e.widths,
     vertex.size = v.sizes,
     edge.arrow.size = 0.1)


# add CDS lengths for all proteins

igraph.options(vertex.size=3, vertex.label.cex=2, vertex.label = V(BRAF)$vertex.names, edge.arrow.size=0.5)
par(mfrow=c(1,3))
plot(BRAF, layout=layout.circle, vertex.color = v.colors, vertex.shape = v.shapes, edge.arrow.size = 0.1, vertex.size = v.sizes, edge.width=e.widths)
plot(BRAF, layout=layout.fruchterman.reingold, vertex.color = v.colors, vertex.shape = v.shapes, edge.arrow.size = 0.1,vertex.size = v.sizes, edge.width=e.widths)
plot(BRAF, layout=layout.kamada.kawai, vertex.color = v.colors, vertex.shape = v.shapes, edge.arrow.size = 0.1,vertex.size = v.sizes, edge.width=e.widths)
dev.off()


## PHOSPHO ATLAS FULL NETWORK
e.widths <- sqrt(E(g)$Peptide)
v.shapes <- "circle"
v.sizes <- 12*sqrt(graph.strength(g))  
v.sizes <- sqrt(V(g)$Gene_CDS_length+1)/10
v.colors <- colbar[as.factor(V(g)$mode)]
V(g)$label.cex <- 0.6
V(g)$vertex.names <- V(g)

igraph.options(vertex.size=3, vertex.label.cex=2, vertex.label = V(g)$vertex.names, edge.arrow.size=0.5)
par(mfrow=c(1,3))
plot(g, layout=layout.circle, vertex.color = v.colors, vertex.shape = v.shapes, edge.arrow.size = 0.1, vertex.size = v.sizes, edge.width=e.widths)
plot(g, layout=layout.fruchterman.reingold, vertex.color = v.colors, vertex.shape = v.shapes, edge.arrow.size = 0.1,vertex.size = v.sizes, edge.width=e.widths)
plot(g, layout=layout.kamada.kawai, vertex.color = v.colors, vertex.shape = v.shapes, edge.arrow.size = 0.1,vertex.size = v.sizes, edge.width=e.widths)
dev.off()

# visualizing large networks

set.seed(42)
layout = layout.kamada.kawai(BRAF)
plot(BRAF, layout=l, vertex.color = v.colors, vertex.size = 3, edge.arrow.size = 0.1)

l <- layout.drl(BRAF)
plot(BRAF, layout = l, vertex.size = 3, vertex.label=NA, vertex.color = v.colors, edge.arrow.size = 0.1)

# describing the network

# degree distribution (for weighted networks)

hist(degree(g), col = "lightblue", xlab = "Vertex Degree", ylab = "Frequency", xlim = c(0, 120), breaks = 120,  main = "Degree Distribution of Phospho Atlas Vertices")

# vertex strength - summing up the weights of edges incident to a given vertex
hist(graph.strength(g), col = "pink", xlab = "Vertex Strength", xlim = c(0, 120), breaks = 120, ylab = "Frequency",  main = "Strength of Phospho Atlas Vertices")

#####BRAF
# describing the network
# degree distribution (for weighted networks)
hist(degree(BRAFn), col = "lightblue", xlab = "Vertex Degree", ylab = "Frequency", xlim = c(0, 120), breaks = 120,  main = "Degree Distribution of Phospho Atlas Vertices")
# vertex strength - summing up the weights of edges incident to a given vertex
hist(graph.strength(BRAF), col = "pink", xlab = "Vertex Strength", xlim = c(0, 120), breaks = 120, ylab = "Frequency",  main = "Strength of Phospho Atlas Vertices")
d.g <- degree(BRAFn)
dd.g <- degree.distribution(BRAFn)
d <- 1:max(d.g) - 1
ind <- (dd.g != 0)
plot(d[ind], dd.g[ind], log="xy",
     col = "blue", xlab = c("Log-Degree"),
     ylab = c("Log - Intensity"),
     main = "Log - Log Degree Distribution")
# understanding the manner in which vertices of different degrees are linked with each other
# degree of neighbors of a given vertex
g2 <- simplify(BRAF)
a.nn.deg.g2 <- graph.knn(g2, V(g2))$knn
d.g2 <- degree(BRAFn)
plot(d.g2, a.nn.deg.g2, log="xy",
     col = "goldenrod", xlab = c("Log Vertex Degree"),
     ylab = c("Log Average Neighbor Degree"))


###########


### CUTTING TOP OFF 
hist(degree(g), col = "lightblue", xlab = "Vertex Degree", ylab = "Frequency", ylim = c(0, 200), 
     xlim = c(0, 120), breaks = 120,  main = "Degree Distribution of Phospho Atlas Vertices")

# vertex strength - summing up the weights of edges incident to a given vertex
hist(graph.strength(g), col = "pink", xlab = "Vertex Strength", ylim = c(0,200),
     xlim = c(0, 120), breaks = 120, ylab = "Frequency",  main = "Strength of Phospho Atlas Vertices")

# given the nature of decay of this distribution , a log - log scale would be more usefulfor summarizin gthe degree information
d.g <- degree(g)
dd.g <- degree.distribution(g)
d <- 1:max(d.g) - 1
ind <- (dd.g != 0)
plot(d[ind], dd.g[ind], log="xy",
     col = "blue", xlab = c("Log-Degree"),
     ylab = c("Log - Intensity"),
     main = "Log - Log Degree Distribution")

# understanding the manner in which vertices of different degrees are linked with each other
# degree of neighbors of a given vertex
g2 <- simplify(g)
a.nn.deg.g2 <- graph.knn(g2, V(g2))$knn
d.g2 <- degree(g2)
plot(d.g2, a.nn.deg.g2, log="xy",
     col = "goldenrod", xlab = c("Log Vertex Degree"),
     ylab = c("Log Average Neighbor Degree"))

#

l <- layout.kamada.kawai(BRAF)
plot(BRAF, layout = l, main = "Hubs",
     vertex.size = 10 * sqrt(hub.score(BRAF)$vector))

plot(BRAF, layout = l, main = "Authorities", vertex.size=10*sqrt(authority.score(BRAF)$vector))

# for network of small to moderate size:

A <- get.adjacency(BRAF, sparse = F)
library(network)
g <- network::as.network.matrix(A)
library(sna)
sna::gplot.target(g, degree(g), main = "Degree",
                  circ.lab=F, circ.col = "skyblue",
                  usearrows = F, vertex.label.cex <- 0.6,
                  edge.col="darkgrey")

sna::gplot.target(g, closeness(g), main = "Closeness",
                  circ.lab=F, circ.col = "skyblue",
                  usearrows = F, 
                  vertex.col = c("blue", rep("red", 32), "yellow"),
                  edge.col="darkgrey")

sna::gplot.target(g, betweenness(g), main = "Betweenness",
                  circ.lab=F, circ.col = "skyblue",
                  usearrows = F, 
                  vertex.col = c("blue", rep("red", 32), "yellow"),
                  edge.col="darkgrey")

sna::gplot.target(g, evcent(g), main = "Eigenvector Centrality",
                  circ.lab=F, circ.col = "skyblue",
                  usearrows = F, 
                  vertex.col = c("blue", rep("red", 32), "yellow"),
                  edge.col="darkgrey")


## edges - betweenness
eb <- edge.betweenness(BRAF)
E(BRAF)[order(eb, decreasing=T)[1:10]]


## subgraphs and censuses

table(sapply(cliques(BRAF),length))

# 198 cliques of size 1 and 217 cliques of size 2

cliques(BRAF)[sapply(cliques(BRAF), length)==5]
table(sapply(maximal.cliques(BRAF), length))
clique.number(BRAF)

#coreness
A <- get.adjacency(BRAF, sparse = F)
library(network)
g <- network::as.network.matrix(A)
library(sna)
cores <- graph.coreness(BRAF)
sna::gplot.target(g, cores, circl.lab=F,
                  circ.col='skyblue', usearrows=F,
                  vertex.col=cores, edge.col="darkgrey")
detach("package:network")
detach("package:sna")

BRAFs <- simplify(BRAF)
dyad.census(BRAFs)

# 2 pairs with mutual connections in BRAF
# 215 pairs with non-mutual connections
# 19286 pairs with no connections between them


#measure of global clustering - frequency with which connected triples close to form triangles (percentage of triplets that form triangles)
transitivity(BRAF)
# 0 

#extent to which reciprocation exists among ties in a directed network; repciprocity: a total number of reciprocated edges divided by the total number of edges

reciprocity(BRAF, mode="default")
reciprocity(BRAF, mode="ratio")


transitivity(g)
reciprocity(g, mode="default")
reciprocity(g, mode="ratio")

table(sapply(cliques(g),length))
# 1028 cliques of 1, 1851 cliques of 2, 180 cliques of 4, 4 cliques of 4
gs <- simplify(g)
dyad.census(gs)
# 9 pairs with mutual connections in phospho atlas
# 1842 pairs with non-mutual connections
# 526027 pairs with no connections between them

is.connected(g)
comps <- decompose.graph(g)
table(sapply(comps, vcount))

# one cluster of 956, three clusters of 2, 20 clusters of 2 and 26 clusters of 1

g.gc <- decompose.graph(g)[[1]]
average.path.length(g.gc)
diameter(g.gc)
transitivity(g.gc)

# small clustering, and the shortest path distance is small

vertex.connectivity(g.gc)
edge.connectivity(g.gc)
# it requires the removal of only a single well chosen vertex or edge in order to break the subgraph into additional components


#a single vertwx that disconnects the graph is called a cut vertex or an articulatin point - where is th enetwork vulnerable
g.cut.vertices <- articulation.points(g.gc)
length(g.cut.vertices)

V(g)[g.cut.vertices]

BRAF.gc <- decompose.graph(BRAF)[[1]]
BRAF.cut.vertices <- articulation.points(BRAF.gc)
length(BRAF.cut.vertices)

V(BRAF)[BRAF.cut.vertices]

# what are the strongly and weakly connected components?

is.connected(g)
is.connected(g, mode=c("weak"))
is.connected(g, mode=c("strong"))

is.connected(BRAF)
is.connected(BRAF, mode=c("weak"))
is.connected(BRAF, mode=c("strong"))

BRAF.scc <- clusters(BRAF, mode=c("strong"))
table(BRAF.scc$csize)

# There are two stronly connected clusters in BRAF network and theyre of size 2


## GRAP PARTITIONING

# partitioning is a tool for finding in an unsupervised fashion
# subset of vertices that demonstrate a cohesiveness with respect to the underlying relational patterns

# hierarchical clustering

BRAF2 <- graph.data.frame(net, directed=F, vertices = braf)
bc <- fastgreedy.community(BRAF2)
### NOT WORKING

#spectral partitioning

b.lap <- graph.laplacian(BRAF)
eig.anal <- eigen(b.lap)
plot(eig.anal$values, col="blue",
     ylab="Eigenvalues of Graph Laplacian")
# there are many eigenvalues exactly equal to zero (bc network is not connected); 

f.vec <- eig.anal$vectors[,33]

plot(f.vec, pch=16, xlab = "Actor Number", ylab = "Fiedler Vector Entry")
abline(0, 0, lwd=2, col="lightgrey")

# validation

func.class <- get.vertex.attribute(g.gc, "mode")
table(func.class)

# requires undirected network
gc <- fastgreedy.community(g.gc)
c.m <- membership(gc)

table(c.m, func.class, useNA=c("no"))



assortativity.nominal(BRAF, (V(BRAF)$mode==1)+1, directed=F)
#for the class  1 representing proteins that are only kinase we see an assortativity coefficient -0.009174

# Pearson correlation coefficient for when the vertex characteristic of interest is continuous

# degree-degree correlation of adjacent vertices

assortativity.degree(BRAF)
# degree correlation positive and very large


## math modeling

nv <- vcount(BRAF)
ne <- ecount(BRAF)
degs <- degree(BRAF)
ntrials <- 1000

num.comm.rg <- numeric(ntrials)
for( i in (1:ntrials)){
  g.rg <- erdos.renyi.game(nv, ne, type = "gnm")
  c.rg <- fastgreedy.community(g.rg)
}

num.comm.grg <- numeric(ntrials)
for( i in (1:ntrials)){
  g.grg <- degree.sequence.game(degs, method = "vl")
  c.grg <- fastgreedy.community(g.grg)
  num.comm.grg[i] <- length(c.grg)
}

rslts <- c(num.comm.rg, num.comm.grg)
indx <- c(rep(0,ntrials), rep(1,ntrials))
counts <- table(indx,rslts)/ntrials
barplot(counts, beside=T, col=c("blue", "red"),
        xlab="Number of Communities",
        ylab="Relative Frequency",
        legend=c("Fixed Size", "Fixed Degree Sequence"))
# distribution of number of communities detected for randm graphs of the same size (blue) and degree sequence(red) as th eBRAF network


### CHAPTER 6 Statistical modeling

A <- get.adjacency(BRAF)
v.attrs <- get.data.frame(BRAF, what="vertices")
e.attrs <- get.data.frame(BRAF, what="edges")

library(ergm)
BRAF.s <- network::as.network(as.matrix(A), directed=F)

list.vertex.attributes(BRAF)
list.edge.attributes(BRAF)


network::set.vertex.attribute(BRAF.s, "Mode", v.attrs$mode)
network::set.vertex.attribute(BRAF.s, "Size", v.attrs$Gene_CDS_length)
network::set.vertex.attribute(BRAF.s, "Mutation_CDS", v.attrs$Mutation_CDS)
network::set.vertex.attribute(BRAF.s, "Mutation_AA", v.attrs$Mutation_AA)
network::set.edge.attribute(BRAF.s, "Peptide", e.attrs$Peptide)

my.ergm.bern <- formula(BRAF.s~edges)
my.ergm.bern
summary.statistics(my.ergm.bern)
# 217 edges

my.ergm <- formula(BRAF.s ~ edges + kstar(2) + kstar(3) + triangle)
summary.statistics(my.ergm)
# 217 edges, 9313 kstar2, 3181376 kstar3

my.ergm <- formula(BRAF.s ~ edges +
                     gwesp(1, fixed=T))
summary.statistics(my.ergm)


BRAF.ergm <- formula(BRAF.s ~ edges + 
                         gwesp(log(3), fixed=T) +
                         nodemain("Mode")+
                         nodemain("Size")+
                         match("Size")+
                         match("Mode")+
                         match("Mutation_CDS"))



summary.statistics(BRAF.ergm)


set.seed(42)
BRAF.ergm.fit <- ergm(BRAF.ergm)
mcmc.diagnostics(BRAF.ergm.fit)
summary.ergm(BRAF.ergm.fit)

# goodnes of fit
gof.BRAF.ergm <- gof(BRAF.ergm.fit)
par(mfrow=c(1,3))
plot(gof.BRAF.ergm)

#stochastic block model
library(mixer)
setSeed(42)
BRAF.sbm <- mixer(as.matrix(get.adjacency(BRAF)),
                   qmin=2, qmax=15)

BRAF.sbm.output <- getModel(BRAF.sbm)
names(BRAF.sbm.output)
BRAF.sbm.output$q
BRAF.sbm.output$alphas
BRAF.sbm.output$Taus[,1:3]
my.ent <- function(x) {-sum(x*log(x,2))}
apply(BRAF.sbm.output$Taus[,1:3],2,my.ent)
log(BRAF.sbm.output$q,2)
summary(apply(BRAF.sbm.output$Taus, 2, my.ent))

plot(BRAF.sbm, classes=as.factor(V(BRAF)$mode))




V(BRAF)[mode==1]$color <- "orangered"
V(BRAF)[mode==2]$color <- "dodgerblue"
V(BRAF)[mode==3]$color <- "goldenrod"

plot(BRAF, vertex.size=5, vertex.label=NA)

# nearest neighbor prediction

clu <- clusters(BRAF)
BRAF.gc <- induced.subgraph(BRAF, clu$membership == which.max(clu$csize))
nn.ave <- sapply(V(BRAF.gc), function(x) mean(V(BRAF.gc)[nei(x)]$mode))

par(mfrow=c(3,1))
hist(nn.ave[V(BRAF.gc)$mode == 1], col = "orangered",
     xlab = "Proportion neighbors that are a kinase",
     main = "Egos that are a kinase")
hist(nn.ave[V(BRAF.gc)$mode == 2], col = "dodgerblue",
     xlab = "Proportion neighbors that are a substrate",
     main = "Egos that are a substrate")
hist(nn.ave[V(BRAF.gc)$mode == 3], col = "goldenrod",
     xlab = "Proportion neighbors that are both a kinase and a substrate",
     main = "Egos that are a both a kinase and a substrate")


nn.pred <- as.numeric(nn.ave > 0.5)
mean(as.numeric(nn.pred != V(BRAF.gc)$mode))
# ^ error rate of roughly 98% given a treshold of 0.5
