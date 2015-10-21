# notes based on "Statistical Analysis of Network Data with R" by Kolaczek and Csardi

# setwd

library(igraph)

g <- graph.formula(1-2,
                   1-3,
                   2-3,
                   2-4,
                   3-5,
                   4-5,
                   4-6,
                   4-7,5-6,6-7)

V(g)
E(g)
str(g)

plot(g)

#directed

dg <- graph.formula(1-+2, 1-+3, 2++3)
plot(dg)

dg <- graph.formula(Sam-+Mary, Sam-+Tom, Mary++Tom)
str(dg)

V(dg)$name <- c("Sam", "Mary", "Tom")

E(dg)
get.edgelist(dg)


h <- induced.subgraph(g, 1:5)
str(h)

h <- g - vertices(c(6,7))

h <- h + vertices(c(6,7))
g <- h + edges(c(4,6), c(4,7), c(5,6), c(6,7))

h1 <- h
h2 <- graph.formula(4-6, 4-7, 5-6, 6-7)
g <- graph.union(h1,h2)

V(dg)$name
V(dg)$gender <- c("M", "F", "M")
V(dg)$color <- "red"

is.weighted(g)

wg <-g
E(wg)$weight <- runif(ecount(wg))
is.weighted(wg)
g$name <- "Toy Graph"


library(sand)
g.lazega <- graph.data.frame(elist.lazega, directed="FALSE", vertices=v.attr.lazega)
g.lazega$name <- "Lazega Lawyers"
vcount(g.lazega)
ecount(g.lazega)
list.vertex.attributes(g.lazega)
is.simple(g)

mg <- g + edge(2,3)
str(mg)
is.simple(mg)

E(mg)$weight <- 1
wg2 <- simplify(mg)
is.simple(wg2)

str(wg2)
E(wg2)$weight

### DESCRIBING CONNECTIVITY

# 1. adjacency

neighbors(g,5)

# degree of a vertex - number of edges incident on the vertex

degree(g)

degree(dg, mode="in")
degree(dg, mode="out")

# a walk on a graph is an alternating sequence of vertix-edge, where the endpoints of current edge are {previous edge, current edge}
# refinements of a walk include trails, which are walks without repeated edges and paths which are trails without repeated vertices
# a trail for which the beginning and ending vertices are the same is called a circuit
# a walk of length at least three for which th ebeginning and ending vertices are the same, but for which all other vertices are distinct from each other is called a cycle
# graphs containing no cycles are acyclic
# a vertext is said to be reachable from another vertex if there exists a walk
# the graph is said to be connected if every verted is reachable from every other
# a component of a graph is a maximally connected subgraph
is.connected(g)
clusters(g)
is.connected(dg, mode="weak")
is.connected(dg, mode="strong")

# distance between vertices is defined as the length of the shortest path between the vertices -- geodesic distance
# the value of the longest distance in a graph is called a diameter
diameter(g, weights=NA)

#### TYPES OF GRAPHS

g.full <- graph.full(7)
g.ring <- graph.ring(7)
g.tree <- graph.tree(7, children=2, mode="undirected")
g.star <- graph.star(7, mode="undirected")

plot(g.full)
plot(g.ring)
plot(g.tree)
plot(g.star)

# a complete graph is a graph where every vertex is joined to evert other vertex by an edge
# a clique is a complete subgraph
# a regular graph is a graph in which every vertex has the same degree
# a connected graph with no cycles is called a tree
# a disjoint union of such graphs is called a forest
# (...)

# bipartite graph has a vertex set that may be partitioned into two disjoint sets and edge has one endpoint in one set and the other in the other set

g.bip <- graph.formula(actor1:actor2:actor3,
                       movie1:movie2, actor1:actor2 - movie1,
                       actor2:actor3 - movie2)
V(g.bip)
str(g.bip, v=T)

plot(g.bip)

## CHAPTER 3

# layouts

library(sand)
g.l <- graph.lattice(c(5,5,5))

data(aidsblog)
summary(aidsblog)

# circular layout
igraph.options(vertex.size=3, vertex.label=NA, edge.arrow.size=0.5)
par(mfrow=c(1,2))
plot(g.l, layout=layout.circle)
title("5x5x5 Lattice")
plot(aidsblog, layout=layout.circle)
title("Blog Network")
dev.off()

# Fruchterman and Reingold method
par(mfrow=c(1,2))
plot(g.l, layout=layout.fruchterman.reingold)
title("5x5x5 Lattice")
plot(aidsblog, layout=layout.fruchterman.reingold)
title("Blog Network")
dev.off()

# multidimensional scaling

par(mfrow=c(1,2))
plot(g.l, layout=layout.kamada.kawai)
title("5x5x5 Lattice")
plot(aidsblog, layout=layout.kamada.kawai)
title("Blog Network")
dev.off()


g.tree <- graph.formula(1-+2, 1-+3, 1-+4, 2-+5, 2-+6, 2-+7, 3-+8, 3-+9, 4-+10)
par(mfrow=c(1,3))
igraph.options(vertex.size=30, edge.arrow.size=0.5, vertex.label=NULL)
plot(g.tree, layout=layout.circle)
plot(g.tree, layout=layout.reingold.tilford(g.tree, circular=T))
plot(g.tree, layout=layout.reingold.tilford)



plot(BRAFg, layout=layout.circle)
plot(BRAFg, layout=layout.reingold.tilford)
# 
# plot(g.bip, layout=-layout.bipartite(g.bip)[,2:1],
#      vertex.size=30, vertex.shape=ifelse(V(g.bip)$type, "rectangle", "circle"),
#      vertex.color=ifelse(V(g.bip)$type,"red", "cyan"))

library(igraphdata)
data(karate)
set.seed(42)
l <- layout.kamada.kawai(karate)
igraph.options(vertex.size=10)

par(mfrow=c(1,1))
plot(karate,layout=l, vertex.label=V(karate))
V(karate)$label <- sub("Actor", "", V(karate)$name)
V(karate)$shape <- "circle"
V(karate)[c("Mr Hi", "John A")]$shape <- "rectangle"
V(karate)[Faction == 1]$color <- "red"
V(karate)[Faction == 2]$color <- "dodgerblue"
V(karate)$size <- 4*sqrt(graph.strength(karate))
V(karate)$size2 <- V(karate)$size * .5
E(karate)$width <- E(karate)$weight
F1 <- V(karate)[Faction ==1]
F2 <- V(karate)[Faction ==2]
E(karate)[ F1 %--% F1 ]$color <- "pink"
E(karate)[ F2 %--% F2 ]$color <- "lightblue"
E(karate)[ F1 %--% F2 ]$color <- "yellow"
V(karate)$label.dist <- ifelse(V(karate)$size >= 10, 0, 0.75)
plot(karate, layout=l)



library(sand)
data(lazega)
colbar <- c("red", "dodgerblue", "goldenrod")
v.colors <- colbar[V(lazega)$Office]
v.shapes <- c("circle", "square")[V(lazega)$Practice]
v.size <- 3.5*sqrt(V(lazega)$Years)
v.label <- V(lazega)$Seniority
set.seed(42)
l <- layout.fruchterman.reingold(lazega)
plot(lazega, layout = l, vertex.color = v.colors, 
     vertex.shape = v.shapes, 
     vertex.size = v.size,
     vertex.label = v.label)


## Visualizing large networks

library(sand)
library(igraph)
summary(fblog)
party.names <- sort(unique(V(fblog)$PolParty))
party.names
set.seed(42)
l <- layout.kamada.kawai(fblog)
party.nums.f <- as.factor(V(fblog)$PolParty)
party.nums <- as.numeric(party.nums.f)
plot(fblog, layout = l, vertex.label=NA, vertex.color = party.nums, vertex.size=3)


set.seed(42)
l <- layout.drl(fblog)
plot(fblog, layout = l, vertex.size = 5, vertex.label = NA, vertex.color = party.nums)

# contract vertices

# fblog.c <- contract.vertices(fblog, party.nums)
# E(fblog.c)$weight <- 1
# fblog.c <- simplify(fblog.c)
# party.size <- as.vector(table(V(fblog)$PolParty))
# plot(fblog.c, vertex.size=5*sqrt(party.size),
#      vertex.label = party.names,
#      vertex.color = V(fblog.c),
#      edge.width = sqrt(E(fblog.c)$weight),
#      vertex.lavel.dist = 1.5,
#      edge.arrow.size = 0)



data(karate)
k.nbhds <- graph.neighborhood(karate, order=1)
sapply(k.nbhds, vcount)
k.1 <- k.nbhds[[1]]
k.34 <- k.nbhds[[34]]
par(mfrow=c(1,2))
plot(k.1, vertex.label=NA,
     vertex.color = c("red", rep("lightblue", 16)))
plot(k.34, vertex.label=NA,
     vertex.color = c(rep("lightblue", 17), "red"))


# Chapter 4

library(sand)
data(karate)
hist(degree(karate), col = "lightblue", xlim = c(0,50),
     xlab="Vertex Degree", ylab = "Frequency", main = "")

hist(graph.strength(karate), col = "pink", 
     xlab="Vertex Strength", ylab = "Frequency", main = "")

library(igraphdata)
data(yeast)
d.yeast <- degree(yeast)
hist(d.yeast, col = "blue",
     xlab = "Degree", ylab = "Frequency",
     main = "Degree Distribution")

dd.yeast <- degree.distribution(yeast)
d <- 1:max(d.yeast) - 1
ind <- (dd.yeast != 0)
plot(d[ind], dd.yeast[ind], log="xy",
     col = "blue", xlab = c("Log-Degree"),
     ylab = c("Log - Intensity"),
     main = "Log - Log Degree Distribution")


# understanding the manner in which vertices of different degrees are linked with each other
# degree of neighbors of a given vertex

a.nn.deg.yeast <- graph.knn(yeast, V(yeast))$knn
plot(d.yeast, a.nn.deg.yeast, log="xy",
     col = "goldenrod", xlab = c("Log Vertex Degree"),
     ylab = c("Log Average Neighbor Degree"))
# there is a tendency for vertices of higher degrees to link with similar vertices but those of lower degrees tend to connect to both low and high degree vertices

# vertex centrality

# describe importance of vertex in a network; vertex degree
# closeness, betweenness and eigenvector centrality

## closeness centrality - vertex is central if its close to many other vertices
## betweenness centrality - summarize the extent to which a vertex is located between other pairs of vertices

# for network of small to moderate size:

A <- get.adjacency(karate, sparse = F)
library(network)
g <- network::as.network.matrix(A)
library(sna)
sna::gplot.target(g, degree(g), main = "Degree",
                  circ.lab=F, circ.col = "skyblue",
                  usearrows = F, 
                  vertex.col = c("blue", rep("red", 32), "yellow"),
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


l <- layout.kamada.kawai(aidsblog)
plot(aidsblog, layout = l, main = "Hubs", vertex.label = "",
     vertex.size = 10 * sqrt(hub.score(aidsblog)$vector))

plot(aidsblog, layout = l, main = "Authorities",
     vertex.label="", vertex.size=10*sqrt(authority.score(aidsblog)$vector))



## edges - betweenness
eb <- edge.betweenness(karate)
E(karate)[order(eb, decreasing=T)[1:3]]

## subgraphs and censuses

table(sapply(cliques(karate),length))

cliques(karate)[sapply(cliques(karate), length)==5]
table(sapply(maximal.cliques(karate), length))
clique.number(yeast)

cores <- graph.coreness(karate)
sna::gplot.target(g, cores, circl.lab=F,
                  circ.col='skyblue', usearrows=F,
                  vertex.col=cores, edge.col="darkgrey")


detach("package:network")
detach("package:sna")


aidsblog <- simplify(aidsblog)
dyad.census(aidsblog)

# density and related notions of relative frequency

ego.instr <- induced.subgraph(karate, neighborhood(karate,1,1)[[1]])
ego.admin <- induced.subgraph(karate, neighborhood(karate,1,34)[[1]])
graph.density(karate)
graph.density(ego.instr)
graph.density(ego.admin)

#measure of global clustering - frequency with which connected triples close to form triangles (percentage of triplets that form triangles)
transitivity(karate)

transitivity(karate, "local", vids = c(1,34))

#extent to which reciprocation exists among ties in a directed network; repciprocity: a total number of reciprocated edges divided by the total number of edges

reciprocity(aidsblog, mode="default")
reciprocity(aidsblog, mode="ratio")

# connectivity, cuts and flows

# does a given graph separate into distinct subgraphs

is.connected(yeast)
comps <- decompose.graph(yeast)
table(sapply(comps, vcount))

yeast.gc <- decompose.graph(yeast)[[1]]
average.path.length(yeast.gc)
diameter(yeast.gc)
transitivity(yeast.gc)

vertex.connectivity(yeast.gc)
edge.connectivity(yeast.gc)
# it requires the removal of only a single well chosen vertex or edge in order to break the subgraph into additional components

yeast.cute.vertices <- articulation.points(yeast.gc)
length(yeast.cute.vertices)


## GRAP PARTITIONING

# partitioning is a tool for finding in an unsupervised fashion
# subset of vertices that demonstrate a cohesiveness with respect to the underlying relational patterns

# hierarchical clustering

kc <- fastgreedy.community(karate)
length(kc)
sizes(kc)
membership(kc)
plot(kc, karate)

library(ape)
dendPlot(kc, mode="phylo")


#spectral partitioning

k.lap <- graph.laplacian(karate)
eig.anal <- eigen(k.lap)
plot(eig.anal$values, col="blue",
     ylab="Eigenvalues of Graph Laplacian")
# there is only one eigenvalue exactly equal to zero (bc network is connected); 
# second smallest eigenvalue is quite close to 0

f.vec <- eig.anal$vectors[,33]

faction <- get.vertex.attribute(karate, "Faction")
f.colors <- as.character(length(faction))
f.colors[faction ==1] <- "red"
f.colors[faction ==2] <- "cyan"
plot(f.vec, pch=16, xlab = "Actor Number", ylab = "Fiedler Vector Entry", col = f.colors)
abline(0, 0, lwd=2, col="lightgrey")

# validation

func.class <- get.vertex.attribute(yeast.gc, "Class")
table(func.class)

yc <- fastgreedy.community(yeast.gc)
c.m <- membership(yc)

table(c.m, func.class, useNA=c("no"))

assortativity.nominal(yeast, (V(yeast)$Class=="P")+1, directed=F)
#for the class "P" representing proteins that are known to play a role in protein synthesis we see an ssortativity coefficient of nearly 0.5

# Pearson correlation coefficient for when the vertex characteristic of interest is continuous

# degree-degree correlation of adjacent vertices

assortativity.degree(yeast)
# degree correlation positive and very large


## CHAPTER 5
## Mathematical Models for Network Graphs

library(sand)
set.seed(42)
g.er <- erdos.renyi.game(100, 0.02)
plot(g.er, layout = layout.circle, vertex.label = NA)
is.connected(g.er)
table(sapply(decompose.graph(g.er), vcount))
mean(degree(g.er))

hist(degree(g.er), col = "lightblue",xlab="Degree", ylab="Frequency", main="")
#degree distribution is quite homogenous
average.path.length(g.er)
diameter(g.er)
transitivity(g.er)

# generalized random graph models


degs <- c(2,2,2,2,3,3,3,3)
g1 <- degree.sequence.game(degs, method = "vl")
g2 <- degree.sequence.game(degs, method="vl")

plot(g1, vertex.label=NA)
plot(g2, vertex.label=NA)

graph.isomorphic(g1, g2)


data(karate)
nv <- vcount(karate)
ne <- ecount(karate)
degs <- degree(karate)
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


## macaques

library(igraphdata)
data(macaque)
summary(macaque)

#average over all vertices, of the vertex-specific clustering coefficients:
clust.coef.dir <- function(graph){
  A <- as.matrix(get.adjacency(graph))
  S <- A + t(A)
  deg <- degree(graph, mode=c("total"))
  num <- diag(S %*% S %*% S)
  denom <- diag(A %*% A)
  denom <- 2*(deg*(deg-1)-2*denom)
  cl <- mean(num/denom)
  return(cl)
}

# ERROR here:
# path lengths
ntrials <- 1000
nv <- vcount(macaque)
ne <- ecount(macaque)
cl.rg <- numeric(ntrials)
apl.rg <- numeric(ntrials)
for(i in (1:ntrials)){
  g.rg <- erdos.renyi.game(nv, ne, type="gnm", directed=T)
  cl.rg[i] <- clust.coef.dir(g.rg)
  apl.rg[i] <- average.path.length(g.rg)
}

# summarizing the resulting distributions of clustering coefficient:
summary(cl.rg)
#average path length
summary(apl.rg)
clust.coef.dir(macaque)

### CHAPTER 6 Statistical modeling

library(sand)
data(lazega)
A <- get.adjacency(lazega)
v.attrs <- get.data.frame(lazega, what="vertices")

library(ergm)
lazega.s <- network::as.network(as.matrix(A), directed=F)
network::set.vertex.attribute(lazega.s, "Office", v.attrs$Office)
network::set.vertex.attribute(lazega.s, "Practice", v.attrs$Practice)
network::set.vertex.attribute(lazega.s, "Gender", v.attrs$Gender)
network::set.vertex.attribute(lazega.s, "Seniority", v.attrs$Seniority)

my.ergm.bern <- formula(lazega.s~edges)
my.ergm.bern
summary.statistics(my.ergm.bern)

my.ergm <- formula(lazega.s ~ edges + kstar(2) + kstar(3) + triangle)
summary.statistics(my.ergm)
 

my.ergm <- formula(lazega.s ~ edges +
                   gwesp(1, fixed=T))
summary.statistics(my.ergm)

lazega.ergm <- formula(lazega.s ~ edges + 
                         gwesp(log(3), fixed=T) +
                         nodemain("Seniority")+
                         nodemain("Practice")+
                         match("Practice")+
                         match("Gender")+
                         match("Office"))

# it allows to assess effect on the formation of collaborative ties among lawyers that is had by
# seniority, type of practice and commonality of practice, gender and location of office

set.seed(42)
lazega.ergm.fit <- ergm(lazega.ergm)
summary.ergm(lazega.ergm.fit)

gof.lazega.ergm <- gof(lazega.ergm.fit)
par(mfrow=c(1,3))
plot(gof.lazega.ergm)

library(mixer)
setSeed(42)
fblog.sbm <- mixer(as.matrix(get.adjacency(fblog)),
                   qmin=2, qmax=15)

fblog.sbm.output <- getModel(fblog.sbm)
names(fblog.sbm.output)
fblog.sbm.output$q
fblog.sbm.output$alphas
fblog.sbm.output$Taus[,1:3]
my.ent <- function(x) {-sum(x*log(x,2))}
apply(fblog.sbm.output$Taus[,1:3],2,my.ent)
log(fblog.sbm.output$q,2)
summary(apply(fblog.sbm.output$Taus, 2, my.ent))

plot(fblog.sbm, classes=as.factor(V(fblog)$PolParty))

library(eigenmodel)
set.seed(42)
A <- get.adjacency(lazega, sparse=F)
lazega.leig.fit1 <- eigenmodel_mcmc(A, R=2, S=11000, burn=1000)
same.prac.op <- v.attr.lazega$Practice %o%
  v.attr.lazega$Practice
same.prac <- matrix(as.numeric(same.prac.op %in% c(1,4,9)), 36, 36)
same.prac <- array(same.prac, dim = c(36, 36, 1))
lazega.leig.fit2 <- eigenmodel_mcmc(A, same.prac, R=2, S=11000, burn = 10000)
same.off.op <- v.attr.lazega$Office %o%
  v.attr.lazega$Office
same.off <- matrix(as.numeric(same.off.op %in% c(1,4,9)), 36, 36)
same.off <- array(same.off, dim=c(36,36,1))
lazega.leig.fit3 <- eigenmodel_mcmc(A, same.off, R=2, S=11000, burn = 10000)

lat.sp.1 <- eigen(lazega.leig.fit1$ULU_postmean)$vec[,1:2]
lat.sp.2 <- eigen(lazega.leig.fit2$ULU_postmean)$vec[,1:2]
lat.sp.3 <- eigen(lazega.leig.fit3$ULU_postmean)$vec[,1:2]

colbar <- c("red", "dodgerblue", "goldenrod")
v.colors <- colbar[V(lazega)$Office]
v.shapes <- c("circle", "square")[V(lazega)$Practice]
v.size <- 3.5*sqrt(V(lazega)$Years)
v.label <- V(lazega)$Seniority

plot(lazega, layout=lat.sp.1, vertex.color = v.colors, vertex.shape = v.shapes, vertex.size = v.size, vertex.label = v.label)

apply(lazega.leig.fit1$L_postsamp,2,mean)
apply(lazega.leig.fit2$L_postsamp,2,mean)
apply(lazega.leig.fit3$L_postsamp,2,mean)


# Chapter 7  Network topology inference

library(sand)
nv <- vcount(fblog)
ncn <- numeric()
A <- get.adjacency(fblog)

# not working:
for(i in (1:(nv-1))){
  ni <- neighborhood(fblog,1,i)
  nj <- neighborhood(fblog, 1, (i+1):nv)
  nbhd.ij <- mapply(intersect, ni, nj, SIMPLIFY=FALSE)
  temp <- unlist(lapply(nbhd.ij, length)) - 2*A[i, (i+1):nv]
  ncn <- c(ncn, temp)
}


data(Ecoli.data)
ls()

heatmap(scale(Ecoli.expr), Rowv=NA)

library(igraph)
g.regDB <- graph.adjacency(regDB.adj, "undirected")
summary(g.regDB)
plot(g.regDB, vertex.size=3, vertex.label=NA)
mycorr <- cor(Ecoli.expr)
z <- 0.5 * log((1 + mycorr) / (1-mycorr))
z.vec <- z[upper.tri(z)]
n <- dim(Ecoli.expr)[1]
corr.pvals <- 2 * pnorm(abs(z.vec),0, sqrt(1/(n-3)), lower.tail=F)
length(corr.pvals)
corr.pvals.adj <- p.adjust(corr.pvals, "BH")
length(corr.pvals.adj[corr.pvals.adj < 0.05])

library(fdrtool)
mycorr.vec <- mycorr[upper.tri(mycorr)]
fdr <- fdrtool(mycorr.vec, statistic = "correlation")

pcorr.pvals <- matrix(0, dim(mycorr)[1], dim(mycorr)[2])
for(i in seq(1, 153)){
  for(j in seq(1,153)){
    rowi <- mycorr[i, -c(i,j)]
    rowj <- mycorr[j, -c(i,j)]
    tmp <- (mycorr[i,j] - rowi*rowj)/sqrt((1-rowi^2) * (1-rowj^2))
    tmp.zvals <- (0.5)*log((1+tmp) / (1-tmp))
    tmp.s.zvals <- sqrt(n-4)*tmp.zvals
    tmp.pvals <- 2* pnorm(abs(tmp.s.zvals),
                          0,1,lower.tail=F)
    pcorr.pvals[i,j] <- max(tmp.pvals)
  }
}

pcorr.pvals.vec <- pcorr.pvals[lower.tri(pcorr.pvals)]
pcorr.pvals.adj <- p.adjust(pcorr.pvals.vec, "BH")
pcorr.edges <- (pcorr.pvals.adj < 0.05)
length(pcorr.pvals.adj[pcorr.edges])

pcorr.A <- matrix(0,153,153)
pcorr.A[lower.tri(pcorr.A)] <- as.numeric(pcorr.edges)
g.pcorr <- graph.adjacency(pcorr.A, "undirected")

str(graph.intersection(g.regDB, g.pcorr, byname=F))

fdr <- fdrtool(pcorr.pvals.vec, statistic="pvalue")
pcorr.edges.2 <- (fdr$qval < 0.05)
length(fdr$qval[pcorr.edges.2])



data(sandwichprobe)
head(delaydata)
SSDelayDiff <- with(delaydata, by(DelayDiff^2, list(SmallPktDest, BigPktDest), sum))
head(SSDelayDiff)

x <- as.dist(1/sqrt(SSDelayDiff))
myclust <- hclust(x, method = "average")
plot(myclust, labels=host.locs, axes=F, ylab=NULL, ann=F)





set.seed(42)
library(sand)
data(ppi.CC)

summary(ppi.CC)
V(ppi.CC)$ICSC[1:10]
V(ppi.CC)[ICSC==1]$color <- "yellow"
V(ppi.CC)[ICSC==0]$color <- "blue"
plot(ppi.CC, vertex.size=5, vertex.label=NA)

# nearest neighbor prediction

clu <- clusters(ppi.CC)
ppi.CC.gc <- induced.subgraph(ppi.CC, clu$membership == which.max(clu$csize))
nn.ave <- sapply(V(ppi.CC.gc), function(x) mean(V(ppi.CC.gc)[nei(x)]$ICSC))

par(mfrow=c(2,1))
hist(nn.ave[V(ppi.CC.gc)$ICSC == 1], col = "yellow",
     ylim = c(0,30), 
     xlab = "Proportion neighbors w/ ICSC",
     main = "Egos w/ ICSC")
hist(nn.ave[V(ppi.CC.gc)$ICSC == 0], col = "blue",
     ylim = c(0,30), 
     xlab = "Proportion neighbors w/ ICSC",
     main = "Egos w/out ICSC")
nn.pred <- as.numeric(nn.ave > 0.5)
mean(as.numeric(nn.pred != V(ppi.CC.gc)$ICSC))
# ^ error rate of roughly 25% given a treshold of 0.5


source("http://bioconductor.org/biocLite.R")
biocLite("GOstats", suppressAutoUpdate=T, suppressUpdates=T)
library(GOstats)
library(GO.db)
biocLite("org.Sc.sgd.db", suppressAutoUpdate = T, suppressUpdates=T)
library(org.Sc.sgd.db)

x <- as.list(org.Sc.sgdGO2ALLORFS)
current.icst <- x[names(x) == "GO:0035556"]
ev.code <- names(current.icst[[1]])
icst.ida <- current.icst[[1]][ev.code == "IDA"]

orig.icsc <- V(ppi.CC.gc)[ICSC == 1]$name
candidates <- intersect(icst.ida, V(ppi.CC.gc)$name)
new.icsc <- setdiff(candidates, orig.icsc)
new.icsc

nn.ave[V(ppi.CC.gc)$name %in% new.icsc]


library(ngspatial)

X <- V(ppi.CC.gc)$ICSC
A <- get.adjacency(ppi.CC.gc, sparse = F)
formula1 <- X~1
