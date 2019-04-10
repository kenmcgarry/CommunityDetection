# dm_randomwalk.R
# https://sna.stanford.edu/lab.php?l=3
# https://www.r-bloggers.com/matrix-factorization/


library(igraph)
library(diffusr)

#making the graph
m1 <- sample(LETTERS, 15000, replace = T)
m2 <- sample(LETTERS, 15000, replace = T)
w  <- runif(15000, .1, 3)

graph <- data.frame(m1, m2, w, stringsAsFactors = F)
graph <- graph[graph$m1 != graph$m2,]

nodes <- list()

#parse the structure
for(node in unique(graph$m1)){
  temp <- list(target = c(), weight = c())
  df <- graph[graph$m1 == node,]
  
  if(nrow(df) == 0){
    next
  }
  
  for(i in 1:nrow(df)){
    temp$target <- c(temp$target, df[i,]$m2)
    temp$weight <- c(temp$weight, df[i,]$w)
  }
  
  nodes[[node]] <- temp
}

#Performs a single random walk
randomWalk <- function(nodes, start = 'A', steps = 5){
  step <- 0:steps
  target <- c(start)
  weight <- c(0)
  
  if(!start %in% names(nodes)){
    return(paste(start, "not in list of nodes"))
  }
  
  for(i in 1:steps){
    options<-nodes[[start]]$target
    if(length(options) < 1){
      break
    }
    
    option_weights<-nodes[[start]]$weight
    
    rs<-sample(1:length(options), 1)
    target<-c(target, options[rs])
    weight<-c(weight, option_weights[rs])
    
    start<-options[rs]
    
  }
  
  return(data.frame(step, target, weight))
}

#Testing
randomWalk(nodes, 'K', 20)
randomWalk(nodes, "X", 20)

many_walks <- function(nodes, steps){
  walks<-list()
  
  for(st in names(nodes)){
    walks[[st]]<-randomWalk(nodes, st, steps)
  }
  return(walks)
}

#Not too slow, but depends on the number of nodes you have.
many_walks(nodes, 50)

####################################################### 
# let's generate two networks and merge them into one graph.
g2 <- barabasi.game(50, p=2, directed=F)
g1 <- watts.strogatz.game(1, size=100, nei=5, p=0.05)
g <- graph.union(g1,g2)

# let's remove multi-edges and loops
g <- simplify(g)

# let's see if we have communities here using the 
# Grivan-Newman algorithm
# 1st we calculate the edge betweenness, merges, etc...
ebc <- edge.betweenness.community(g, directed=F)

# Now we have the merges/splits and we need to calculate the modularity
# for each merge for this we'll use a function that for each edge
# removed will create a second graph, check for its membership and use
# that membership to calculate the modularity
mods <- sapply(0:ecount(g), function(i){
  g2 <- delete.edges(g, ebc$removed.edges[seq(length=i)])
  cl <- clusters(g2)$membership
  # March 13, 2014 - compute modularity on the original graph g 
  # (Thank you to Augustin Luna for detecting this typo) and not on the induced one g2. 
  modularity(g,cl)
})

# we can now plot all modularities
plot(mods, pch=20)

# Now, let's color the nodes according to their membership
g2<-delete.edges(g, ebc$removed.edges[seq(length=which.max(mods)-1)])
V(g)$color=clusters(g2)$membership
# Let's choose a layout for the graph
g$layout <- layout.fruchterman.reingold
# plot it
plot(g, vertex.label=NA)

# if we wanted to use the fastgreedy.community agorithm we would do
fc <- fastgreedy.community(g)
com <- igraph::cutat(fc, steps= which.max(fc$modularity)-1)

c <- walktrap.community(g)
b <- cutat(c, steps=2)


V(g)$color <- b
g$layout <- layout.fruchterman.reingold
plot(g, vertex.label=NA)

############################################################
library(igraph)
set.seed(1)
resample <- function(x, ...) x[sample.int(length(x), ...)]
n <- 1000
tm <- matrix(sample(0:1, n^2, prob = c(0.95, 0.05), replace = TRUE), n, n)
tm <- (tm == 1 | t(tm) == 1) * 1
diag(tm) <- 0

start <- 100 # Random walk starting vertex
len <- 10 # Walk length
path <- c(start, rep(NA, len))
for(i in 2:(len + 1)) {
  idx <- tm[path[i - 1], ] != 0
  if(any(idx)) {
    path[i] <- resample(which(idx), 1, prob = tm[path[i - 1], idx])
  } else {
    break # Stoping if we get stuck
  }
}
path

############################################################
# diffusr
# count of nodes
n <- 5
# starting distribution (has to sum to one)
p0    <- as.vector(rmultinom(1, 1, prob=rep(.2, n)))
# adjacency matrix (either normalized or not)
graph <- matrix(abs(rnorm(n*n)), n, n)
# computation of stationary distribution
pt    <- random.walk(p0, graph)
print(t(pt))



