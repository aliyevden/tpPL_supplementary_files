# This example reads a set of datapoints from a file
# It uses path length (TSP), tpPL, and OLO criteria to find orderings.
# It plots the paths and the heat maps all in one graph.

closeAllConnections()  # helps avoid "Error in file(con, "w") : all connections are in use"
library(lattice)      # for the levelplot function
library(latticeExtra) # allow adding plotxy objects
library(MASS)
library(seriation)
library(viridis)
library(doParallel)
registerDoParallel() # for parallel processing in TSP

### function that calculates the tree penalty matrix
penaltymatrix <-function(d,method){ 
  # d = dissimilarity matrix as a "dist" object
  # method = type of linkage: single, complete, average, Ward
  D<-as.matrix(d) # distance matrix as a "matrix" object
  if (method == "ward") {
    method = "ward.D2"
  }
  Tree1 <- hclust(d,method=method)   # tree object from hierarchical clustering
  n  <- length(Tree1$height)+1 # number of leaves in a tree
  a1 <- matrix(cutree(Tree1,1:n),ncol=n,nrow=n,byrow=F,dimnames=NULL) # matrix of cluster identities for different cuts
  h1 <- rev(Tree1$height)# vector of heights of cuts in decreasing order of length n-1
  PM <- matrix(rep(0,n^2),nrow=n) # create a matrix to save tree penalties
  for (i in 1:(n-1)) {
    for(j in (i+1):n) {
      e11 <- c(a1[i,])     # vector of cluster indices for point i
      e12 <- c(a1[j,])     # vector of cluster indices for point j
      lca1 <- max(which(e11==e12))  # least common ancestor node of i and j
      PM[i,j] <- h1[lca1]
    }
  }
  NM <- PM * mean(D) / mean(PM)           # normalized penalty matrix
  return(list(PM=PM,NM=NM))
}

### Reading in data matrix from a file
### This code produces a matrix A of points
A <- read.table("six_normal_clusters_4.txt",header=F)
dim(A)
head(A)
dim <- ncol(A) # number of dimensions
dim

### plot the points, if two dimensional
if (dim == 2) {
  plot(A, pch=20, main="Data points")
}

### define distance matrix and penalty matrix
n <- nrow(A) # number of points
d <- dist(A, method="euclidean") # distance matrix as a "dist" object
D <- as.matrix(d) # distance matrix as a "matrix" object

mylinkage <- "average" ## Choose a hierarchical tree type for tpPL and OLO seriation methods
TP <- penaltymatrix(d,method=mylinkage)$NM # tree penalty matrix (normalized)
tp <- as.dist(TP) # tree penalty matrix as a "dist" object
b <- 0.5 # tuning parameter for tpPL method (0.3-0.8 range recommended)
tppl <- d+b*tp # tree-penalized path length matrix as a "dist" object

### obtain TSP, OLO, and tpPL seriation orderings
rp <- 50 # number of iterations for TSP and tpPL methods
methods <- c("TSP","tpPL","OLO")
ord <- list()
cat("Computing",rp,"TSP orderings\n")
ord$tsp <- get_order(seriate(d,"TSP",control=list(rep=rp,"two_opt"),parallel=T)) 
cat("Computing",rp,"tpPL orderings\n")
ord$tppl <- get_order(seriate(tppl,"TSP",control=list(rep=rp,"two_opt"),parallel=T)) 
cat("Computing OLO ordering\n")
ord$olo <- get_order(seriate(d,method="OLO_average"))

### plot the three orderings
cat("Plotting the orderings\n")
dev.off()
for(i in 1:3){
  score <- criterion(d,ord[[i]],"Path_length")
  score2 <- criterion(tppl,ord[[i]],"Path_length")

  ### start ordering with smaller D values; helps orderings match each other
  BestPath = ord[[i]]
  j = round(n/2)
  if (sum(sum(D[BestPath[1:(j+1)],BestPath[1:(j+1)]])) > sum(sum(D[BestPath[(n-j):n],BestPath[(n-j):n]]))) {
    BestPath = rev(BestPath)  # reverse ordering; not a substantive change
  }
  imagetitle = paste(methods[i],"path length",round(score,4))
  
  ### plot points and line segments to represent the path
  ### starting point is green, ending point is black
  pa <-xyplot(A[BestPath,2]~A[BestPath,1],
              aspect = 1:1,type="p",cex=c(1,rep(0.6,n-2),1),
              pch=20,scales=list(x=list(draw=F),y=list(draw=F)),
              col=rep(c("green","red","black"),c(1,n-2,1)),
              xlab="",ylab="",main=imagetitle
              ) +xyplot(A[BestPath,2]~A[BestPath,1],type="l",
                        cex=0.5, pch=20, col="blue")
  
  print(pa,split=c(1,i,2,3),more=TRUE)
  
  ### create and plot heat map
  f <- function (mm) t(mm)[ ,nrow(mm):1] # rotates matrix 90 deg. clockwise
  ga <-levelplot(f(D[BestPath,BestPath]), xlab=" ", ylab=" ",scales=list(draw=F),
                 col.regions =viridis,
                 colorkey = list(space= "right"))
  
  if (i < 3) {
    print(ga,split=c(2,i,2,3),more=TRUE)
  } else {
    print(ga,split=c(2,i,2,3),more=FALSE)
  }
}
