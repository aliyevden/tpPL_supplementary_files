##### All Datasets, Seriation Article #####
### Basic Diagnostic Datasets
library(lattice)
library(MASS)

n1 <- 100 # number of points per block
w <- 0.9    # width of rectangle
h <- 1/3*w  # height of rectangle


## SET1: 3 rectangles, below each other on longer sides
rec1 <- matrix(c(runif(n1,0, w),runif(n1,0,0+h)),byrow=F,ncol=2)
rec2 <- matrix(c(runif(n1,0, w), runif(n1, 2*h,3*h)),byrow=F,ncol=2) 
rec3 <- matrix(c(runif(n1,0, w), runif(n1,-2*h,-h)),byrow=F,ncol=2)
A1 <- rbind(rec1,rec2,rec3) 
A1 <- A1*1/max(dist(A1))
rand <- sample(nrow(A1))
A1 <- A1[rand,]
n <- n1*3
eqscplot(A1,main="Basic Diagnostic 1",pch=20)


## SET2: 3 rectangles, 2 next to each other and one on their top
rec1 <- matrix(c(runif(n1,0, h),runif(n1,0,0+w)),byrow=F,ncol=2)
rec2 <- matrix(c(runif(n1,2*h, 3*h), runif(n1, 0, 0+w)),byrow=F,ncol=2) 
rec3 <- matrix(c(runif(n1,0,0+w), runif(n1, w+h, w+2*h)),byrow=F,ncol=2) 
A2 <- rbind(rec1,rec2,rec3) 
A2 <- A2*1/max(dist(A2))
rand <- sample(nrow(A2))
A2 <- A2[rand,]
n <- n1*3
eqscplot(A2,ratio=1,main="Basic Diagnostic 2",pch=20)


# ## SET3: 3 rectangles, below each other on shorter sides
rec1 <- matrix(c(runif(n1,0, w),runif(n1,0,0+h)),byrow=F,ncol=2)
rec2 <- matrix(c(runif(n1,w/3, w+w/3), runif(n1, -2*h,-h)),byrow=F,ncol=2) 
rec3 <- matrix(c(runif(n1,-w/3, 0), runif(n1,-h-w, -h)),byrow=F,ncol=2)
A3 <- rbind(rec1,rec2,rec3) 
A3 <- A3*1/max(dist(A3))
rand <- sample(nrow(A3))
A3 <- A3[rand,]
n <- n1*3
eqscplot(A3,ratio=1,main="Basic Diagnostic 3",pch=20)


## SET4: 3 rectangles, Greek letter Pi
rec1 <- matrix(c(runif(n1,0, h),runif(n1,0,0+w)),byrow=F,ncol=2)
rec2 <- matrix(c(runif(n1,2*h, 3*h), runif(n1, 0, 0+w)),byrow=F,ncol=2) 
rec3 <- matrix(c(runif(n1,0,0+w), runif(n1, w, w+h)),byrow=F,ncol=2) 
A4 <- rbind(rec1,rec2,rec3) 
A4 <- A4*1/max(dist(A4))
rand <- sample(nrow(A4))
A4 <- A4[rand,]
n <- n1*3
eqscplot(A4,ratio=1,main="Basic Diagnostic 4",pch=20)


## SET5: Ring Data set
n1<-100
n <- n1*3   # number of points
r1 =w/2
r2 = w*sqrt(1/pi+1/4)
A5 <- matrix(rep(0,2*n),ncol=2)
to <-0
while (to < n) { 
  x <-runif(1,-r2,r2)
  y <- runif(1,-r2,r2)
  br <- sqrt(x^2+y^2)
  if (br<= r2 & br>=r1) {
    to <- to+1
    A5[to, 1] <-x
    A5[to,2] <-y}
}
A5 <- A5*1/max(dist(A5))
rand <- sample(nrow(A5))
A5 <- A5[rand,]
eqscplot(A5,ratio=1,main="Basic Diagnostic 5",pch=20) 


## SET6: 5 multivariate samples in a square
R <- 5 # length of a side
n2 <-30 # number of points per center
n <- n2*5 # total number of points
sigma <- 0.3 # variance
v1 <- c(0,0)
v2 <- c(1,0)
v3 <- c(0,1)
v4 <- c(1,1)
OO <- c(1/2,1/2)
VV <-matrix(rbind(v1,v2,v3,v4,OO),ncol=2)
VV <- R*VV
sigma <-diag(rep(sigma,2))
obs1 <-mvrnorm(n2,VV[1,],sigma)
obs2 <- mvrnorm(n2,VV[2,],sigma)
obs3 <- mvrnorm(n2,VV[3,],sigma)
obs4 <- mvrnorm(n2,VV[4,],sigma)
obs5 <- mvrnorm(n2,VV[5,],sigma)
A6 <- matrix(rbind(obs1,obs2,obs3,obs4,obs5),ncol=2)
A6 <- A6*1/max(dist(A6))
rand <- sample(nrow(A6))
A6 <- A6[rand,]
plot(A6, main="Basic Diagnostic 6",pch=20)


### Advanced dataset 1
## OFFLINE POINTS DATA
# Generating n1+n2=N points in 2D
n1 <- 48   ## number of points on a line
line_points <- matrix(c(runif(n1,1,4),rep(2,n1)),byrow=F,ncol=2)
n2 <-6   ## number of random points off line
n <- n1+n2  ## total number of points
rand_points <- matrix(c(runif(n2,1,4),runif(n2,0,4)),byrow=F,ncol=2)
A1 <- rbind(line_points,rand_points)
A1 <- A1*1/max(dist(A1))
rand <- sample(nrow(A1))
A1 <- A1[rand,]
plot(A1,main="Advanced Diagnostic 1",pch=20)


### Advanced dataset 2
## Smiling man
## eye 1
n1 <- 40
r1 <- sqrt(runif(n1,0,0.1))
theta <- runif(n1,0,2*pi)
x1 <- r1*cos(theta)+2
y1 <- r1*sin(theta)+1
## eye2
n2 <- n1
r2 <- sqrt(runif(n2,0,0.1))
x2 <- r2*cos(theta)-2
y2 <- r2*sin(theta)+1
## face
n3 <- 180
r3 <- 3.8
xx3 <- rnorm(n3,0,0.1)
yy3 <- rnorm(n3,0,0.1)
beta <- runif(n3,0,2*pi)
x3 <- r3*cos(beta)+xx3
y3 <- r3*sin(beta)+yy3

## smile
n4 <-80
r4 <- 1.5
gamma <- runif(n4,-pi,0)
xx4 <- rnorm(n4,0,0.1)
yy4 <- rnorm(n4,0,0.1)
x4 <- r4*cos(gamma)+xx4
y4 <- r4*sin(gamma)+yy4-0.7
n <- n1+n2+n3+n4
A2 <- matrix(c(x1,x2,x3,x4,y1,y2,y3,y4),nrow=n)
A2 <- A2*1/max(dist(A2))
rand <- sample(nrow(A2))
A2 <- A2[rand,]
plot(A2,main="Advanced Diagnostic 2",pch=20)


### Advanced dataset 3
# 3 Circles data set
n <- 150    # number of points
r = sample(c(1,2,3),n,replace=T,prob=c(1/6,1/3,1/2))      # radii lengths
theta = runif(n,0,2*pi) # random angles between 0 to 2*pi
xx <-   rnorm(n,0,0) #rnorm(n,0,0.1)
yy <-   rnorm(n,0,0) #rnorm(n,0,0.05)
x = r * cos(theta)+xx
y = r * sin(theta)+yy
A3 <- matrix(cbind(x,y),nrow=length(x))
A3<- A3*1/max(dist(A3))
rand <- sample(nrow(A3))
A3 <- A3[rand,]
plot(A3,main="Advanced Diagnostic 3",pch=20)


### Advanced dataset 4
# pinwheel data
n <- 270    # number of points
theta = (2*pi/3)*floor(runif(n)*3) + (pi/3)*runif(n)
rr = sample.int(n,n)
x = rr * cos(theta)
y = rr * sin(theta)
A4 <- matrix(cbind(x,y),nrow=length(x))
A4 <- A4*1/max(dist(A4))
rand <- sample(nrow(A4))
A4 <- A4[rand,]
plot(A4,main="Advanced Diagnostic 4",pch=20)


### Advanced dataset 5
### A cross through dense cluster + a few outliers
n1 <-33 # cluster points
x1 <- runif(n1,-0.1,0.1)
y1 <- runif(n1,-0.1,0.1)
n2 <- 10 # noise points
x2 <- runif(n2,-1,1)
y2 <- runif(n2,-1,1)
n3 <-31 # horizontal line
x3 <- runif(n3,-1,1)  
y3 <- runif(n3,-0.1,0.1)
x4 <- runif(n3,-0.1,0.1) # vertical line
y4 <- runif(n3,-1,1)
n <-  n1+n2+2*n3
A5 <- matrix(c(x1,x2,x3,x4,y1,y2,y3,y4),nrow=n)
A5 <- A5*1/max(dist(A5))
rand <- sample(nrow(A5))
A5 <- A5[rand,]
plot(A5,main="Advanced Diagnostic 5",pch=20)


### Advanced dataset 6
## Letter I
n1<-100
w <-1
h <- w/3
rec1 <- matrix(c(runif(n1,0, w),runif(n1,0,0+h)),byrow=F,ncol=2)
rec2 <- matrix(c(runif(n1,0, w), runif(n1, w+h, w+2*h)),byrow=F,ncol=2) 
rec3 <- matrix(c(runif(n1,h, 2*h), runif(n1, h, w+h)),byrow=F,ncol=2) 
A6 <- rbind(rec1,rec2,rec3) 
A6 <- A6*1/max(dist(A6))
rand <- sample(nrow(A6))
A6 <- A6[rand,]
n <- n1*3
plot(A6,main="Advanced Diagnostic 6",pch=20)


########## 6 Normal Clusters
### cen multivariate clusters
dim <- 2 # number of dimensions
cen <- 6 # number of clusters

s <- matrix(0.004,cen)
n1 <- matrix(300/cen,cen)
mu <- matrix(runif(dim*cen),ncol=dim) # matrix of centers

C <- list()

for (i in 1:cen) {
  C[[i]] <- mvrnorm(n1[i],mu[i,], diag(rep(s[i],dim),dim,dim))
}

A <- C[[1]]
for (i in 2:cen) {
  A <- rbind(A,C[[i]])
}
plot(A, main="Six Normal Clusters",pch=20)


###### 2 Parallel Lines ######
n <- 20 # number of points per line
x = c(rep(0,n), rep(1,n))
y = rep(seq(0,5,length.out=n),2)
A <- matrix(c(x,y),ncol=2)
plot(A,main="Two parallel lines",pch=20)


####### 5 by two clusters ##############
xcenters = c(0,0,0,0,0,1,1,1,1,1)
ycenters = c(0,1,2,3,4,0,1,2,3,4)

n = 30    # points per cluster
dim = 2
s = 0.05

C <- list()

mu <- matrix(rep(0, dim*length(xcenters)), nrow=length(xcenters))

for (i in 1:length(xcenters)) {
  mu[i,1] = xcenters[i]
  mu[i,2] = ycenters[i]
  C[[i]] <- mvrnorm(n,mu[i,], diag(rep(s,dim),dim,dim))
}

A <- C[[1]]
for (i in 2:length(xcenters)) {
  A <- rbind(A,C[[i]])
}

plot(A, main="Five by two clusters",pch=20)


#### single bivariate normal cluster
n <-300 # number of points per center
sigma <- 1 # variance
v1 <- c(0,0)
VV <-matrix(rbind(v1),ncol=2)
sigma <-diag(rep(sigma,2))
A8 <-matrix(mvrnorm(n,VV[1,],sigma),ncol=2)
A8 <- A8*1/max(dist(A8))
rand <- sample(nrow(A8))
A8 <- A8[rand,]
eqscplot(A8,ratio=1,main="Single normal cluster",pch=20)


### One normal cluster 2 by 1
dim <- 2 # number of dimensions
cen <- 1 # number of clusters

s <- matrix(1,cen)
n1 <- matrix(300/cen,cen)
mu <- matrix(0,1,2)

C <- list()

for (i in 1:cen) {
  C[[i]] <- mvrnorm(n1[i],mu[i,], diag(c(1,4),dim,dim))
}

A <- C[[1]]
if (cen > 1) {
  for (i in 2:cen) {
    A <- rbind(A,C[[i]])
  }
}

plot(A, main="One Normal Cluster 2 by 1",pch=20)


### Hypertriangle Normal Clusters
## Generate dim dimensional normal clusters at the vertices of equilateral hyper triangles
dim <- 2 # number of dimensions
cen <- dim+1 # # of clusters= # of vertices
N <- 300 # total # of points
n <- ceiling(N/cen) # # of points per cluster
mu <- matrix(rep(0, dim*cen), nrow=cen)
mu <- diag(ncol=dim, nrow=cen) # matrix of cluster centers
mu[cen,] <- rep((1+sqrt(1+dim))/dim,dim)
mu <- 1/sqrt(2)*mu # scaled matrix of cluster centers

s <- 0.05 # std
C <- list()

for (i in 1:cen) {
  C[[i]] <- mvrnorm(n,mu[i,], diag(rep(s,dim),dim,dim))
}

A <- C[[1]]
for (i in 2:cen) {
  A <- rbind(A,C[[i]])
}
plot(A,main=paste("Hypertriangle Normal Clusters",dim,"Dim"),pch=20)
