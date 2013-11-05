library(quadrupen)

beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
cor <- 0.75
Soo <- toeplitz(cor^(0:(25-1))) ## Toeplitz correlation for irrelevant variables
Sww  <- matrix(cor,10,10) ## bloc correlation between active variables
Sigma <- bdiag(Soo,Sww,Soo,Sww,Soo)
diag(Sigma) <- 1
n <- 50
x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
y <- 10 + x %*% beta + rnorm(n,0,10)

out <- ridge(x,y)
plot(out)

p <- ncol(x)
C <- bandSparse(p,k=0:1,diagonals=list(rep(1,p),rep(-1,p-1)))
L <- diag(rep(10,p)) %*% t(C) %*% C %*% diag(rep(10,p))

plot(ridge(x,y, struct=L, lambda.max=1000))
