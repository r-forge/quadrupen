context("Consistency of the Elastic-net solution path")

test_that("testing warm start", {

  require(quadrupen)
  
  get.coef <- function(x,y) {
    lambda1 <- .25
    
    enet.ref <- elastic.net(x,y,lambda1=lambda1, control=list(timer=TRUE))
    enet.ref.bot <- elastic.net(x,y,lambda1=lambda1*2)

    enet.cg  <- elastic.net(x,y,lambda1=lambda1, control=list(timer=TRUE,usechol=FALSE,threshold=1e-3))
    enet.cg.warm <- elastic.net(x,y,lambda1=lambda1, beta0 = enet.ref.bot@coefficients, control=list(timer=TRUE,usechol=FALSE,threshold=1e-3))
    
    cat("\nTimings with warm-restart along the path")
    cat("\n\tfrom stratch (cholesky): ",enet.ref@monitoring$internal.timer)
    cat("\n\tfrom stratch (conjugate-gradient): ",enet.cg@monitoring$internal.timer)
    cat("\n\tCG starting from sparser solution: ",enet.cg.warm@monitoring$internal.timer)
        
    return(list(
      coef.ref=as.matrix(enet.ref@coefficients),
      coef.cg.warm=as.matrix(enet.cg.warm@coefficients),
      coef.cg =as.matrix(enet.cg@coefficients)))
  }
  
  ## PROSTATE DATA SET
  prostate <- read.table("http://www-stat.stanford.edu/~tibs/ElemStatLearn/datasets/prostate.data")
  x <- as.matrix(prostate[,1:8])
  y <- prostate[,9]

  ## Run the tests...
  out <-get.coef(x,y)
  expect_that(out$coef.cg.warm,is_equivalent_to(out$coef.ref))
  expect_that(out$coef.cg ,is_equivalent_to(out$coef.ref))

  ## RANDOM DATA
  seed <- sample(1:10000,1)
  ## cat(" #seed=",seed)
  set.seed(seed)

  beta <- rep(rep(c(0,1,0,-1,0), c(25,10,25,10,25)),5)
  n <- 300
  p <- length(beta)

  mu <- 3 # intercept
  sigma <- 30 # huge noise
  Sigma <- matrix(0.95,p,p) # huge correlation
  diag(Sigma) <- 1

  x <- as.matrix(matrix(rnorm(p*n),n,p) %*% chol(Sigma))
  y <- 10 + x %*% beta + rnorm(n,0,10)

  ## Run the tests...
  out <-get.coef(x,y)
  expect_that(out$coef.cg.warm,is_equivalent_to(out$coef.ref))
  expect_that(out$coef.cg ,is_equivalent_to(out$coef.ref))
  
})

