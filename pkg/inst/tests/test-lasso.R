context("Consistency of the Lasso solution path")

test_that("Consistency between 'quadrupen' and 'lars' packages", {

  require(lars)

  get.coef <- function(x,y,intercept,normalize) {
      lasso.larsen <- lars(x,y,intercept=intercept,normalize=normalize)
      iols <- nrow(lasso.larsen$beta) ## remove last entry corresponding to the OLS estimator
      lambda1 <-  lasso.larsen$lambda ## usde the lars lambda grid
      lasso.quadru <- elastic.net(x,y, intercept=intercept, normalize=normalize,
                                  lambda1=lambda1, lambda2=0, control=list(method="quadra"))
      return(list(coef.quad=as.matrix(lasso.quadru@coefficients),
                  coef.lars=lasso.larsen$beta[-iols, ]))
  }

  ## PROSTATE DATA SET
  prostate <- read.table("http://www-stat.stanford.edu/~tibs/ElemStatLearn/datasets/prostate.data")
  x <- as.matrix(prostate[,1:8])
  y <- prostate[,9]

  ## Run the tests...
  with.intercept <-get.coef(x,y,TRUE,TRUE)
  expect_that(with.intercept$coef.quad,
              is_equivalent_to(with.intercept$coef.lars))

  with.intercept.unnormalized <-get.coef(x,y,TRUE,FALSE)
  expect_that(with.intercept.unnormalized$coef.quad,
              is_equivalent_to(with.intercept.unnormalized$coef.lars))

  without.intercept <-get.coef(x,y,FALSE,TRUE)
  expect_that(without.intercept$coef.quad,
              is_equivalent_to(without.intercept$coef.lars))

  without.intercept.unnormalized <-get.coef(x,y,FALSE,FALSE)
  expect_that(without.intercept.unnormalized$coef.quad,
              is_equivalent_to(without.intercept.unnormalized$coef.lars))

  ## RANDOM DATA
  seed <- sample(1:10000,1)
  ## cat("\n#seed=",seed)
  set.seed(seed)

  beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
  n <- 100
  p <- length(beta)

  mu <- 3 # intercept
  sigma <- 30 # huge noise
  Sigma <- matrix(0.95,p,p) # huge correlation
  diag(Sigma) <- 1

  x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
  y <- 10 + x %*% beta + rnorm(n,0,10)

  ## Run the tests...
  with.intercept <-get.coef(x,y,TRUE,TRUE)
  expect_that(with.intercept$coef.quad,
              is_equivalent_to(with.intercept$coef.lars))

  with.intercept.unnormalized <-get.coef(x,y,TRUE,FALSE)
  expect_that(with.intercept.unnormalized$coef.quad,
              is_equivalent_to(with.intercept.unnormalized$coef.lars))

  without.intercept <-get.coef(x,y,FALSE,TRUE)
  expect_that(without.intercept$coef.quad,
              is_equivalent_to(without.intercept$coef.lars))

  without.intercept.unnormalized <-get.coef(x,y,FALSE,FALSE)
  expect_that(without.intercept.unnormalized$coef.quad,
              is_equivalent_to(without.intercept.unnormalized$coef.lars))

})

test_that("Consistency between 'quadrupen' and  'glmnet'", {

  require(glmnet)

  ## SECOND CHECK: compare to glmnet with prescaling of x
  x <- matrix(rnorm(100*50),100,50)
  y <- rnorm(100)
  y <- y-mean(y)
  n <- nrow(x)
  p <- ncol(x)

  ## If thresh is set to the default, the test won't pass!!!
  ## This is beacause coordinate descent is fast yet not extremely accurate
  lasso.glmn <- glmnet(x,y, lambda.min.ratio=1e-2, thresh=1e-20)
  lasso.quad <- elastic.net(x,y, lambda1=lasso.glmn$lambda*sqrt(n), lambda2=0)

  expect_that(as.matrix(lasso.quad@coefficients), is_equivalent_to(as.matrix(t(lasso.glmn$beta))))
})
