##' Cross-validation function for quadrupen fitting methods.
##'
##' Function that computes K-fold cross-validated error of a
##' \code{quadrupen} fit, possibly on a grid of
##' \code{lambda1,lambda2}.
##'
##' @param penalty a string for the fitting procedure used for
##' cross-validation. Either \code{"lasso"}, \code{"elastic.net"},
##' \code{"bounded.reg"}, at the moment. Default is
##' \code{elastic.net}.
##'
##' @param x matrix of features, possibly sparsely encoded
##' (experimental). Do NOT include intercept.
##'
##' @param y response vector.
##'
##' @param K integer indicating the number of folds. Default is 10.
##'
##' @param folds list of \code{K} vectors that describes the folds to
##' use for the cross-validation. By default, the folds are randomly
##' sampled with the specified K. The same folds are used for each
##' values of \code{lambda2}.
##'
##' @param lambda2 tunes the \eqn{\ell_2}{l2}-penalty (ridge-like) of
##' the fit. If none is provided, a vector of values is generated and
##' a CV is performed on a grid of \code{lambda2} and \code{lambda1},
##' using the same folds for each \code{lambda2}). Ignored when
##' \code{penalty} equals \code{"lasso"}. CV is only performed on
##' \code{lambda2} when the \code{ridge} penalty is used.
##'
##' @param verbose logical; indicates if the progression (the current
##' lambda2) should be displayed. Default is \code{TRUE}.
##'
##' @param mc.cores the number of cores to use. The default uses all
##' the cores available.
##'
##' @param ... additional parameters to overwrite the defaults of the
##' fitting procedure identified by the \code{'penalty'} argument. See
##' the corresponding documentation (\code{\link{lasso}},
##' \code{\link{elastic.net}} or \code{\link{bounded.reg}}).
##'
##' @note If the user runs the fitting method with option
##' \code{'bulletproof'} set to \code{FALSE}, the algorithm may stop
##' at an early stage of the path. Early stops are handled internally,
##' in order to provide results on the same grid of penalty tuned by
##' \eqn{\lambda_1}{lambda1}.  This is done by means of \code{NA}
##' values, so as mean and standard error are consistently
##' evaluated. If, while cross-validating, the procedure experiences
##' too much early stoppings, a warning is sent to the user, in which
##' case you should reconsider the grid of \code{lambda1} used for the
##' cross-validation.  If \code{bulletproof} is \code{TRUE} (the
##' default), there is nothing to worry about, except a possible slow
##' down when any switching to the proximal algorithm is required.
##'
##' @return An object of class "cvpen" for which a \code{plot} method
##' is available.
##'
##' @seealso \code{\linkS4class{quadrupen}}, \code{\link{plot.cvpen}}
##' and \code{\linkS4class{cvpen}}.
##'
##' @examples \dontrun{
##' ## Simulating multivariate Gaussian with blockwise correlation
##' ## and piecewise constant vector of parameters
##' beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
##' cor  <- 0.75
##' Soo  <- toeplitz(cor^(0:(25-1))) ## Toeplitz correlation for irrelevant variable
##' Sww  <- matrix(cor,10,10) ## bloc correlation between active variables
##' Sigma <- bdiag(Soo,Sww,Soo,Sww,Soo) + 0.1
##' diag(Sigma) <- 1
##' n <- 100
##' x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
##' y <- 10 + x %*% beta + rnorm(n,0,10)
##'
##' ## Use fewer lambda1 values by overwritting the default parameters
##' ## and cross-validate over the sequences lambda1 and lambda2
##' cv.grid <- crossval(x,y, lambda2=10^seq(2,-2,len=50), nlambda1=50)
##' ## Rerun simple cross-validation with the appropriate lambda2
##' cv.10K <- crossval(x,y, lambda2=slot(cv.double, "lambda2.min"))
##' ## Try leave one out also
##' cv.loo <- crossval(x,y, K=n, lambda2=slot(cv.double, "lambda2.min"))
##'
##' plot(cv.grid)
##' plot(cv.10K)
##' plot(cv.loo)
##'
##' ## Performance for selection purpose
##' beta.min.10K <- slot(cv.10K, "beta.min")
##' beta.min.loo <- slot(cv.loo, "beta.min")
##'
##' cat("\nFalse positives with the minimal 10-CV choice: ", sum(sign(beta) != sign(beta.min.10K)))
##' cat("\nFalse positives with the minimal LOO-CV choice: ", sum(sign(beta) != sign(beta.min.loo)))
##' }
##'
##' @keywords models, regression
##' @name crossval
##' @aliases crossval
##' @rdname crossval
##'
##' @export
crossval <- function(x,
                     y,
                     penalty  = c("elastic.net", "lasso", "bounded.reg", "ridge"),
                     K        = 10,
                     folds    = split(sample(1:nrow(x)), rep(1:K, length=nrow(x))),
                     lambda2  = switch(penalty, "ridge" = 10^seq(-2,2,len=100), 10^seq(-2,0,len=20) ),
                     verbose  = TRUE,
                     mc.cores = detectCores(),
                     ...) {

  ## =============================================================
  ## INITIALIZATION & PARAMETERS RECOVERY
  if (Sys.info()[['sysname']] == "Windows") {
    warning("\nWindows does not support fork, enforcing mc.cores to '1'.")
    mc.cores <- 1
  }
  penalty <- match.arg(penalty)
  fit.func <- switch(penalty,
                     "elastic.net" = elastic.net,
                     "lasso"       = lasso      ,
                     "bounded.reg" = bounded.reg,
                     "ridge"       = ridge)
  get.lambda1 <- switch(penalty,
                        "elastic.net" = get.lambda1.l1,
                        "lasso"       = get.lambda1.l1,
                        "bounded.reg" = get.lambda1.li,
                        "ridge"       = NULL)

  if (penalty == "lasso") {lambda2 <- NULL}


  user <- list(...)
  defs <- default.args(penalty,nrow(x)-max(sapply(folds,length)),ncol(x),user)
  args <- modifyList(defs, user)
  ## Compute a grid of lambda1 (the same for each fold)
  if (is.null(args$lambda1) & penalty != "ridge") {
    input <- standardize(x,y,args$intercept,args$normalize,args$penscale)
    args$lambda1 <- get.lambda1(input$xty,args$nlambda1,args$min.ratio)
    rm(input)
  }

  if (penalty=="ridge") {
    major.lambda <- lambda2
  } else {
    major.lambda <- args$lambda1
  }

  ## =============================================================
  if (length(lambda2) > 1 & penalty != "ridge") {
    ## CROSS-VALIDATION WORK ON A GRID
    if (verbose){
      cat("\nCROSS-VALIDATION FOR ",penalty," REGULARIZER \n\n")
      cat(length(folds),"-fold CV on the lambda1 grid for each lambda2\n", sep="")
    }
    cv <- sapply(1:length(lambda2), function(i) {
      if(verbose){
        cat(round(lambda2[i],3),"\t")
        if (i %% 5 == 0) {cat("\n")}
      }
      simple.cv(folds, x, y, fit.func, args, lambda2[i], mc.cores)
    }, simplify=FALSE)
    if(verbose){cat("\n")}

    ## Recovering the best lambda1 and lambda2
    lambda1.cv <- sapply(1:length(lambda2), function(j) {
      cv.min <- min(cv[[j]]$mean)
      lambda1.min   <- max(major.lambda[cv[[j]]$mean <= cv.min], na.rm=TRUE)
      lambda1.1se   <- max(major.lambda[cv[[j]]$mean <(cv[[j]]$mean+cv[[j]]$serr+1e-5)[match(lambda1.min,major.lambda)]], na.rm=TRUE)
      return(c(cv.min, lambda1.min, lambda1.1se))
    })
    ind.lb2.min <- which.min(lambda1.cv[1,])
    lambda2.min <- lambda2[ind.lb2.min]
    lambda2.1se <- NULL
    lambda1.min <- lambda1.cv[2, ind.lb2.min]
    lambda1.1se <- lambda1.cv[3, ind.lb2.min]

    ## formatting cv.error for ggplot
    cv <- data.frame(do.call(rbind,cv),
                     lambda1=rep(major.lambda, length(lambda2)),
                     lambda2=rep(lambda2,rep(length(major.lambda),length(cv))))

    ## Apply the fitting procedure with these best lambda2 parameter
    args$lambda2 <- lambda2.min
    args$checkargs <- FALSE ## enforcing checkargs to FALSE to save some time
    best.fit <- do.call(fit.func, c(list(x=x,y=y),args))

    ind.max <- nrow(best.fit@coefficients)
    ind.min <- min(match(lambda1.min, major.lambda),ind.max)
    ind.1se <- min(match(lambda1.1se, major.lambda),ind.max)

  } else {
    ## SIMPLE CROSS-VALIDATION WORK
    if (penalty == "ridge") {
      if (verbose) {
        cat("\nCROSS-VALIDATION FOR ",penalty," REGULARIZER \n\n")
        cat(length(folds),"-fold CV on the lambda2 grid.\n", sep="")
      }
      cv <- simple.cv(folds, x, y, fit.func, args, lambda2, mc.cores)

      ## Recovering the best lambda1 and lambda2
      lambda1.min <- 0
      lambda1.1se <- 0
      lambda2.min <- max(major.lambda[cv$mean <= min(cv$mean)], na.rm=TRUE)
      lambda2.1se <- max(major.lambda[cv$mean <(cv$mean+cv$serr+1e-5)[match(lambda2.min,major.lambda)]], na.rm=TRUE)

      ## Apply the fitting procedure with these best lambda2 parameter
      args$lambda2 <- lambda2.min
      args$checkargs <- FALSE ## enforcing checkargs to FALSE to save some time
      best.fit <- do.call(fit.func, c(list(x=x,y=y),args))

      ind.max <- nrow(best.fit@coefficients)
      ind.min <- min(match(lambda2.min, major.lambda),ind.max)
      ind.1se <- min(match(lambda2.1se, major.lambda),ind.max)

    } else {
      if (verbose) {
        cat("\nCROSS-VALIDATION FOR ",penalty," REGULARIZER \n\n")
        cat(length(folds),"-fold CV on the lambda1 grid, lambda2 is fixed.\n", sep="")
      }
      cv <- simple.cv(folds, x, y, fit.func, args, lambda2, mc.cores)

      ## Recovering the best lambda1 and lambda2
      lambda1.min <- max(major.lambda[cv$mean <= min(cv$mean)], na.rm=TRUE)
      lambda1.1se <- max(major.lambda[cv$mean <(cv$mean+cv$serr+1e-5)[match(lambda1.min,major.lambda)]], na.rm=TRUE)
      lambda2.min <- lambda2
      lambda2.1se <- lambda2

      ## Apply the fitting procedure with these best lambda2 parameter
      args$lambda2 <- lambda2.min
      args$checkargs <- FALSE ## enforcing checkargs to FALSE to save some time
      best.fit <- do.call(fit.func, c(list(x=x,y=y),args))

      ind.max <- nrow(best.fit@coefficients)
      ind.min <- min(match(lambda1.min, major.lambda),ind.max)
      ind.1se <- min(match(lambda1.1se, major.lambda),ind.max)

    }
    cv <- data.frame(cv, lambda=major.lambda)
  }

  ## Finally recover the CV choice (minimum and 1-se rule)
  beta.min <- best.fit@coefficients[ind.min,]
  beta.1se <- best.fit@coefficients[ind.1se,]

  return(new("cvpen",
             lambda1     = args$lambda1,
             lambda1.min = lambda1.min,
             lambda1.1se = lambda1.1se,
             lambda2     = lambda2,
             lambda2.min = lambda2.min,
             lambda2.1se = lambda2.1se,
             cv.error    = cv,
             folds       = folds,
             beta.min    = beta.min,
             beta.1se    = beta.1se))
}

simple.cv <- function(folds, x, y, fit.func, args, lambda2, mc.cores) {

  K <- length(folds)
  n <- length(y)

  ## overwrite irrelevant arguments
  args$control$verbose <- 0
  args$lambda2         <- lambda2
  args$checkargs       <- FALSE
  if (!is.null(args$max.feat))
    args$max.feat        <- ncol(x)

  ## Multicore approach
  one.fold <- function(k) {
    omit <- folds[[k]]
    fit <- do.call(fit.func, c(list(x=x[-omit, ], y=y[-omit]), args))
    fold.err <- sweep(matrix(predict(fit,matrix(x[omit,], nrow=length(omit))), nrow=length(omit)), 1L, y[omit], check.margin = FALSE)^2
    if (ncol(fold.err) < length(args$lambda1) & length(args$lambda2 == 1)) {
      NAs <- length(args$lambda1)-ncol(fold.err)
      fold.err <- cbind(fold.err, matrix(NA,nrow(fold.err),NAs))
    }
    return(fold.err)
  }
  ## turn a list to matrix
  err  <- do.call(rbind,mclapply(1:K, one.fold, mc.cores=mc.cores,
                                 mc.preschedule=ifelse(K > 10,TRUE,FALSE)))
  ## efficient computation of means and the standard error
  mean <- colMeans(err, na.rm=TRUE)
  if (any(is.nan(mean))) {
    warning("\nThere have been a lot of early stops along the path: I keep on running, but you really should reconsider 'min.ratio' regarding the n<<p setting.")
  }
  mean[is.nan(mean)] <- NA
  serr <- sqrt((colSums(sweep(err, 2L, mean, check.margin = FALSE)^2,na.rm=TRUE)/(n-1)) /K)

  return(data.frame(mean=mean, serr=serr))
}
