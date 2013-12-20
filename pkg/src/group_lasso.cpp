/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */

#include "group_lasso.h"

using namespace Rcpp;
using namespace arma;

// CONTRUCTOR
GROUP_LASSO::GROUP_LASSO(SEXP X        , // matrix of features
			 SEXP Y        , // vector of response	    
			 SEXP TYPE     , // -1="l1/linf", 2="l1/l2" group-Lasso
			 SEXP PK       , // successive groups length
			 SEXP OMEGA    , // structuring matrix
			 SEXP PENSCALE ,
			 SEXP LAMBDA2  ,
			 SEXP INTERCEPT,
			 SEXP NORMALIZE,
			 SEXP WEIGHTS  ,
			 SEXP NAIVE    ,
			 SEXP EPS      ,
			 SEXP MAXITER  ,
			 SEXP MAXFEAT  ,
			 SEXP FUN_OPTIM,
			 SEXP VERBOSE  ,
			 SEXP SPARSE   ) :

  y         (as<vec>    (Y)        ), // responses
  pk        (as<vec>    (PK)       ), // group numbers
  weights   (as<vec>    (WEIGHTS)  ), // observation weights (not use at the moment)
  penscale  (as<vec>    (PENSCALE) ), // penalty scales for l1-norm
  lambda2   (as<double> (LAMBDA2)  ), // the smooth (ridge) penality
  intercept (as<bool>   (INTERCEPT)), // boolean for intercept mode
  normalize (as<bool>   (NORMALIZE)), // boolean for standardizing the predictor
  naive     (as<bool>   (NAIVE)    ), // naive elastic-net or not
  eps       (as<double> (EPS)      ), // precision required
  fun       (as<int>    (FUN_OPTIM)), // solver (0=quadra, 1=pathwise, 2=fista)
  verbose   (as<int>    (VERBOSE)  ), // int for verbose mode (0/1/2)
  sparse    (as<bool>   (SPARSE)   ), // boolean for sparse mode
  max_iter  (as<int>    (MAXITER)  ), // max # of iterates of the active set
  max_feat  (as<int>    (MAXFEAT)  )  // max # of variables activated

{

  // construct the structuring matrix associated to the OMEGA SEXP
  omega = get_struct(OMEGA, lambda2, penscale);

  // handling sparse / dense predictors and standardization work
  if ( sparse ) {
    sp_x  =  as<sp_mat> (X) ;
    standardize(sp_x, y, intercept, normalize, penscale, xty, normx, normy, xbar, ybar) ;
    sp_xt = sp_x.t()        ;
  } else {
    x  =  as<mat>       (X) ;
    standardize(x, y, intercept, normalize, penscale, xty, normx, normy, xbar, ybar) ;
    x  = x - sqrt(weights) * trans(xbar) ;
    xt = x.t()              ;
  }
  
  
  
}

SEXP group_lasso(SEXP BETA0    ,
		 SEXP X        , // matrix of features				    
		 SEXP Y        , // vector of response				    
		 SEXP TYPE     , // -1="l1/linf", 2="l1/l2" group-Lasso
		 SEXP PK       , // successive groups length
		 SEXP OMEGA    , // structuring matrix
		 SEXP LAMBDA1  , 
		 SEXP NLAMBDA  ,
		 SEXP MIN_RATIO,
		 SEXP PENSCALE ,
		 SEXP LAMBDA2  ,
		 SEXP INTERCEPT,
		 SEXP NORMALIZE,
		 SEXP WEIGHTS  ,
		 SEXP NAIVE    ,
		 SEXP EPS      ,
		 SEXP MAXITER  ,
		 SEXP MAXFEAT  ,
		 SEXP FUN_OPTIM,
		 SEXP VERBOSE  ,
		 SEXP SPARSE   ) {

  // disable messages being printed to the err2 stream
  std::ostream nullstream(0);
  set_stream_err2(nullstream);
  
  // ==============================================================
  // INSTANTIATE THE GROUP LASSO PATH
  GROUP_LASSO grp_lasso(X, Y, TYPE, PK, OMEGA, PENSCALE, LAMBDA2, INTERCEPT, NORMALIZE, WEIGHTS, NAIVE, EPS, MAXITER, MAXFEAT, FUN_OPTIM, VERBOSE, SPARSE);
  
//   // ==============================================================
//   // GET LAMBDA VALUES
//   grp_lasso.get_lambda(LAMBDA, NLAMBDA, LAMBDAMIN, LAMBDAMAX) ;

//   // ==============================================================
//   // COMPUTE THE PATH OF SOLUTIONS
//   grp_lasso.LetsRoll();

   return List::create(Named("coefficients") = grp_lasso.get_coef(),
		       Named("mu")           = grp_lasso.get_mu()    ,
		       Named("normx")        = grp_lasso.get_normx() ,
		       Named("lambda2")      = grp_lasso.get_lambda(),
		       Named("df")           = grp_lasso.get_df()    );

}

// void GROUP_LASSO::LetsRoll() {

//   arma::mat cinv = inv(trimatu(c)) ; // inverting the Cholesky decomp. of the structuring matrix

//   // // SVD DECOMPOSITION OF ( X * C^-1)
//   // arma::vec eta ; // eigen value of X cinv
//   // arma::mat U   ; // left singular vectors of X
//   // arma::mat V   ; // right singular vectors of X
//   // svd_econ(U, eta, V, x*cinv) ;
  
//   // arma::mat cinvV = cinv * V ;
//   // arma::mat Uty = trans(U) * y ;

//   // beta.resize(lambda.n_elem, x.n_cols);
//   // df.resize(lambda.n_elem, 1);

//   // for (int i; i<lambda.n_elem; i++) {
//   //   // computing the structured ridge estimate
//   //   beta.row(i) = trans(cinvV * diagmat(eta/(square(eta) + lambda(i))) * Uty) / normx;
//   //   // computing the estimated degrees of freedom
//   //   df(i) = sum(eta/(square(eta) + lambda(i)));
//   // }

//   // // estimating the intercept term
//   // mu = ybar - beta * xbar;
// }

// void GROUP_LASSO::standardize(SEXP INTERCEPT, SEXP NORMALIZE) {

//   bool intercept  = as<bool>   (INTERCEPT) ; // boolean for intercept mode
//   bool normalize  = as<bool>   (NORMALIZE) ; // boolean for standardizing the predictor
//   uword n = x.n_rows;
//   uword p = x.n_cols;

//   if (intercept == 1) {
//     xbar = trans(rowvec(mean(x, 0)));
//     ybar = mean(y) ;
//     for (int i=0; i<p; i++) {
//       x.col(i) -= xbar(i);
//     }
//   } else {
//     xbar = zeros(p) ;
//     ybar = 0;
//   }

//   if (normalize == 1) {
//     normx = sqrt(trans(sum(square(x),0)));
//     for (int i=0; i<p; i++) {
//       x.col(i) /= normx(i);
//     }
//   } else {
//     normx = ones(p);
//   }
//   normy = sqrt(sum(square(y))) ;

//   if (intercept == 1) {
//     xty = trans(trans(y-ybar) * x) ;
//   } else {
//     xty = trans(y.t()*x) ;
//   }
// }

// void GROUP_LASSO::get_lambda(SEXP LAMBDA, SEXP NLAMBDA, SEXP LAMBDAMIN, SEXP LAMBDAMAX) {

//   if (LAMBDA != R_NilValue) {
//     lambda  = as<vec>(LAMBDA)  ;
//   } else {
//     lambda = exp10(linspace(log10(as<double>(LAMBDAMIN)), log10(as<double>(LAMBDAMAX)), as<uword>(NLAMBDA))) ;
//   }

// }
