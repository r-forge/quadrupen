/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */
#ifndef _quadrupen_GROUP_LASSO_H
#define _quadrupen_GROUP_LASSO_H

#define ARMA_NO_DEBUG
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS

#ifndef ARMA_HAVE_GETTIMEOFDAY
#define ARMA_HAVE_GETTIMEOFDAY
#endif

#include "quadrupen_headers.hpp"

using namespace Rcpp;
using namespace arma;

RcppExport SEXP group_lasso(SEXP BETA0    ,
			    SEXP X        , // matrix of features
			    SEXP Y        , // vector of response	    
			    SEXP TYPE     , // -1="l1/linf", 2="l1/l2" group-Lasso
			    SEXP PK       , // successive groups length
			    SEXP OMEGA    , // structuring matrix
			    SEXP LAMBDA1  , 
			    SEXP NLAMBDA1 ,
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
			    SEXP SPARSE   ) ;

class GROUP_LASSO {

 public:
  GROUP_LASSO(SEXP X        , // matrix of features
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
	      SEXP SPARSE   );
  
  void get_lambda1(SEXP LAMBDA1, SEXP NLAMBDA, SEXP MIN_RATIO);

  void LetsRoll();

  // vec grp_norm(vec x, uvec pk, int rep) {};

  // various function to acces private members
  const arma::mat & get_coef() const { return beta; }
  const arma::vec & get_lambda() const { return lambda1; }
  const arma::vec & get_df() const { return df; }
  const arma::vec & get_normx() const { return normx; }
  const arma::vec & get_mu() const { return mu; }

 private:
  arma::mat       x  ; // matrix of predictors (dense coding)
  arma::sp_mat sp_x  ; // matrix of predictors (sparse coding)
  arma::mat xt       ; // transpose once and keep it, since this
  arma::sp_mat sp_xt ; // costs much for sparse matrix)

  vec y           ; // vector of reponses
  vec pk          ; // vector of group numbers
  vec    weights  ; // observation weights (not use at the moment)
  vec    penscale ; // penalty weights
  double lambda2  ; // vector of tuning parameters for l2-penalty
  bool intercept  ; // boolean for intercept mode
  bool normalize  ; // boolean for standardizing the predictor
  bool   naive    ; // naive elastic-net or not
  double eps      ; // precision required
  uword  fun      ; // solver (0=quadra, 1=pathwise, 2=fista)
  int    verbose  ; // int for verbose mode (0/1/2)
  bool   sparse   ; // boolean for sparse mode
  uword  max_iter ; // max # of iterates of the active set
  uword  max_feat ; // max # of variables activated

  arma::mat omega    ; // matrix structuring the l2-penalty
  arma::vec xty      ; // reponses to predictors vector
  arma::vec xbar     ; // mean of the predictors
  arma::vec normx    ; // norm of the predictors
  double normy       ; // norm of the response
  double ybar        ; // mean of the response

  arma::vec lambda1  ; // vector of tuning parameters for l1-penalty
  arma::mat beta     ; // matrix of coefficients
  arma::vec df       ; // degrees of freedom along the path
  arma::vec mu       ; // vector of intercept term

};

#endif
