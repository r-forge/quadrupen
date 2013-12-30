/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */
#ifndef _quadrupen_DATA_REG_H
#define _quadrupen_DATA_REG_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// TODO: use template to handle dense or sparse encoding (mat/sp_mat in armadillo) 
class REGRESSION_DATA {

protected:
  // DATA VARIABLES FOR REGRESSION PURPOSE
  uword  n          ; // sample size
  uword  p          ; // # of features
  mat    x          ; // matrix of predictors
  vec    xbar       ; // mean of the predictors
  vec    normx      ; // norm of the predictors
  vec    y          ; // vector of response
  double ybar       ; // mean of the response
  double normy      ; // norm of the response
  bool   intercept  ; // should intercept be considered?
  bool   normalize  ; // should predictors be standardized?
  vec    obsweights ; // observation weights
  vec    preweights ; // predictor weights
  mat    xt         ; // transpose matrix of predictor once and keep it to save time
  vec    xty        ; // reponses to predictors vector

public:
  REGRESSION_DATA(SEXP X          , // matrix of features
		  SEXP Y          , // vector of response
		  SEXP INTERCEPT  , // should intercept be considered?
		  SEXP NORMALIZE  , // should predictors be standardized?
		  SEXP PREWEIGHTS , // prediction weights
		  SEXP OBSWEIGHTS); // observation weights
  
  // DATA NORMALIZATION
  // intercept treatment, predictor standardization, predictor weihgts and observation weights
  void standardize();

  // various function to acces private members
  const mat & get_x()   const { return x   ; }
  const vec & get_y()   const { return y   ; }
  const vec & get_xty() const { return xty ; }
  
};

#endif
