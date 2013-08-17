#include "utils.h"

using namespace Rcpp;
using namespace arma;

void standardize_dense(mat &x, vec &y, bool &intercept, bool &normalize, vec &weights,
		       vec &xty, rowvec &normx, double &normy, rowvec &xbar, double &ybar) {

  // uword n = x.n_rows;
  // uword p = x.n_cols;

  // mat meanx ;
  // mat xnorm ;

  // if (intercept == 1) {
  //   meanx = mean(x, 0);
  //   ybar = mean(y) ;
  // } else {
  //   meanx = zeros(p) ;
  //   ybar = 0;
  // }

  // if (normalize == 1) {
  //   xnorm = sqrt(sum(square(x),0) - n * square(meanx));
  //   x = x * diagmat(pow(xnorm,-1)) ;
  //   meanx = meanx/xnorm ;
  // } else {
  //   xnorm = ones(p, 1);
  // }
  // normy = sqrt(sum(square(y))) ;

  // if (any(weights != 1)) {
  //   x = x * diagmat(pow(weights,-1));
  //   meanx = meanx/weights;
  // }

  // if (intercept == 1) {
  //   xty = trans(x-ones(n,1) * meanx) * (y-ybar) ;
  // } else {
  //   xty = trans(x) * y;
  // }

  // xbar = as<vec>(meanx);
  // normx = trans(xnorm);

}

void standardize_sparse(sp_mat &x, vec &y, bool &intercept, bool &normalize, vec &weights,
			vec &xty, rowvec &normx, double &normy, rowvec &xbar, double &ybar) {

  // uword n = x.n_rows;
  // uword p = x.n_cols;

  // rowvec meanx ;
  // rowvec xnorm ;

  // if (intercept == 1) {
  //   meanx = mean(x, 0);
  //   ybar = mean(y) ;
  // } else {
  //   meanx = zeros(p) ;
  //   ybar = 0;
  // }

  // if (normalize == 1) {
  //   xnorm = sqrt(sum(square(x),0) - n * square(meanx));
  //   x = x * diagmat(pow(xnorm,-1)) ;
  //   meanx = meanx/xnorm ;
  // } else {
  //   xnorm = ones(p, 1);
  // }
  // normy = sqrt(sum(square(y))) ;

  // if (any(weights != 1)) {
  //   x = x * diagmat(pow(weights,-1));
  //   meanx = meanx/weights;
  // }

  // if (intercept == 1) {
  //   xty = trans(x-ones(n,1) * meanx) * (y-ybar) ;
  // } else {
  //   xty = trans(x) * y;
  // }

  // xbar = trans(meanx);
  // normx = trans(xnorm);
}

vec cg(mat A, vec b, vec x, double tol) {

  vec r = b - A * x;
  vec p = r ;
  double rs_old = sum(square(r)) ;

  double rs_new = rs_old ;
  int i = 0;
  double alpha ;
  mat Ap ;

  while (sqrt(rs_new) > tol & i < 1e3) {
    Ap = A * p;
    alpha = rs_old/dot(p,Ap) ;
    x += alpha * p ;
    r -= alpha * Ap ;
    // Polak-Ribière update
    rs_new = dot(r,-alpha * Ap);
    p = r + rs_new/rs_old*p;
    rs_old = rs_new;
    i++;
  }

  // Rprintf("\n nb of iterate %d",i);
  return(x);
}

// Can't find a reasonable Preconditioner that does not
// require a computational burden equivalent to a Cholesky decomposition
vec pcg(mat A, mat P, vec b, vec x, double tol) {

  vec r = b - A * x;
  vec z = P * r;
  vec p = z ;
  //double rs_old = sum(square(r)) ;
  double rs_old = dot(r,z) ;

  double rs_new = rs_old ;
  int i = 0;
  double alpha ;
  mat Ap ;

  while (sqrt(rs_new) > tol & i < 1e3) {
    Ap = A * p;
    alpha = rs_old/dot(p,Ap) ;
    x += alpha * p ;
    r -= alpha * Ap ;
    // Polak-Ribière update
    z = P * r ;
    rs_new = dot(z,-alpha * Ap);
    p = z + rs_new/rs_old*p;
    rs_old = rs_new;
    i++;
  }

  // Rprintf("\n nb of iterate %d",i);
  return(x);
}
