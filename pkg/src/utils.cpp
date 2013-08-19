#include "utils.h"

using namespace Rcpp;
using namespace arma;

void standardize(mat &x, vec &y, bool &intercept, bool &normalize, vec &weights,
		 vec &xty, vec &normx, double &normy, vec &xbar, double &ybar) {

  uword n = x.n_rows;
  uword p = x.n_cols;

  if (intercept == 1) {
    xbar = trans(rowvec(mean(x, 0)));
    ybar = mean(y) ;
  } else {
    xbar = zeros(p) ;
    ybar = 0;
  }

  if (normalize == 1) {
    normx = sqrt(trans(sum(square(x),0)) - n * square(xbar));
    for (int i=0; i<p; i++) {
      x.col(i) /= normx(i);
    }
    xbar /= normx ;
  } else {
    normx = ones(p);
  }
  normy = sqrt(sum(square(y))) ;

  if (any(weights != 1)) {
    for (int i=0; i<n; i++) {
       x.row(i) /= weights ;
    }
    xbar /= weights;
  }

  if (intercept == 1) {
    xty = trans(trans(y-ybar) * x) ;
    for (int i=0;i<p;i++) {
       xty(i) -=  sum(y-ybar) * xbar(i);
    }
  } else {
    xty = trans(y.t()*x) ;
  }
}

void standardize(sp_mat &x, vec &y, bool &intercept, bool &normalize, vec &weights,
		 vec &xty, vec &normx, double &normy, vec &xbar, double &ybar) {

  uword n = x.n_rows;
  uword p = x.n_cols;

  if (intercept == 1) {
    xbar = trans(rowvec(mean(x, 0)));
    ybar = mean(y) ;
  } else {
    xbar = zeros(p) ;
    ybar = 0;
  }

  if (normalize == 1) {
    normx = sqrt(trans(sum(square(x),0)) - n * square(xbar));
    for (int i=0; i<p; i++) {
      x.col(i) /= normx(i);
    }
    xbar /= normx ;
  } else {
    normx = ones(p);
  }
  normy = sqrt(sum(square(y))) ;

  if (any(weights != 1)) {
    for (int i=0; i<n; i++) {
       x.row(i) /= weights ;
    }
    xbar /= weights;
  }

  if (intercept == 1) {
    xty = trans(trans(y-ybar) * x) ;
    for (int i=0;i<p;i++) {
       xty(i) -=  sum(y-ybar) * xbar(i);
    }
  } else {
    xty = trans(y.t()*x) ;
  }
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
