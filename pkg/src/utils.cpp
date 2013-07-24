#include "utils.h"

using namespace Rcpp;
using namespace arma;

// from Rcpp Gallery -- Dirk Eddelbuettel — Dec 25, 2012
sp_mat convertSparse(S4 mat) {
  IntegerVector dims = mat.slot("Dim");
  IntegerVector i = mat.slot("i");
  IntegerVector p = mat.slot("p");
  NumericVector x = mat.slot("x");

  int nrow = dims[0], ncol = dims[1];
  arma::sp_mat res(nrow, ncol);

  // create space for values, and copy
  arma::access::rw(res.values) = arma::memory::acquire_chunked<double>(x.size() + 1);
  arma::arrayops::copy(arma::access::rwp(res.values), x.begin(), x.size() + 1);

  // create space for row_indices, and copy -- so far in a lame loop
  arma::access::rw(res.row_indices) = arma::memory::acquire_chunked<arma::uword>(x.size() + 1);
  for (int j=0; j<i.size(); j++)
    arma::access::rwp(res.row_indices)[j] = i[j];

  // create space for col_ptrs, and copy -- so far in a lame loop
  arma::access::rw(res.col_ptrs) = arma::memory::acquire<arma::uword>(p.size() + 2);
  for (int j=0; j<p.size(); j++)
    arma::access::rwp(res.col_ptrs)[j] = p[j];

  // important: set the sentinel as well
  arma::access::rwp(res.col_ptrs)[p.size()+1] = std::numeric_limits<arma::uword>::max();

  // set the number of non-zero elements
  arma::access::rw(res.n_nonzero) = x.size();

  return(res);
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