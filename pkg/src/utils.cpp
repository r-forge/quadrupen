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
