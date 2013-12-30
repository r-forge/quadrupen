#include "active_set.hpp"

using namespace Rcpp;
using namespace arma;

ACTIVE_SET::ACTIVE_SET(SEXP BETA0)
{
  vec beta0 = as<vec>(BETA0) ;
  active_var = find(beta0 != 0) ;
  beta_active = beta0.elem(active_var) ;
  nbrIn = active_var.n_elem;
  is_active_var.zeros(beta0.n_elem) ;
  for (int i=0; i<nbrIn;i++) {
    is_active_var(active_var(i)) = 1;
  }
}

void ACTIVE_SET::add_var(uword varIn) {
  active_var.resize(nbrIn+1) ; // update the active set
  active_var[nbrIn] = varIn  ;
  beta_active.resize(nbrIn+1); // update the vector of active parameters
  beta_active[nbrIn]  = 0.0  ;
  nbrIn++                    ; // update the number of active variable
  is_active_var[varIn] = 1   ;
}

void ACTIVE_SET::del_var(uword indVarOut) {
  is_active_var[active_var[indVarOut]] = 0 ; // update the active set
  active_var.shed_row(indVarOut)           ;
  beta_active.shed_row(indVarOut)          ; // update the vector of active parameters
  nbrIn--                                  ; // update the number of active variable
}

ACTIVE_SET_GROUPWISE::
ACTIVE_SET_GROUPWISE(SEXP BETA0, SEXP PK):
  ACTIVE_SET(BETA0),
  pk (as<uvec>(PK)) // a vector with the group numbers
{
  K = pk.n_elem          ; // number of groups
  gk_start.zeros(K)      ; // starting index for each group
  gk_end.zeros  (K)      ; // starting index for each group
  group.zeros(sum(pk))   ; // indicator of group belonging
  is_active_grp.zeros(K) ;
  nbrGrpIn = 0 ;

  gk_start[0] = 0 ;
  gk_end[0]   = pk[0];
  for (int k = 1; k < K; k++) {
    gk_start[k] = gk_start[k-1] + pk[k-1];
    gk_end[k]   = gk_start[k] + pk[k];
  }
  for (int k = 0; k < K; k++) {
    group.subvec(gk_start[k], gk_end[k]-1) = k*arma::ones<arma::uvec>(pk[k]);
  }
}

void ACTIVE_SET_GROUPWISE::add_var(uword grpIn) {

  active_grp.resize(nbrGrpIn+1) ; // update the active set
  active_grp[nbrGrpIn] = grpIn  ;
  nbrGrpIn++                    ; // update the number of active variable
  is_active_grp[grpIn] = 1      ;

  for (int k=gk_start[grpIn]; k<gk_end[grpIn]; k++) {
    ACTIVE_SET::add_var(k);
  }
}

void ACTIVE_SET_GROUPWISE::del_var(uword indGrpOut) {
  is_active_grp[active_grp[indGrpOut]] = 0; // update the active set
  active_grp.shed_row(indGrpOut)            ;
  nbrGrpIn--                                ; // update the number of active variable

  // Update set of active variables
  uvec remove = sort(find(group.elem(active_var) == active_grp[indGrpOut]),1) ;
  for (int k=0; k<remove.n_elem; k++) {
    ACTIVE_SET::del_var(remove[k]);
  }
}

// ACTIVE_ALGORITHM
// vec GROUP_LASSO::get_dist2opt(double lambda) {

//   vec normGrad = max(zeros(K), grp_norm(smooth_grad) - lambda) ;
//   // normGrad.elem(active_group) = grp_norm(smooth_grad + lambda * betaA);

//   return(normGrad) ;
// }

