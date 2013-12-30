/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */
#ifndef _quadrupen_ACTIVE_SET_H
#define _quadrupen_ACTIVE_SET_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

class ACTIVE_SET {

protected:
  // VARIABLES FOR HANDLING THE ACTIVE SET
  uvec active_var    ; // set of currently activated variables
  uvec is_active_var ; // a vector to check if a variable is already in the active set
  uword nbrIn        ; // number of active variables
  vec  beta_active   ; // vector of currently activated features

public:
  ACTIVE_SET (SEXP BETA0);

  // ACTIVE SET HANDLING
  void add_var(uword varIn)     ; // add a variable in the active set
  void del_var(uword indVarOut) ; // remove the variable activated in postition ind_var_out

  const uvec & get_active_var() const { return active_var ; }
  const vec  & get_beta_active() const { return beta_active ; }

};

class ACTIVE_SET_GROUPWISE: public ACTIVE_SET {

protected:
  // VARIABLES FOR HANDLING THE ACTIVE SET
  uvec active_grp    ; // set of currently activated groups
  uvec is_active_grp ; // a vector to check if a variable is already in the active set
  uword nbrGrpIn     ; // number of active variables
  uvec  pk           ; // vector of group numbers
  uword K            ; // number of groups
  uvec  gk_start     ; // starting index for each group
  uvec  gk_end       ; // starting index for each group
  uvec  group        ; // indicator of group belonging for each variable

public:
  ACTIVE_SET_GROUPWISE (SEXP BETA0, SEXP PK);

  // ACTIVE SET HANDLING
  void add_var(uword grpIn)   ; // add a group of variables in the active set
  void del_var(uword indGrpOut) ; // remove the group of variable activated in postition ind_var_class

  const uvec   & get_active_grp() const { return active_grp ; }

};

// class ACTIVE_ALGORITHM_GROUPWISE {

// protected:
//   uword maxIter   ;
//   uword itActive  ;
//   double thres    ;
//   double dualGap  ;
//   vec smooth_grad ;

//   ACTIVE_SET_GROUPWISE active_set ;

// public:
//   ACTIVE_ALGORITHM(SEXP BETA0, SEXP PK, SEXP MAXITER);

// }

#endif

