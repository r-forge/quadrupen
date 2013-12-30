/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */

#include "group_lasso_main.hpp"

using namespace Rcpp;
using namespace arma;

SEXP group_lasso(SEXP BETA0    ,
		 SEXP X        , // matrix of features
		 SEXP Y        , // vector of response
		 SEXP TYPE     , // -1="l1/linf", 2="l1/l2" group-Lasso
		 SEXP PK       , // successive groups length
		 SEXP OMEGA    , // structuring matrix
		 SEXP LAMBDA1  ,
		 SEXP NLAMBDA1 ,
		 SEXP MINRATIO ,
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
  // INSTANTIATION OF THE REQUIRED OBJECTS

  // data
  REGRESSION_DATA data(X, Y, INTERCEPT, NORMALIZE, PENSCALE, WEIGHTS) ;
  // path
  PATH path(MAXFEAT) ;
  // group_Lasso object: group structure, active set (groupwise), norm and dual norm.  
  // L1_NORM lasso ;
  // LINF_NORM bounded ;
  L1L2_NORM grpLasso1(PK) ;
  // L1LINF_NORM grpLasso2(PK) ;
  // COOP_NORM coopLasso(PK) ;

  PENALTY lasso(PK) ;
  lasso.setPenalty("l1") ;
    
  // ==============================================================
  // DATA NORMALIZATION AND BASICS PRETREATMENTS
  data.standardize();

  lasso.elt_norm(data.get_xty()).print();

  // ==============================================================
  // GET VALUES OF THE PENALTY LEADING THE PATH
  path.grid_penLevels(LAMBDA1, NLAMBDA1, MINRATIO, grpLasso1.dual_norm(data.get_xty()));
  
  // ==============================================================
  // COMPUTE THE PATH OF SOLUTIONS
  // _____________________________________________________________
  //
  // START THE LOOP OVER THE LEADING PENALTY
  // timer.tic();
  uword npen = path.get_penLevels().n_elem;
  // uword maxIter = as<uword>(MAXITER) ;
  uword grpIn ;
  bool verbose = as<bool>(VERBOSE) ;
  // vec dualGap = zeros(npen) ;
  // vec itActive = zeros(npen) ;
  // double eps = as<double>(EPS) ;

  ACTIVE_SET_GROUPWISE active_set(BETA0, PK);

  // vec dualGap.zeros(get_penlevels().n_elem) ;
  for (int m=0; m<npen; m++) {
    if (verbose == 1) {Rprintf("\n penalty level: %f",path.get_penLevels()[m]);}
    // _____________________________________________________________
    //
    // START THE ACTIVE SET ALGORITHM
    // _____________________________________________________________
    //
    // OPTIMALITY TESTING AT INITIALIZATION

    // group associated with the highest optimality violation
    // dualGap[m] = grpLasso.get_dist2opt(path.get_penLevels()[m]).max(grpIn) ;

    // while (dualGap[m] > eps & itActive[m] < maxIter) {
      // _____________________________________________________________
      //
      // (1) GROUP ACTIVATION IF APPLICABLE
      // ____________________________________________________________

    //   if (is_active_group[grpIn] == 0) {
    // 	add_group(grpIn) ;
    // 	if (verbose == 2) {Rprintf("Group %i newly added.\n",grpIn);}
    //   } else {
    // 	if (verbose == 2) {Rprintf("Group %i already in.\n",grpIn);}
    //   }
    //   // _____________________________________________________________
    //   //
    //   // (2) OPTIMIZATION OVER THE CURRENTLY ACTIVATED VARIABLES
    //   // _____________________________________________________________
    //   //

    //   // _____________________________________________________________
    //   //
    //   // (3) GROUP DELETION IF APPLICABLE
    //   // _____________________________________________________________
    //   //


    //   // _____________________________________________________________
    //   //
    //   // (4) OPTIMALITY TESTING
    //   // _____________________________________________________________

    //   itActive[m]++;
    //   R_CheckUserInterrupt(); // let the possibility to interrupt from R
    // } // ending the active set algorihm

  } // ending the loop over the penalty

  // return List::create(Named("coefficients") = grp_lasso.get_coef(),
  // 		      Named("mu")           = grp_lasso.get_mu()    ,
  // 		      Named("normx")        = grp_lasso.get_normx() ,
  // 		      Named("lambda1")      = grp_lasso.get_penalty(),
  // 		      Named("df")           = grp_lasso.get_df()    );

}


//   // _____________________________________________________________
//   //
//   // NEW DECLARATION
//   uword grpIn     ; // currently added group
//   uword varIn     ; // currently added variable

//   // _____________________________________________________________
//   //
//   // SOME INITIALIZATION
//   itActive.zeros(penalty.n_elem) ;
//   smooth_grad = -xty             ;
//   timing.zeros(penalty.n_elem)   ; // successive timings
//   df.zeros(penalty.n_elem)       ; // degrees of freedom
//   dualGap.zeros(penalty.n_elem)  ; // distance to optimality
//   is_active_group.zeros(K)       ;

