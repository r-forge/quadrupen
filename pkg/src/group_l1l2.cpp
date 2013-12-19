// /*
//  * Author: Julien CHIQUET
//  *         Statistique et Génome
//  */

// #include "group_l1l2.h"

// using namespace Rcpp ;
// using namespace arma ;

// SEXP group_l1l2(SEXP N, SEXP X, SEXP XTY, SEXP PK, SEXP S, SEXP LAMBDA, SEXP WEIGHTS, SEXP XBAR, SEXP YBAR, SEXP NORMX, SEXP EPS, SEXP ZERO, SEXP MAXITER, SEXP MAXFEAT, SEXP FUN, SEXP VERBOSE, SEXP SPARSE) {

//   // Reading input variables
//   arma::uword n        = Rcpp::as<int>       (N)         ; // sample size
//   Rcpp::List  SX       = Rcpp::List(X)                   ; // sparsely encoded SCALED design matrix
//   arma::vec   xty      = Rcpp::as<arma::vec> (XTY)       ; // reponses to predictors vector
//   arma::uvec  pk       = Rcpp::as<arma::uvec> (PK)       ; // # of element in each group
//   arma::vec   lambda   = Rcpp::as<arma::vec> (LAMBDA)    ; // penalty levels
//   arma::vec   weights  = Rcpp::as<arma::vec> (WEIGHTS)   ; // penality weights
//   arma::vec   xbar     = Rcpp::as<arma::vec> (XBAR)      ; // mean of the predictors
//   double      ybar     = Rcpp::as<double>    (YBAR)      ; // mean of the response
//   arma::vec   normx    = Rcpp::as<arma::vec> (NORMX)     ; // norm of the predictors
//   double      eps      = Rcpp::as<double>    (EPS)       ; // precision required
//   double      zero     = Rcpp::as<double>    (ZERO)      ; // practical zero
//   arma::uword fun      = Rcpp::as<int>       (FUN)       ; // solver (0=quadra, 1=pathwise, 2=fista)
//   bool        verbose  = Rcpp::as<bool>      (VERBOSE)   ; // boolean for verbose mode
//   bool        sparse   = Rcpp::as<bool>      (SPARSE )   ; // boolean for sparse mode
//   arma::uword max_iter = Rcpp::as<int>       (MAXITER)   ; // max # of iterates of the active set
//   arma::uword max_feat = Rcpp::as<int>       (MAXFEAT)   ; // max # of variables activated

//   // Managing the sparse encoding of the structuring matrix
//   Rcpp::List SS                                          ; // sparsely encoded structuring matrix
//   arma::uvec Si                                          ; // row indices of nonzeros
//   arma::uvec Sp                                          ; // col indices of nonzeros
//   arma::vec  Sx                                          ; // values of nonzeros
//   arma::vec  Sxnzero                                     ; // temporary variables used for
//   arma::uvec Sjnzero                                     ; // updating the Gram matrix
//   double gamma = 0                                       ; // the smooth (ridge) penality
//   if (S != R_NilValue) { // Check if their is any structure to read
//     SS    = Rcpp::List(S)                                ;
//     Si    = Rcpp::as<arma::uvec>(SS[0])                  ;
//     Sp    = Rcpp::as<arma::uvec>(SS[1])                  ;
//     Sx    = Rcpp::as<arma::vec> (SS[2])                  ;
//     gamma = Sx[0]                                        ;
//   }

//   // Managing the design matrix is both cases of sparse or dense coding
//   arma::mat  x                                           ; // dense encoding of the design matrix
//   arma::uvec Xi                                          ; // row indices of nonzeros
//   arma::uvec Xp                                          ; // indices of first nonzero of each column
//   arma::uvec Xnp                                         ; // # of nonzero in each column
//   arma::vec  Xx                                          ; // values of nonzeros
//   arma::uvec j_nz                                        ; // consider just the column of X which are non zero
//   arma:: vec col_Vx                                      ;
//   arma::uvec col_Xi                                      ;
//   arma:: vec col_Xx                                      ;
//   if (sparse == 1) { // Check how x is encoded for reading
//     Xi        = Rcpp::as<arma::uvec>(SX[0]);
//     Xp        = Rcpp::as<arma::uvec>(SX[1]);
//     Xnp       = Rcpp::as<arma::uvec>(SX[2]);
//     Xx        = Rcpp::as<arma::vec> (SX[3]);
//     j_nz = find(Xnp > 0);
//   } else {
//     x = Rcpp::as<arma::mat>(SX[0]) - arma::ones(n,1) * trans(xbar) ;
//   }

//   // Initializing "first level" variables (outside of the lambda loop)
//   arma::uword p        = xty.n_elem                      ; // problem size
//   arma::uword K        = pk.n_elem                       ; // # of group
//   arma::uword n_lambda = lambda.n_elem                   ; // # of penalty
//   arma::vec  mu        = arma::zeros<arma::vec>(n_lambda)  ; // vector of intercept
//   arma::uvec GA                                            ; // set of currently activated groups
//   arma::uvec A                                             ; // set of currently activated variables
//   arma::vec  betaA                                         ; // vector of currently activated parameters
//   arma::vec  new_col                                       ; // column currently added to xtxA
//   arma::mat  xtxA                                          ; // t(x) * x_A  covariance matrix
//   arma::mat  xAtxA                                         ; // t(x_A) * x_A covariance matrix of the activated variable
//   arma::vec  xtxw                                          ; // t(x_A) * x_A * beta(A)
//   arma::vec  grd       = -xty                              ; // smooth part of the gradient
//   arma::vec  max_grd   = arma::zeros<arma::vec>(n_lambda)  ; // a vector with the successively reach duality gap
//   arma::vec  converge  = arma::zeros<arma::vec>(n_lambda)  ; // a vector indicating if convergence occured (0/1/2)
//   arma::uvec it_active = arma::zeros<arma::uvec>(n_lambda) ; // # of loop in the active set for each lambda
//   arma::uvec it_optim                                      ; // # of loop in the optimization process for each loop of the active se
//   double L0            = 1.0 + gamma                       ; // Lipschitz constant for proximal methods
//   arma::vec  timing      (n_lambda)                        ; // succesive timing in
//   arma::wall_clock timer                                   ; // clock

//   // Initializing "second level" variables (within the active set - for a fixed value of lamdba)
//   int        iter     = 0   ; // current iterate
//   arma::uword grp_in        ; // currently added group
//   arma::uvec  gk_start (K)  ; // starting index for each group
//   arma::uvec  gk_end   (K)  ; // starting index for each group
//   gk_start[0] = 0 ;
//   gk_end[0]   = pk[0];
//   for (int k = 1; k < K; k++) {
//     gk_start[k] = gk_start[k-1] + pk[k-1];
//     gk_end[k]   = gk_start[k] + pk[k];
//   }
//   arma::uvec  group   (p)  ; // indicator of group belonging
//   for (int k = 0; k < K; k++) {
//     group.subvec(gk_start[k], gk_end[k]-1) = k*arma::ones<arma::uvec>(pk[k]);
//   }

//   int        nbr_grp_in = 0 ; // # of currently added groups
//   arma::uword var_in        ; // currently added variables
//   int        nbr_in   = 0   ; // # of currently added variables
//   int        nbr_opt  = 0   ; // # of current calls to the optimization routine

//   arma::uvec are_in  = arma::zeros<arma::uvec>(K) ; // a vector to check if a group is already in the active set
//   Rcpp::List out_optim      ; // the list of output of the optimization function
//   arma::uvec null           ; // stores the variables which go to zero during optimization
//   arma::vec grd_norm (p)    ; // current value of the grd_norm for each variable
//   arma::vec pen      (p)    ; // current vector of penalties
//   arma::mat nonzeros        ; // contains non-zero value of beta
//   arma::mat iA              ; // contains row indices of the non-zero values
//   arma::mat jA              ; // contains column indices of the non-zero values

//   // _____________________________________________________________
//   //
//   // START THE LOOP OVER LAMBDA
//   timer.tic();
//   for (int m=0; m<n_lambda; m++) {
//     if (verbose == 1) {Rprintf("\n lambda = %f",lambda(m));}
//     // _____________________________________________________________
//     //
//     // START THE ACTIVE SET ALGORITHM
//     // _____________________________________________________________
//     //

//     // current vector of weighted penalities
//     pen = lambda[m]*weights;
//     // dual norm of gradient for each group
//     grd_norm = grp_norm(grd, pk, 2, 1) - pen;
//      // variable associated with the highest optimality violation
//     max_grd[m]  = grd_norm.max(var_in) ;
//     // group associated with the highest optimality violation
//     grp_in = group(var_in) ;

//     while (max_grd[m] > eps & it_active[m] < max_iter) {
//       // _____________________________________________________________
//       //
//       // (1) GROUP ACTIVATION IF APPLICABLE
//       // ____________________________________________________________

//       // Check if the group is already in the active set
//       if (are_in[grp_in] == 0) {
// 	// If this is a newly added group, then
// 	GA.resize(nbr_grp_in+1) ; // update the set of active groups
// 	GA[nbr_grp_in] = grp_in ;
//  	A.resize(nbr_in+pk[grp_in]) ; // resize the active set
// 	betaA.resize(nbr_in+pk[grp_in]) ; // resize the vector of active paramters

// 	// Update the gram matrix and Cholesky decomposition elementwise
//  	for (int k=gk_start[grp_in]; k<gk_end[grp_in]; k++) {
// 	  betaA[nbr_in] = 0 ;
// 	  A[nbr_in] = k;

// 	  // Update the xtxA and xAtxA matrices
// 	  if (nbr_in > 0) {
// 	    xAtxA = join_cols(xAtxA, xtxA.row(k)) ;
// 	  }

// 	  if (sparse == 1) {
// 	    // if any nonzero in the X[, var] column, do the sparse product
// 	    new_col = arma::zeros<arma::vec>(p) ;
// 	    if (Xnp[k] > 0) {
// 	      col_Vx = arma::zeros<arma::vec>(n) ;
// 	      col_Vx.elem(Xi.subvec(Xp[k],Xp[k+1]-1)) = Xx.subvec(Xp[k],Xp[k+1]-1);
// 	      // loop along each column of X
// 	      for (int j=0; j<j_nz.n_elem; j++) {
// 		col_Xx = Xx.subvec(Xp[j_nz[j]],Xp[j_nz[j]+1]-1) ;
// 		col_Xi = Xi.subvec(Xp[j_nz[j]],Xp[j_nz[j]+1]-1) ;
// 		new_col[j_nz[j]] = dot(col_Xx, col_Vx.elem(col_Xi));
// 	      }
// 	    }
// 	    new_col = new_col - n * xbar * arma::as_scalar(xbar(k));
// 	  } else {
// 	    new_col = trans(x) * x.col(k);
// 	  }

// 	  if (gamma > 0) {
// 	    // Sparse product with the structurating matrix
// 	    Sxnzero = Sx.subvec(Sp[k],Sp[k+1]-1) ;
// 	    Sjnzero = Si.subvec(Sp[k],Sp[k+1]-1) ;
// 	    new_col.elem(Sjnzero) +=  new_col.elem(Sjnzero) % Sxnzero ;
// 	  }

// 	  xtxA  = join_rows(xtxA, new_col) ;
// 	  xAtxA = join_rows(xAtxA, trans(xtxA.row(k))) ;
// 	  if (fun == 1) {
// 	    xtxw.resize(nbr_in+1) ;
// 	    xtxw(nbr_in) = dot(xAtxA.col(nbr_in),betaA.subvec(0,nbr_in));
// 	  }
// 	  nbr_in++;
// 	}
//  	if (verbose == 1) {Rprintf("\nnewly added group %i\n",grp_in);}
// 	are_in[grp_in] = 1;
//  	nbr_grp_in++;
//       }

//       // _____________________________________________________________
//       //
//       // (2) OPTIMIZATION OVER THE CURRENTLY ACTIVATED VARIABLES
//       // _____________________________________________________________
//       //
// //       if (fun == 1) {
// // 	  out_optim = pathwise_grpl1l2(betaA, pk.elem(GA), xAtxA, xty.elem(A), xtxw, pen.elem(A), gamma, eps);
// // 	  xtxw = as<arma::vec>(out_optim[2]);
// //       }
//       if (fun == 2) {
// 	out_optim = fista_grp(betaA, pk.elem(GA), xAtxA, xty.elem(A), pen.elem(A), 2, L0, eps, zero);
// 	L0 = as<double>(out_optim[2]);
//       }
//       // save the number of iterates performed along the optimization process
//       it_optim.reshape(nbr_opt + 1,1) ;
//       it_optim[nbr_opt] = Rcpp::as<int>(out_optim[0]) + 1 ;
//       nbr_opt++;

//       // _____________________________________________________________
//       //
//       // (3) VARIABLE DELETION IF APPLICABLE
//       // _____________________________________________________________
//       //
//       betaA  = as<arma::vec>(out_optim[1]);
//       grd = -xty + xtxA * betaA;
//       null = arma::sort(arma::find(abs(grp_norm(betaA, pk.elem(GA), 2, 0) + grp_norm(grd.elem(A) - weights.elem(A) * lambda(m), pk.elem(GA), 2, 0))/p < zero),1) ;
//       if (!null.is_empty()) {
//  	for (int k=0; k<null.n_elem; k++) {
// 	  arma::uvec remove  = sort(find(group.elem(A) == null[k]),1) ;
// 	  for (int j = 0; j<remove.n_elem; j++) {
// 	    if (verbose == 1) {Rprintf("removing group %i\n",remove[j]);}
// 	    A.shed_row(remove[j])     ;
// 	    betaA.shed_row(remove[j]) ;
// 	    xtxA.shed_col(remove[j])  ;
// 	    xAtxA.shed_col(remove[j]) ;
// 	    xAtxA.shed_row(remove[j]) ;
// 	    if (fun == 1) {
// 	      xtxw.shed_row(remove[j]);
// 	    }
// 	    nbr_in--;
// 	  }
// 	  GA.shed_row(null[k]);
// 	  nbr_grp_in--;
// 	}
//       }

//       // _____________________________________________________________
//       //
//       // (4) OPTIMALITY TESTING
//       // _____________________________________________________________

//       // dual norm of gradient for each group
//       grd_norm = grp_norm(grd, pk, 2, 1) - pen;
//       // variable associated with the highest optimality violation
//       max_grd[m]  = grd_norm.max(var_in) ;
//       // group associated with the highest optimality violation
//       grp_in = group(var_in) ;

//       // Moving to the next iterate
//       it_active[m]++;
//       R_CheckUserInterrupt();
//     }

//     // Preparing next value of the penalty
//     nonzeros = join_cols(nonzeros, betaA/normx.elem(A))    ;
//     iA = join_cols(iA, m*arma::ones(betaA.n_elem,1) )      ;
//     jA = join_cols(jA, arma::conv_to<arma::mat>::from(A) ) ;

//     mu[m] = ybar - dot(xbar.elem(A), betaA);

//     if (it_active[m] >= max_iter) {
//       converge[m] = 1;
//     }

//     timing(m) = timer.toc();
//     if (nbr_grp_in > max_feat) {
//       converge[m] = 2 ;
//       lambda      =    lambda.subvec(0,m) ;
//       converge    =  converge.subvec(0,m) ;
//       mu          =        mu.subvec(0,m) ;
//       max_grd     =   max_grd.subvec(0,m) ;
//       it_active   = it_active.subvec(0,m) ;
//       timing      =    timing.subvec(0,m) ;
//       break;
//     }

//   }
//   // END OF THE LOOP OVER LAMBDA

//   return Rcpp::List::create(Rcpp::Named("mu")        = mu       ,
// 			    Rcpp::Named("nzeros")    = nonzeros ,
// 			    Rcpp::Named("iA")        = iA       ,
// 			    Rcpp::Named("jA")        = jA       ,
//  			    Rcpp::Named("lambda")    = lambda   ,
//  			    Rcpp::Named("nbr.in")    = nbr_in   ,
//  			    Rcpp::Named("it.active") = it_active,
//  			    Rcpp::Named("it.optim")  = it_optim ,
//  			    Rcpp::Named("max.grd")   = max_grd ,
//  			    Rcpp::Named("timing")    = timing   ,
//  			    Rcpp::Named("converge")  = converge );

// }
