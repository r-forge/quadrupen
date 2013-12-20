/***
 * @file interface.hpp
 *
 */
#ifndef __QUADRUPEN_HEADERS_HPP
#define __QUADRUPEN_HEADERS_HPP

// Include Armadillo / Rcpp / R to C/C++ basics
#include <sys/time.h>
#include <RcppArmadillo.h>

#define ZERO 2e-16 // practical zero

// Include utils and optimization routines
#include "utils/utils.h"
#include "optimization/quadratic.h"
#include "optimization/first_order.h"

#endif

