#ifndef EPEE_H
#define EPEE_H

#include "covariance.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

namespace epee {

class Epee
{
  private:
    const mat       obs;
    const uword     mx,
                    mz;
          mat       kappa,
                    mu,
                    sigma,
                    alpha,
                    beta,
                    tau_cavity,
                    mu_cavity,
                    tau_hat,
                    mu_hat,
                    Z_hat,
                    tau_tilde,
                    nu_tilde;
    const uvec      m;

  public:
          double    tol         = 1e-4;
          bool      converged   = false,
                    do_gradient = true,
                    fix_start   = false;
          uword     maxiter     = 100,
                    iter;
          vec       Z;
          mat       dL,
                    dS,
                    dm;

  Epee (const Covariance& K) :
    obs        (K.storage.yx),
    mx         (K.storage.tx.n_elem),
    mz         (K.storage.tz.n_elem),
    Z          (maxiter, arma::fill::zeros),
    kappa      (arma::size(obs)), 
    mu         (arma::size(obs)),
    sigma      (arma::size(obs)), 
    alpha      (arma::size(obs)),
    beta       (arma::size(obs)),
    tau_cavity (arma::size(obs)),
    mu_cavity  (arma::size(obs)),
    tau_hat    (arma::size(obs)),
    mu_hat     (arma::size(obs)),
    Z_hat      (arma::size(obs)),
    tau_tilde  (arma::size(obs)),
    nu_tilde   (arma::size(obs)),
    m          (arma::find_finite(obs)),
    dL         (0,0),
    dS         (0,0),
    dm         (arma::size(obs), arma::fill::zeros)
  {}

  void run      (mat&, Covariance&);
  void abort    (void);
  void start    (mat&, Covariance&);
  void update   (const mat&, Covariance&, Covariance&);
  void gradient (Covariance&, Covariance&);
  vec  erf      (vec);
  void options  (const double, const uword, const bool);
  /* overloads:
   * b/c the Covariance class isn't as fast as LAPACK for smallish matrices (<=100)
   * would be good to have a variant that runs with a typical matrix
   */
  //        void run      (mat&, CovarianceMatrix&);
  // inline void start    (mat&, CovarianceMatrix&);
  // inline void update   (mat&, CovarianceMatrix&);
  // inline void gradient (mat&, CovarianceMatrix&);
};

} /* end namespace epee */
#endif
