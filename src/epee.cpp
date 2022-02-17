#include "epee.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

namespace epee {

void Epee::run (mat& eta, Covariance& K)
{
  if (mx == 0)
  {
    abort ();
    return;
  }
  start (eta, K);
  Covariance Kd (K); //must have a copy ...
  for(iter=0; iter<maxiter; ++iter)
  {
    update (eta, K, Kd);
    if (std::isnan(Z(iter)))
    {
      Z.resize(iter+1);
      break;
    }
    if (iter > 0 && fabs(Z(iter) - Z(iter-1)) < tol) 
    {
      converged = true;
      Z.resize(iter+1);
      break;
    }
  }
  gradient (K, Kd);
}

void Epee::abort (void)
{
  dL = arma::zeros<mat>(mx+mz, mx+mz);
  dS = arma::zeros<mat>(mx+mz, mx+mz);
  dm = arma::zeros<mat>(mx, obs.n_cols);
  Z(0) = 0.;
  Z.resize(1);
  iter = 0;
  converged = true;
}

void Epee::start (mat& eta, Covariance& K)
{
  eta.elem(arma::find_nonfinite(obs)).fill(arma::datum::nan);//use .replace
  K.traverse_tree (eta, arma::zeros<mat>(arma::size(obs)));
  kappa.fill(arma::datum::nan);
  kappa.elem(m) = K.Y.elem(m);
  
  mu.fill(arma::datum::nan);
  mu.elem(m) = eta.elem(m);
  sigma.elem(m) = 1./K.diag.elem(m);
  if (fix_start)
  {
    mu.elem(m).zeros();
    sigma.elem(m).ones();
  }

  Z.set_size(maxiter);
  Z.fill(arma::datum::nan);
  tau_cavity.fill(arma::datum::nan);
  tau_tilde.fill(arma::datum::nan);
  tau_hat.fill(arma::datum::nan);
  mu_cavity.fill(arma::datum::nan);
  mu_hat.fill(arma::datum::nan);
  nu_tilde.fill(arma::datum::nan);
  Z_hat.fill(arma::datum::nan);
  beta.fill(arma::datum::nan);
  alpha.fill(arma::datum::nan);
  tau_tilde.elem(m).zeros();
  nu_tilde.elem(m).zeros();
  converged = false;
}

void Epee::update (const mat& eta, Covariance& K, Covariance& Kd)
{
  // TODO: arguably don't need all the elem shit as NaN should just stay that way
  /* normal EP updates */
  tau_cavity.elem(m)  = 1./sigma.elem(m) - tau_tilde.elem(m); // remember to use elem
  mu_cavity.elem(m)   = 1./tau_cavity.elem(m) % (mu.elem(m)/sigma.elem(m) - nu_tilde.elem(m));
  alpha.elem(m) = mu_cavity.elem(m)/arma::sqrt(2./tau_cavity.elem(m));
  Z_hat.elem(m) = 0.5*(1. - erf(-alpha.elem(m))) - (1. - obs.elem(m));
  beta.elem(m) = 1./Z_hat.elem(m) % arma::sqrt(1./(2.*arma::datum::pi*tau_cavity.elem(m))) % arma::exp(-alpha.elem(m)%alpha.elem(m));
  mu_hat.elem(m)      = mu_cavity.elem(m) + beta.elem(m); 
  tau_hat.elem(m)     = 1./(mu_cavity.elem(m)%mu_cavity.elem(m) + 1./tau_cavity.elem(m) + 
                            mu_cavity.elem(m)%beta.elem(m) - mu_hat.elem(m)%mu_hat.elem(m)); 
  tau_tilde.elem(m)   = tau_hat.elem(m) - tau_cavity.elem(m);
  nu_tilde.elem(m)    = mu_hat.elem(m)%tau_hat.elem(m) - mu_cavity.elem(m)%tau_cavity.elem(m);

  /* update sigma and mu */
  beta.elem(m)        = 1./tau_tilde.elem(m);
  alpha.elem(m)       = beta.elem(m) % (kappa.elem(m) + nu_tilde.elem(m));
  Kd.traverse_tree (alpha, beta);
  sigma.elem(m)       = beta.elem(m) % (1. - beta.elem(m) % Kd.diag.elem(m));
  mu.elem(m)          = alpha.elem(m) - beta.elem(m) % Kd.Y.elem(m);

  /* normalizing constant */
  Z(iter) = 0.5 * (
             Kd.logdet -
             arma::accu(arma::log(tau_tilde.elem(m))) -
             K.Ayy +
             arma::accu(alpha.elem(m) % (kappa.elem(m) + nu_tilde.elem(m))) -
             Kd.Ayy +
             arma::accu(2. * arma::log(arma::abs(Z_hat.elem(m))) + 
                        arma::log(1. + tau_tilde.elem(m)/tau_cavity.elem(m)) +
                        (tau_tilde.elem(m) % mu_cavity.elem(m) % mu_cavity.elem(m) -
                         2. * mu_cavity.elem(m) % nu_tilde.elem(m) -
                         nu_tilde.elem(m) % nu_tilde.elem(m) / tau_cavity.elem(m))/
                        (1. + tau_tilde.elem(m)/tau_cavity.elem(m))));
}

void Epee::gradient (Covariance& K, Covariance& Kd)
{
  mat dY = mu;
  dY.elem(m) = 2. * mu.elem(m);
  dY(arma::find_nonfinite(dY)).zeros();//TODO use replace method
  Kd.traverse_reverse(arma::zeros<mat>(arma::size(mu)), 1., -1.);
  K.traverse_reverse(dY, 0., -1.);
  dL = 0.5 * (Kd.dL + K.dL);
  dS = 0.5 * (Kd.dS + K.dS);
  dm = 0.5 * (K.dy);
}

//inline mat Epee::erf (mat x)
//{
//  for (uword i=0; i<x.n_rows; ++i)
//    for (uword j=0; j<x.n_cols; ++j)
//      x(i,j) = 2.0 * R::pnorm(x(i,j) * sqrt(2.0), 0.0, 1.0, true, false) - 1.0;
//  return x;
//}

vec Epee::erf (vec x)
{
  for (uword i=0; i<x.n_elem; ++i)
      x(i) = 2.0 * R::pnorm(x(i) * sqrt(2.0), 0.0, 1.0, true, false) - 1.0;
  return x;
}

void Epee::options (const double tolerance, const uword maxitt, const bool fixstart)
{
  tol = tolerance;
  maxiter = maxitt;
  fix_start = fixstart;
}

} /* end namespace epee */
