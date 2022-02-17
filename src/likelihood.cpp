#include "likelihood.h"
//#include <gperftools/profiler.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

namespace epee {

double Regression::loglikelihood (const vec& L, const vec& S, const std::vector<vec>& b)
{ /* calculate loglikelihood */
  mat Lambda  = constraints (L, constraints.Lambda); // Choleski factors
  mat Sigma   = constraints (S, constraints.Sigma);
  covariance.parameters (Lambda * Lambda.t(), Sigma * Sigma.t(), design(b));
  double ll = 0;

  // I think this is OK b/c shouldn't update and all quantities shuld be 0. But Z will be correct.
  covariance.traverse_tree_marginal ();
  ll += mvn(design.nz) + 
      0.5 * covariance.logdet_marginal - 
      0.5 * covariance.Dyy;
  approximation.run (covariance.Z, covariance);
  if (!approximation.converged)
  {
    approximation.Z.t().print("Z");
    if (approximation.Z.has_nan())
      Rcpp::stop ("EP approximation failed; perhaps region of integration lacks sufficient density?");
    Rcpp::warning ("EP approximation did not converge in maximum number of iterations");
  }
  ll += approximation.Z(approximation.Z.n_elem-1);

  mat dm (covariance.dL.n_rows, covariance.Y.n_cols, arma::fill::zeros);//sloppy poppy
  covariance.traverse_reverse_marginal (approximation.dm, 0.5, -0.5);
  dL = constraints (Lambda, approximation.dL + covariance.dL, 
                    constraints.Lambda);
  dS = constraints (Sigma, approximation.dS + covariance.dS, 
                    constraints.Sigma);
  dm.rows(design.cont) -= covariance.dr; //sloppy ... is NEGATIVE correct???
  dm.rows(design.cat)  += approximation.dm;
  db = design (dm);

  return ll;
}

double Regression::mvn (const uword k)
{
  /* "constant" for MVN, not really necessary but useful in checking accuracy */
  return -0.5 * double(k) * log(2.0 * arma::datum::pi);
}

void Regression::options (const double tol, const uword maxiter, const bool fix_start)
{
  approximation.options (tol, maxiter, fix_start);
}

void Regression::profile (uword n, const vec& L, const vec& S, const std::vector<vec>& b)
{ /* perftools profiling */
  //ProfilerStart("./epee2_profile.log");
  //for (uword i=0; i<n; ++i)
  //{
  //  vec Ll = L; vec Ss = S; std::vector<vec> Bb = b;
  //  double l = loglikelihood (Ll, Ss, Bb);
  //}
  //ProfilerStop();
}


} /* end namespace epee */
