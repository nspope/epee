#ifndef EPEE_H
#define EPEE_H

#include <RcppArmadillo.h> 
#include "Matrices.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

namespace epee {

  /* inputs including design matrices, observations, and phylogenetic tree */
  class Inputs 
  {
    public:
      std::vector<mat>        design;
      mat                     observations;
      uvec                    type;
      Tree                    tree;

      Inputs (mat o, uvec t, umat e, vec l, vec h, std::vector<mat> d) :
        observations (o),
        type (t),
        tree (o, e, l, h),
        design (d)
      {
        if (observations.n_rows != type.n_elem)
          Rcpp::stop ("Dimension mismatch: must define type for every trait");
        if (observations.n_rows != design.size())
          Rcpp::stop ("Dimension mismatch: must have same number of traits and design matrices");
        for (uword i=0; i<observations.n_rows; ++i)
          if (design.at(i).n_rows != observations.n_cols)
            Rcpp::stop ("Dimension mismatch: every tip must have a row in each design matrix");
      }

      inline rowvec mean (const uword i, const vec& coef)
      {
        return arma::trans(design.at(i) * coef);
      }

      inline rowvec residual (const uword i, const vec& coef)
      {
        return observations.row(i) - mean(i, coef);
      }
  };

  /* calculate the EP approximation */
  template <typename Form>
  struct Epee
  {
    public:
      bool            converged = false;
      uword           iter;
      double          Z = arma::datum::nan,
                      sigma_delta = 0.,
                      mu_delta = 0.;
      mat             kappa,
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
                      nu_tilde,
                      delta;

    Epee (const mat& z, const mat& eta, Precision<Form>& K, const double tol=1e-4, const uword maxiter=100) :
      kappa (size(z), arma::fill::zeros), 
      mu (size(z), arma::fill::zeros),
      sigma (size(z), arma::fill::zeros), 
      alpha (size(z), arma::fill::zeros),
      beta (size(z), arma::fill::zeros),
      tau_cavity (size(z), arma::fill::zeros),
      mu_cavity (size(z), arma::fill::zeros),
      tau_hat (size(z), arma::fill::zeros),
      mu_hat (size(z), arma::fill::zeros),
      Z_hat (size(z), arma::fill::ones),
      tau_tilde (size(z), arma::fill::zeros),
      nu_tilde (size(z), arma::fill::zeros),
      delta (size(z), arma::fill::zeros)
    {
      mu = eta;
      //mu.fill(0.1);
      kappa = K.linear(eta);
      sigma = K.inverse_diagonal();

      for(iter=0; iter<maxiter; ++iter)
      {
        update (z, eta, K);
        if (mu_delta < tol && sigma_delta < tol) //TODO not the best way to monitor convergence!
        {
          converged = true;
          break;
        }
      }

      Z = normalizing_constant (z, eta, K);
    }

    void update (const mat& z, const mat& eta, Precision<Form>& K)
    {
      /* normal EP updates */
      tau_cavity = 1.0/sigma - tau_tilde;
      mu_cavity = 1.0/tau_cavity % (mu/sigma - nu_tilde);
      alpha = mu_cavity/arma::sqrt(2.0/tau_cavity);
      Z_hat = 0.5*(1.0 - erf(-alpha)) - (1.0 - z);
      beta = 1.0/Z_hat % arma::sqrt(1.0/(2.0*arma::datum::pi*tau_cavity)) % arma::exp(-alpha%alpha);
      mu_hat = mu_cavity + beta; 
      tau_hat = 1.0/(mu_cavity%mu_cavity + 1.0/tau_cavity + mu_cavity%beta - mu_hat%mu_hat); 
      tau_tilde = tau_hat - tau_cavity;
      nu_tilde = mu_hat%tau_hat - mu_cavity%tau_cavity;

      /* update sigma */
      delta = sigma;
      beta = 1.0/tau_tilde;
      K.backbone(beta);
      sigma = K.diagonal();
      sigma = beta % (1.0 - beta % sigma);
      delta -= sigma;
      sigma_delta = arma::abs(delta).max();

      /* update mu */
      delta = mu;
      alpha = beta % (kappa + nu_tilde);
      mu = K.linear(alpha);
      mu = alpha - beta % mu;
      delta -= mu;
      mu_delta = arma::abs(delta).max();

    }

    double normalizing_constant (const mat& z, const mat& eta, Precision<Form>& K)
    {
      mat b = 1.0 + tau_tilde/tau_cavity;
      auto kdet = K.determinant();
      K.backbone(); // is there a way to be more judicious here I wonder?
      auto mkm = K.quadratic(eta);
      auto mukmu = K.quadratic(mu);

      double out = 0;
      out -= arma::accu(arma::log(tau_tilde));
      out += kdet; 
      out -= mkm;
      out += mukmu;
      out += arma::accu(tau_tilde%mu%mu + 2.0*arma::log(arma::abs(Z_hat)) + arma::log(b) + 
          (mu_cavity%mu_cavity%tau_tilde - 2.0*mu_cavity%nu_tilde - 
           nu_tilde%nu_tilde/tau_cavity)/b);
      out *= 0.5;

      return out;
    }

    inline mat erf (mat x)
    {
      for (uword i=0; i<x.n_rows; ++i)
        for (uword j=0; j<x.n_cols; ++j)
          x(i,j) = 2.0 * R::pnorm(x(i,j) * sqrt(2.0), 0.0, 1.0, true, false) - 1.0;
      return x;
    }
  };

  template <typename Model>
  struct Regression 
  {
    private:
      Inputs* const   data;
      const uvec      continuous,
                      categorical;
      mat             mean,
                      residuals;

    public:

      Regression (Inputs& d) :
        data (&d),
        continuous (arma::find(d.type == 0)),
        categorical (arma::find(d.type > 0)),
        mean (categorical.n_elem, d.tree.ntips),
        residuals (continuous.n_elem, d.tree.ntips)
      {}

      double loglikelihood (Model& model, const std::vector<vec>& coef, const double tol=1e-4, const uword maxiter=100)
      {
        Marginal<Model> full (model);

        if (data->observations.n_rows != full.dim || data->observations.n_rows != coef.size())
          Rcpp::stop ("Dimension mismatch: parameters don't match dimension of data");

        double likelihood = 0;
        if (continuous.n_elem > 0)
        {
          /* assemble matrices, etc. */
          Marginal<Model> continuous_only = full.marginalize (continuous);
          Precision<Marginal<Model>> K (continuous_only, data->tree);

          /* calculate marginal portion of likelihood */
          for (uword i=0; i<continuous.n_elem; ++i)
            residuals.row(i) = 
              data->residual(continuous(i), coef.at(continuous(i)));

          likelihood += mvn(continuous.n_elem * data->tree.ntips) -  //TODO this would have to change with missing data
            0.5 * K.quadratic(residuals) + 0.5 * K.determinant();

          if (categorical.n_elem > 0)
          {
            /* calculate conditional mean */
            Conditional<Model> categorical_only = continuous_only.condition();
            Precision<Conditional<Model>> W (categorical_only, data->tree);
            for (uword i=0; i<categorical.n_elem; ++i)
              mean.row(i) = 
                data->mean(categorical(i), coef.at(categorical(i)));
            residuals = K.linear (residuals);   
            mean += W.offdiagonal (residuals); 

            /* get EP approximation */
            Epee<Conditional<Model>> problem (data->observations.rows(categorical), mean, W, tol, maxiter);
            likelihood += problem.Z;
          }
        }
        else 
        {
          Marginal<Model> categorical_only = full.marginalize (categorical);
          Precision<Marginal<Model>> K (categorical_only, data->tree);
          for (uword i=0; i<categorical.n_elem; ++i)
            mean.row(i) = 
              data->mean(categorical(i), coef.at(categorical(i)));

          /* get EP approximation */
          Epee<Marginal<Model>> problem (data->observations, mean, K, tol, maxiter);
          likelihood += problem.Z;
        }
        return(likelihood);
      }

      inline double mvn (uword k)
      {
        return -0.5 * double(k) * log(2.0 * arma::datum::pi);
      }
  };
} /* end namespace epee */

#endif
