#ifndef MODELS_H
#define MODELS_H

#include <RcppArmadillo.h> 
#include "Tree.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

namespace epee 
{
	uvec complement (const uvec&, const uvec&); /* set complement */
  uvec indices (uvec&, const uvec&); /* convert labels to indices, return vector with non-matching indices */

  /* base types of models */

  template <typename Model> struct Conditional;

  template <typename Model> 
    struct Marginal 
    {
      public:
        Marginal<Model>* const          parent;
        const uvec                      traits;
        const uword                     dim;
        const Model                     parameters;

        Marginal (Model& m) : 
          parent (this),
          traits (m.traits),
          dim (traits.n_elem),
          parameters (m)
        {}

        Marginal (Marginal<Model>& m, const uvec& i) : 
          parent (&m),
          traits (i),
          dim (traits.n_elem),
          parameters (Model(m.parameters, i))
        {}

        mat                 edgelength_matrix  (const Tree::Edge&) const;
        vec                 tip_variance       (const Tree::Rooter&) const;

        Marginal<Model>     marginalize        (const uvec& i)
        {
          return Marginal<Model>(*this, i);
        }

        Conditional<Model>  condition          (void)
        {
          return Conditional<Model>(*this);
        }
    };

  template <typename Model>
    struct Conditional
    {
      public:
          Marginal<Model>* const        parent;
          const uvec                    traits;
          const uword                   dim;
          const Model                   parameters;

          Conditional (Marginal<Model>& m) :
            parent (&m),
            traits (complement(m.traits, m.parent->traits)),
            dim (traits.n_elem),
            parameters (Model(m.parent->parameters, traits))
          {}

          mat                 edgelength_matrix  (const Tree::Edge&) const;
          mat                 covariance_matrix  (const Tree::Edge&) const;
          vec                 tip_variance       (const Tree::Rooter&) const;
    };

  /* 
     Brownian motion. Define:
     L = between species trait covariance,
     S = within species trait covariance,
     V = phylogenetic covariance matrix,
     Then this model has precision matrix K, where:
     K^{-1} = V \otimes L + I \otimes S
   */
  struct Brownian
  {
    public:
      const uvec  traits;
      const mat   lambda, 
                  sigma, 
                  delta,
                  gamma;

      /* used for initialization */
      mat init_diagonal (uvec i, const uvec& t, const mat& l)
      {
        uvec j = indices(i, t);
        return l.submat(i, i);
      }

      mat init_offdiagonal (uvec i, const uvec& t, const mat& l)
      {
        uvec j = indices(i, t);
        return l.submat(i, j);
      }

      Brownian (const mat& l, const mat& s) : /* initial cnstr */
        traits (arma::regspace<uvec>(0, 1, l.n_rows-1)),
        lambda (l),
        sigma (s),
        delta (arma::zeros<mat>(0, 0)),
        gamma (arma::zeros<mat>(0, 0))
      {}

      Brownian (const Brownian& b, const uvec& i): /* marginal or conditional */
        traits (i),
        lambda (init_diagonal(i, b.traits, b.lambda)),
        sigma (init_diagonal(i, b.traits, b.sigma)),
        delta (init_offdiagonal(i, b.traits, b.lambda)),
        gamma (init_offdiagonal(i, b.traits, b.sigma))
      {} // TODO: throw error if "conditioning" with full set
  };

} /* namepsace epee */

#endif
