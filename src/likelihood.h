#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

#include <RcppArmadillo.h> 
#include "epee.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

namespace epee {

struct Regression 
{
    private:
      const Design          design;
      const Constraints     constraints;
            Covariance      covariance;
            Epee            approximation;

    public:
      vec   dL,
            dS;
      std::vector<vec> db;
            
      Regression (const Design&      des, 
                  const Constraints& cons) :
        design (des), 
        constraints (cons),
        covariance (design),
        approximation (covariance)
      {
        /* check validity */
        if (design.mx + design.mz != constraints.dim)
          Rcpp::stop ("Dimensions of constraints do not match design matrix");
      }

      double loglikelihood (const vec&, const vec&, const std::vector<vec>&);
      double mvn           (const uword);
      void   options       (const double, const uword, const bool);
      void   profile       (uword, const vec&, const vec&, const std::vector<vec>&);
};

} /* end namespace epee */
#endif
