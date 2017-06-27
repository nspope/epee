#include "Models.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

namespace epee
{

	uvec complement (const uvec& partial, const uvec& full) 
	{
		/* SetDiff wrapper; not the prettiest, but works OK */
		Rcpp::IntegerVector complement = 
			Rcpp::setdiff(Rcpp::IntegerVector(full.begin(), full.end()), 
					Rcpp::IntegerVector(partial.begin(), partial.end()));
		return arma::sort(arma::uvec(std::vector<uword>(complement.begin(), complement.end())));
	}

  uvec indices (uvec& subset, const uvec& full)
  {
    /* converts trait "names" to indices, and returns indices of complement */
    uvec other = complement(subset, full);
    uvec tmp (1);
    for(uword i=0; i<subset.n_elem; ++i)
    {
      tmp = find(full==subset.at(i), 1);
      subset.at(i) = tmp.at(0);
    }
    for(uword i=0; i<other.n_elem; ++i)
    {
      tmp = find(full==other.at(i), 1);
      other.at(i) = tmp.at(0); 
    }
    return other;
  }

  /* Brownian */
  template <>
  mat Marginal<Brownian>::edgelength_matrix   (const Tree::Edge& edge) const
  {
    mat lambda = parameters.lambda * edge.length();
    if (edge.is_tip)
      lambda += parameters.sigma;
    return lambda;
  }

  template <>
  vec Marginal<Brownian>::tip_variance           (const Tree::Rooter& edge) const
  {
    return edge.distal_height()*parameters.lambda.diag() + parameters.sigma.diag();
  }
  
  template <>
  mat Conditional<Brownian>::edgelength_matrix   (const Tree::Edge& edge) const
  {
    mat lambda = parameters.lambda * edge.length();
    if (edge.is_tip)
      lambda += parameters.sigma;
    return lambda;
  }
  
  template <>
  mat Conditional<Brownian>::covariance_matrix   (const Tree::Edge& edge) const
  {
    mat delta = parameters.delta * edge.length();
    if (edge.is_tip)
      delta += parameters.gamma;
    return delta;
  }

  template <>
  vec Conditional<Brownian>::tip_variance        (const Tree::Rooter& edge) const
  {
    return edge.distal_height()*parameters.lambda.diag() + parameters.sigma.diag(); //TODO this isn't done!
  }

} /* end namespace epee */
