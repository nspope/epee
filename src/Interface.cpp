#include <RcppArmadillo.h> 
#include "Epee.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace epee;

// [[Rcpp::export]]
void do_test (arma::mat lambda, arma::mat sigma, arma::mat Z, arma::umat E, arma::vec L, arma::vec H, arma::uvec ind)
{
	/* tree */
  Tree tree (Z, E, L, H);
	
	/* models */
  Brownian model (lambda, sigma);
	model.traits.print("model.traits");
  model.lambda.print("model.lambda");
  model.sigma.print("model.sigma");
  Brownian marg_model (model, ind);
	marg_model.traits.print("marg_model.traits");
  marg_model.lambda.print("marg_model.lambda");
  marg_model.sigma.print("marg_model.sigma");
	Brownian cond_model (model, ind);
	cond_model.traits.print("cond_model.traits");
  cond_model.lambda.print("cond_model.lambda");
  cond_model.sigma.print("cond_model.sigma");
  cond_model.delta.print("cond_model.delta");

  Marginal<Brownian> real_model (model);
  real_model.traits.print("real_model.traits");
  real_model.parameters.lambda.print("real_model.parameters.lambda");
  real_model.parameters.sigma.print("real_model.parameters.sigma");
  real_model.parameters.delta.print("real_model.parameters.delta");
  Marginal<Brownian> real_model2 = real_model.marginalize(ind);
  real_model2.traits.print("real_model2.traits");
  real_model2.parameters.lambda.print("real_model2.parameters.lambda");
  real_model2.parameters.sigma.print("real_model2.parameters.sigma");
  real_model2.parameters.delta.print("real_model2.parameters.delta");
  Conditional<Brownian> real_model3 = real_model2.condition();
  real_model3.traits.print("real_model3.traits");
  real_model3.parameters.lambda.print("real_model3.parameters.lambda");
  real_model3.parameters.sigma.print("real_model3.parameters.sigma");
  real_model3.parameters.delta.print("real_model3.parameters.delta");

	/* matrices */
  Precision<Marginal<Brownian>> prec (real_model2, tree);
  std::cout << "--------" << std::endl;
  for (auto i: prec.backbone.xxss)
		i.print();
  std::cout << "--------" << std::endl;
  for (auto i: prec.backbone.xx)
		i.print();
  std::cout << "--------" << std::endl;

	arma::mat D = arma::ones<mat>(ind.n_elem, tree.ntips);

  Precision<Marginal<Brownian>> prec_D (real_model2, tree, D);
  std::cout << "--------" << std::endl;
  for (auto i: prec_D.backbone.xxss)
		i.print();
  std::cout << "--------" << std::endl;
  for (auto i: prec_D.backbone.xx)
		i.print();
  std::cout << "--------" << std::endl;

	auto O = prec_D.linear(D);
	O.print("linear");
  std::cout << "--------" << std::endl;

	auto M = prec_D.quadratic(D);
	std::cout << "quadratic\n" << M << std::endl;
  std::cout << "--------" << std::endl;
	auto A = prec_D.diagonal();
	A.print("diagonal");
  std::cout << "--------" << std::endl;
	auto V = prec_D.determinant();
	std::cout << "determinant\n" << V << std::endl;
  std::cout << "--------" << std::endl;

 // auto prec_D_inv = prec_D.inverse(D);
 // auto meh = prec_D_inv.linear(D);
 // meh.print("covariance linear");
 // std::cout << "--------" << std::endl;
 // for(auto i: prec_D_inv.backbone.L)
 //   i.print();
 // std::cout << "--------" << std::endl;


	arma::mat D2 = arma::ones<mat>(Z.n_rows - ind.n_elem, tree.ntips);
  Precision<Conditional<Brownian>> prec_D2 (real_model3, tree, D2);
  std::cout << "--------" << std::endl;
  for (auto i: prec_D2.backbone.xx)
		i.print();
  std::cout << "--------" << std::endl;

	auto O2 = prec_D2.linear(D2);
	O2.print("linear");
  std::cout << "--------" << std::endl;
  prec_D2.linear.x.print("x");
  std::cout << "--------" << std::endl;
  prec_D2.linear.z.print("z");
  std::cout << "--------" << std::endl;

	auto M2 = prec_D2.quadratic(D2);
	std::cout << "quadratic\n" << M2 << std::endl;
  std::cout << "--------" << std::endl;
  for (auto i: prec_D2.quadratic.yx)
		i.print();
  std::cout << "--------" << std::endl;
  for (auto i: prec_D2.quadratic.yz)
		i.print();
  std::cout << "--------" << std::endl;

	auto A2 = prec_D2.diagonal();
	A2.print("diagonal");
  std::cout << "--------" << std::endl;
	auto V2 = prec_D2.determinant();
	std::cout << "determinant\n" << V2 << std::endl;
  std::cout << "--------" << std::endl;
}

class BrownianMotion
{
  private:
    Inputs                data;
    Regression<Brownian>  problem;

  public:
    BrownianMotion (mat o, uvec t, umat e, vec l, vec h, std::vector<mat> d) :
      data (o, t, e, l, h, d),
      problem (data)
    {}

    double loglikelihood (const mat& between, const mat& within, const std::vector<vec>& coef, const double tol=1e-4, const uword maxiter=100)
    {
      Brownian model (between, within);
      return problem.loglikelihood(model, coef, tol, maxiter);
    }
};

RCPP_EXPOSED_CLASS(BrownianMotion)

RCPP_MODULE(epee) {
  using namespace Rcpp;

  class_<BrownianMotion>("BrownianMotion")
    .constructor<mat, uvec, umat, vec, vec, std::vector<mat>>()
    .method("loglikelihood", &BrownianMotion::loglikelihood)
    ;
}

