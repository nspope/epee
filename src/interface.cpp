#include "inputs.h"
#include "covariance.h"
#include "epee.h"
#include "likelihood.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

RCPP_EXPOSED_CLASS(epee::Covariance)
RCPP_EXPOSED_CLASS(epee::Epee)
RCPP_EXPOSED_CLASS(epee::Design)
RCPP_EXPOSED_CLASS(epee::Constraints)
RCPP_EXPOSED_CLASS(epee::Regression)

RCPP_MODULE(epee) {
  using namespace Rcpp;

  class_<epee::Covariance>("Covariance")
    .constructor<arma::umat, arma::vec, arma::mat, arma::uvec>()
    .field("Ayy", &epee::Covariance::Ayy)
    .field("Dyy", &epee::Covariance::Dyy)
    .field("logdet", &epee::Covariance::logdet)
    .field("logdet_marginal", &epee::Covariance::logdet_marginal)
    .field("Y", &epee::Covariance::Y)
    .field("Z", &epee::Covariance::Z)
    .field("diag", &epee::Covariance::diag)
    .field("dL", &epee::Covariance::dL)
    .field("dS", &epee::Covariance::dS)
    .field("dy", &epee::Covariance::dy)
    .field("dr", &epee::Covariance::dr)
    .method("traverse_tree_marginal", &epee::Covariance::traverse_tree_marginal)
    .method("traverse_tree", &epee::Covariance::traverse_tree)
    .method("traverse_reverse_marginal", &epee::Covariance::traverse_reverse_marginal)
    .method("traverse_reverse", &epee::Covariance::traverse_reverse)
    .method("parameters", &epee::Covariance::parameters)
    ;

  class_<epee::Epee>("Epee")
    .constructor<const epee::Covariance&>()
    .field("Z", &epee::Epee::Z)
    .field("dL", &epee::Epee::dL)
    .field("dS", &epee::Epee::dS)
    .field("dm", &epee::Epee::dm)
    .method("run", &epee::Epee::run)
    .method("options", &epee::Epee::options)
    ;

  class_<epee::Design>("Design")
    .constructor<arma::mat, std::vector<arma::mat>, arma::uvec, arma::umat, arma::vec>()
    .field_readonly("design", &epee::Design::design)
    .field_readonly("edge", &epee::Design::edge)
    .field_readonly("cat", &epee::Design::cat)
    .field_readonly("cont", &epee::Design::cont)
    .field_readonly("traits", &epee::Design::traits)
    .method("linear", &epee::Design::linear)
    .method("linear2", &epee::Design::linear2)
    ;

  class_<epee::Constraints>("Constraints")
    .constructor<arma::umat, arma::umat, arma::uword>()
    .field_readonly("Lambda", &epee::Constraints::Lambda)
    .field_readonly("Sigma", &epee::Constraints::Sigma)
    .method("find_nonzeros", &epee::Constraints::find_nonzeros)
    .method("apply", &epee::Constraints::apply)
    .method("apply_d", &epee::Constraints::apply_d)
    ;

 class_<epee::Regression>("Regression")
    .constructor<const epee::Design&, const epee::Constraints&>()
    .field("dL", &epee::Regression::dL)
    .field("dS", &epee::Regression::dS)
    .field("db", &epee::Regression::db)
    .method("loglikelihood", &epee::Regression::loglikelihood)
    .method("options", &epee::Regression::options)
    .method("profile", &epee::Regression::profile)
    ;
}
