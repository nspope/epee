#include <RcppArmadillo.h> 

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

/* the point is to determine the fastest way to reference non-contiguous views */
using namespace epee {

using cube = arma::cube;
using uvec = arma::uvec;
using uword = arma::uword;
using mat = arma::mat;

struct MyStruct
{
  private:
    cube A,
         B,
         C,
         D;

  public:
    MyStruct (uword m, uword n) :
      A (m, m, n, arma::fill::ones),
      B (m, m, n, arma::fill::ones),
      C (m, m, n, arma::fill::ones),
      D (m, m, n, arma::fill::ones)
    {}

    void doThang (uvec dz, uvec dx, uword n)
    {
      clear ();
      for (uword i=0; i<n; ++i)
      {
      auto Ar = A.slice(i)(dx,dx),
           Br = B.slice(i)(dx,dz),
           Cr = C.slice(i)(dz,dz),
           Dr = D.slice(i)(dz,dz);
      Ar += Br * Cr * Br.t();
      Cr  = Cr * Br.t() * Ar * Br;
      Dr  = Dr * Cr;
      Br += Ar * Br * Dr;
      }
    }

    void doThang2 (uvec dz, uvec dx, uword n)
    {
      clear ();
      for (uword i=0; i<n; ++i)
      {
      A.slice(i)(dx,dx) += B.slice(i)(dx,dz) * C.slice(i)(dz,dz) * B.slice(i)(dx,dz).t();
      C.slice(i)(dz,dz)  = C.slice(i)(dz,dz) * B.slice(i)(dx,dz).t() * A.slice(i)(dx,dx) * B.slice(i)(dx,dz);
      D.slice(i)(dz,dz)  = D.slice(i)(dz,dz) * C.slice(i)(dz,dz);
      B.slice(i)(dx,dz) += A.slice(i)(dx,dx) * B.slice(i)(dx,dz) * D.slice(i)(dz,dz);
      }
    }

    void doThang3 (uvec dz, uvec dx, uword n)
    {
    }

    void clear (void)
    {
      A.ones();
      B.ones();
      C.ones();
      D.ones();
    }
};

} /* end namespace epee */

RCPP_EXPOSED_CLASS (MyStruct)

RCPP_MODULE(RefTest) {
  using namespace Rcpp;

  class_<MyStruct>("RefTest")
    .constructor<uword,uword>()
    .method("doThang", &MyStruct::doThang)
    .method("doThang2", &MyStruct::doThang2)
    ;
}
