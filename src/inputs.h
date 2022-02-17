#ifndef INPUTS_H
#define INPUTS_H

#include <RcppArmadillo.h> 

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

namespace epee {

using submat  = arma::subview<double>;
using subvec  = arma::subview_col<double>;
using umat    = arma::umat;
using uvec    = arma::uvec;
using uword   = arma::uword;
using vec     = arma::vec;
using cube    = arma::cube;
using mat     = arma::mat;
using lvec    = std::vector<arma::vec>;
using luvec   = std::vector<arma::uvec>;
using lmat    = std::vector<arma::mat>;

struct Design
{
  public:
    const uvec   type,
                 cat,
                 cont;
    const uword  mx,
                 mz,
                 tips;
    const mat    traits;
    const umat   edge;
    const vec    edge_length;
    const uword  nx,
                 nz;
    const lmat   design;


    Design (mat  y,               /* trait values; each row is a trait          */
            lmat X,               /* list of design matrices                    */
            uvec typ,             /* vector denoting types of traits (1=binary) */
            umat edges,           /* edge matrix of phylogeny                   */
            vec  edge_lengths) :  /* vector of edge lengths                     */
      type        (typ),
      cat         (arma::find(type>0)),
      cont        (arma::find(type==0)),
      mx          (cat.n_elem),
      mz          (cont.n_elem),
      tips        (y.n_cols),
      traits      (missing_data_transfer(y, X)),
      design      (missing_data_transfer(X, y)),
      edge        (edges.min() == 1 ? edges - 1 : edges),
      edge_length (edge_lengths),
      nx          (arma::numel(arma::find_finite(traits.rows(cat)))),
      nz          (arma::numel(arma::find_finite(traits.rows(cont))))
    {
      /* check validity */
      if (edge_lengths.n_elem != edges.n_rows)
        Rcpp::stop("Must have a length for each edge");
      if (type.n_elem != traits.n_rows)
        Rcpp::stop("Must have a type for each trait");
      if (design.size() != mx+mz)
        Rcpp::stop("Must have design matrix for each trait");
      for (auto& x : design)
        if (x.n_rows != tips)
          Rcpp::stop("Design matrices must have a row for each tip in phylogeny");
        else if (x.n_cols < 1)
          Rcpp::stop("Design matrices should include at least a single column (intercept)");
    }

    /* used for construction */
    mat missing_data_transfer (mat y, lmat& X)
    { /* make sure missing values in design matrix are transferred to traits */
      for (uword i=0; i<mx+mz; ++i)
        for (uword j=0; j<tips; ++j)
          if (X.at(i).row(j).has_nan())
          {
            y(i,j) = arma::datum::nan;
            X.at(i).row(j).fill(arma::datum::nan); // (for now) setting row to NaN in design
          }
      return y;
    }

    lmat missing_data_transfer (lmat X, const mat& y)
    { /* make sure missing values in traits are transferred to design matrix */
      for (uword i=0; i<mx+mz; ++i)
        for (uword j=0; j<tips; ++j)
          if (std::isnan(y(i,j)))
            X.at(i).row(j).fill(arma::datum::nan);
      return X;
    }
    
    /* const */
    mat operator() (const lvec& b) const
    { /* take coefficients and map to matrix of means */
      mat mean (mx + mz, tips);
      for (uword i=0; i<mx+mz; ++i)
        mean.row(i) = b.at(i).t() * design.at(i).t();
      return mean;
    }

    lvec operator() (const mat& mean) const
    { /* adjoint takes matrix of means and maps to coefficients */
      std::vector<vec> b (mx+mz);
      for (uword i=0; i<mx+mz; ++i)
      { /* TODO sloppy poppy */
        mat dd = design.at(i).t();
        mat mn = mean.row(i).t();
        dd.replace(arma::datum::nan, 0.);
        mn.replace(arma::datum::nan, 0.);
        b.at(i) = dd * mn;
      }
      return b;
    }
    // !!!!
    // duplication just for testing, not sure how to export operators to R
    mat linear (const lvec& b) const
    { /* take coefficients and map to matrix of means */
      mat mean (mx + mz, tips);
      for (uword i=0; i<mx+mz; ++i)
        mean.row(i) = b.at(i).t() * design.at(i).t();//TODO missing covariate data? 
      return mean;
    }

    lvec linear2 (const mat& mean) const
    { /* adjoint takes matrix of means and maps to coefficients */
      std::vector<vec> b (mx+mz);
      for (uword i=0; i<mx+mz; ++i)
      { /* TODO sloppy poppy */
        mat dd = design.at(i).t();
        mat mn = mean.row(i).t();
        dd.replace(arma::datum::nan, 0.);
        mn.replace(arma::datum::nan, 0.);
        //b.at(i) = design.replace(arma::datum::nan, 0.).at(i).t() * mean.replace(arma::datum::nan, 0.).row(i).t(); 
        b.at(i) = dd * mn;
      }
      return b;
    }
};

struct Constraints
{
  /* container for zero-constraints on covariance matrix */
  public:
    const uword dim;
    const umat  Lambda,
                Sigma;

    Constraints (umat L, umat S, uword d) :
      dim    (d),
      Lambda (sort(L)),
      Sigma  (sort(S))
    {
    }

    umat sort (umat M)
    { /* sort and remove duplicates, also check validity */
      if (M.n_elem > 0)
      {
        if (M.max() >= dim || M.n_cols != 2)
          Rcpp::stop ("Constraints do not match dimension");
        //for (uword i=0; i<M.n_rows; ++i)
        //  if (M(i,0) == M(i,1))
        //    Rcpp::stop ("Zero constraints not allowed on diagonal");
        M = arma::sort(M, /* largest index comes second in each row */
            "descend", 
            1);               
        uvec Mi = arma::unique(
            arma::sub2ind(arma::size(dim, dim), 
              M.t()));  /* unique linear indices */
        M = arma::ind2sub(arma::size(dim, dim), 
            Mi).t();
      }
      return M;
    }

    uvec find_nonzeros (const umat& M) const
    { /* find nonzero elements (of lower triangle) */
      umat U = arma::trimatl(arma::ones<umat>(dim,dim));
      for (uword i=0; i<M.n_rows; ++i)
        U(M(i,0),M(i,1)) = 0;
      return arma::find(U);
    }

    mat operator() (const vec&  v,          /* free Choleski elements as a vector */
                    const umat& cons) const
    { /* create Choleski factor from free elements + constraints */
      uword i, j;
      mat V = arma::zeros<mat>(dim, dim);
      uvec nonzeros = find_nonzeros(cons);
      V.elem(nonzeros) = v;
      V.diag() = arma::exp(V.diag());
      for (uword z=0; z<cons.n_rows; ++z)//TODO use iterator
        if (cons(z,1) > 0)
        {
          i = cons(z,0); j = cons(z,1);// k = j - 1;
          if (i==j) 
            V.row(j) /= sqrt(arma::dot(V.row(j), V.row(j))); 
          else
            V(i,j) = 
              -arma::sum(V.row(i).head(j) %
                V.row(j).head(j) /
                V(j,j));
        }
      return V;
    }

    vec operator() (      mat   V,          /* Choleski factor as a matrix */
                          mat   dV,         /* Derivatives of covariance matrix */
                    const umat& cons) const
    { /* map derivatives back to free elements of Choleski factor */
      uword i, j;
      uvec nonzeros = find_nonzeros(cons);
      dV = (dV + dV.t()) * V;
      dV = arma::trimatl(dV);
      for (uword z=cons.n_rows; z>0;--z)//TODO use iterator
        if (cons(z-1,1) > 0)
        {
          i = cons(z-1,0); j = cons(z-1,1);// k = j - 1;
          if (i==j)
          { 
            dV.row(j).head(j) = dV.row(j).head(j+1) *
              (arma::eye(j+1,j) - V.row(j).head(j+1).t() * V.row(j).head(j)) * V(j,j);
            V.row(j).head(j+1) /= V(j,j); /* otherwise will fail ? */
          }
          else
          {
            dV.row(j).head(j) -= dV(i,j) * V.row(i).head(j)/V(j,j);
            dV.row(i).head(j) -= dV(i,j) * V.row(j).head(j)/V(j,j);
            dV(j,j)           += dV(i,j) * arma::sum(V.row(i).head(j) % 
                                  V.row(j).head(j))/(V(j,j)*V(j,j));
          }
        }
      dV.diag() %= V.diag();
      return dV.elem(nonzeros);
    }

    // duplication just for testing, not sure how to export operators to R
    mat apply (const vec&  v,     
               const umat& cons) const 
    {
      return operator()(v, cons);
    }

    vec apply_d (      mat   V, 
                       mat   dV,
                 const umat& cons) const
    { 
      return operator()(V, dV, cons);
    }

};

} /* end namespace epee */
#endif

