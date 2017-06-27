#ifndef MATRICES_H
#define MATRICES_H

#include <RcppArmadillo.h> 
#include "Tree.h"
#include "Models.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

namespace epee
{
  inline mat make_symmetric (const mat& m)
  {
    return m.t() + m;
  }

  template <typename T>
  inline void clean_vector (std::vector<T>& v)
  {
    for (auto &i : v) i.zeros();
  }

  template <>
  inline void clean_vector (std::vector<double>& y)
  {
    std::fill (y.begin(), y.end(), 0.0);
  }

  template <typename E, typename T, typename ... Args>
    void do_recursion (T& object, E&& edge, Args... args)
    {
      object.clean ();
      for (; edge<=edge.end(); edge++)
      {
        if (edge.is_tip)
          object.tip (edge, args...);
        else
          object.node (edge, args...);
        object.increment (edge, args...);
      }
    }

 // template <typename Form>
 //   struct Covariance 
 //   {}; // TODO: throw error as should *always* be using specialized case?

  template <typename Form>
    struct Precision
    {}; // TODO: throw error as should *always* be using specialized case?

  /* Marginal case */
 // template <>
 //   template <typename Model>
 //   struct Covariance<Marginal<Model>>
 //   {
 //     private:
 //       Marginal<Model>* const     model;
 //       Tree* const                tree;

 //     public:
 //       struct Backbone
 //       {
 //         private:
 //           Covariance<Marginal<Model>>* const   W;

 //         public:
 //           std::vector<mat> L;

 //           Backbone                  (Covariance<Marginal<Model>>& p) :
 //             W (&p),
 //             L (p.tree->nedges, arma::zeros<mat>(p.model->dim, p.model->dim))
 //           {}

 //           void tip (const Tree::Edge& edge)
 //           {
 //             L.at(edge.position) += W->model->edgelength_matrix(edge);
 //           }

 //           void tip (const Tree::Edge& edge, const mat& d)
 //           {
 //             L.at(edge.position) = arma::diagmat(d.col(edge.descendant));
 //             tip (edge);
 //           }

 //           void node (const Tree::Edge& edge)
 //           {
 //             auto k = double(edge.descent().n_elem);
 //             L.at(edge.position) = W->model->edgelength_matrix(edge) * k;
 //           }

 //           void node (const Tree::Edge& edge, const mat& d)
 //           {
 //             node(edge);
 //           }

 //           void increment (const Tree::Edge& edge)
 //           {}

 //           void increment (const Tree::Edge& edge, const mat& d)
 //           {}

 //           inline void clean (void)
 //           {
 //           }

 //           void operator() (void)
 //           {
 // 						do_recursion<Tree::Edge, Backbone>(*this, W->tree->begin());
 //           }

 //           void operator() (const mat& d)
 //           {
 // 						do_recursion<Tree::Edge, Backbone, const mat&>(*this, W->tree->begin(), d);
 //           }
 //       };

 //       struct Linear
 //       {
 //         private:
 //           Covariance<Marginal<Model>>* const   W;

 //         public:
 //           mat       o;

 //           Linear                    (Covariance<Marginal<Model>>& p) :
 //             W (&p),
 //             o (p.model->dim, p.tree->ntips, arma::fill::zeros)
 //           {}

 //           void tip (const Tree::Edge& edge, const mat& y)
 //           {
 //             o.col(edge.descendant) =
 //               W->backbone.L.at(edge.position) *
 //               y.col(edge.descendant);
 //           }

 //           void node (const Tree::Edge& edge, const mat& y)
 //           {
 //             for(auto i : edge.descent())
 //               o.col(i) += (W->backbone.L.at(edge.position) * y.col(i));
 //           }

 //           void increment (const Tree::Edge& edge, const mat& y) 
 //           {}

 //           inline void clean (void)
 //           {
 //           }

 //           mat operator() (const mat& y)
 //           {
 // 						do_recursion<Tree::Edge, Linear, const mat&>(*this, W->tree->begin(), y);
 // 						return o;
 //           }
 //       };

 //       Backbone backbone;
 //       Linear linear;

 // 			Covariance (Marginal<Model>& m, Tree& t) : /* cnstr */
 // 				model (&m),
 // 				tree (&t),
 //         backbone (*this),
 // 				linear (*this)
 // 			{
 //         backbone();
 //       }

 // 			Covariance (Marginal<Model>& m, Tree& t, const mat& d) : /* cnstr */
 // 				model (&m),
 // 				tree (&t),
 //         backbone (*this),
 // 				linear (*this)
 // 			{
 //         backbone(d);
 //       }
 //   };

  template <>
  template <typename Model>
    struct Precision<Marginal<Model>> 
    {
      /* precision matrix arising from phylogenetic tree. Some notation:
         K is precision matrix
         X is design matrix, mapping observations to traits
       */
      private:
        Marginal<Model>* const     model;
        Tree* const                tree;

      public:
        struct Backbone 
        {
          /* X'KX */ 
          private:
            Precision<Marginal<Model>>* const    K;

          public:
            std::vector<mat>           xx,
                                       xxss;

            Backbone                   (Precision<Marginal<Model>>& p) :
              K (&p),
              xx (p.tree->nnodes, arma::zeros<mat>(p.model->dim, p.model->dim)),
              xxss (p.tree->nedges, arma::zeros<mat>(p.model->dim, p.model->dim))
            {}

						inline void tip            (const Tree::Edge& edge)
						{
              xx.at(edge.descendant) = arma::inv(xx.at(edge.descendant) + 
                K->model->edgelength_matrix(edge));
							xxss.at(edge.position) = xx.at(edge.descendant);
						}

						inline void tip            (const Tree::Edge& edge, const mat& d)
						{
              xx.at(edge.descendant) = arma::diagmat(d.col(edge.descendant));
              tip(edge);
						}

						inline void node           (const Tree::Edge& edge)
						{
							xxss.at(edge.position) = 
								arma::inv(arma::eye(size(xx.at(edge.descendant))) + 
										K->model->edgelength_matrix(edge)*xx.at(edge.descendant)); 
							xx.at(edge.descendant) *= xxss.at(edge.position);
						}

						inline void node           (const Tree::Edge& edge, const mat& d)
						{
							node (edge);
						}

						inline void increment      (const Tree::Edge& edge)
						{
							xx.at(edge.ancestor) += xx.at(edge.descendant);
						}

						inline void increment      (const Tree::Edge& edge, const mat& d)
						{
							increment (edge);
						}

            inline void clean          (void)
            {
              clean_vector<mat> (xx);
            }

						inline mat operator()      (void)
						{
							do_recursion <Tree::Edge, Backbone> (*this, K->tree->begin());
							return xx.at(K->tree->root);
						}

						inline mat operator()      (const mat& d)
						{
							do_recursion <Tree::Edge, Backbone, const mat&> (*this, K->tree->begin(), d);
							return xx.at(K->tree->root);
						}
        };

        struct Linear
        {
          /* Ky */
          private:
            Precision<Marginal<Model>>* const   K;

          public:
            cube                      x;
            std::vector<vec>          xy;

            Linear                    (Precision<Marginal<Model>>& p) :
              K (&p),
              x (p.model->dim, p.model->dim, p.tree->ntips),
              xy (p.tree->nnodes, arma::zeros<vec>(p.model->dim))
            {}

						inline void tip                    (const Tree::Edge& edge, mat& y)
						{
							x.slice(edge.descendant) = K->backbone.xxss.at(edge.position);
							xy.at(edge.descendant) = x.slice(edge.descendant)*y.col(edge.descendant);
							y.col(edge.descendant) = xy.at(edge.descendant);
						}

						inline void node                   (const Tree::Edge& edge, mat& y)
						{                                                                               /* let A = (I + TX'K^{-1}X)^{-1} */
							for (uword i : edge.descent())
							{
								y.col(i) -= x.slice(i) * K->backbone.xxss.at(edge.position) *               /* K^{-1}y - K^{-1}XATX'K^{-1}y */
									K->model->edgelength_matrix(edge) * xy.at(edge.descendant);  
								x.slice(i) *= K->backbone.xxss.at(edge.position);                           /* K^{-1}XA */
							}
							xy.at(edge.descendant) = K->backbone.xxss.at(edge.position).t() * 
								xy.at(edge.descendant);                                                     /* A'X'K^{-1}y */
						}

						inline void increment              (const Tree::Edge& edge, mat& y)
						{
							xy.at(edge.ancestor) += xy.at(edge.descendant);
						}

            inline void clean                  (void)
            {
              clean_vector<vec> (xy);
            }

						inline mat operator()              (mat y)
						{
							do_recursion<Tree::Edge, Linear, mat&>(*this, K->tree->begin(), y);
							return y;
						}
        };

        struct Quadratic
        {
          /* y'Ky */
          private:
            Precision<Marginal<Model>>* const    K;

          public:
            std::vector<vec>           xy;
            std::vector<double>        yy;

            Quadratic                  (Precision<Marginal<Model>>& p) :
              K (&p),
              xy (p.tree->nnodes, arma::zeros<vec>(p.model->dim)),
              yy (p.tree->nnodes, 0.0)
            {}

						inline void tip            (const Tree::Edge& edge, const mat& y)
						{
							xy.at(edge.descendant) = 
								K->backbone.xxss.at(edge.position) * y.col(edge.descendant);
							yy.at(edge.descendant) = 
								arma::dot(y.col(edge.descendant), xy.at(edge.descendant));
						}

						inline void node           (const Tree::Edge& edge, const mat& y)
						{
							yy.at(edge.descendant) -= 
								arma::dot(xy.at(edge.descendant), K->backbone.xxss.at(edge.position) *
										K->model->edgelength_matrix(edge) * xy.at(edge.descendant));
							xy.at(edge.descendant) = K->backbone.xxss.at(edge.position).t() * 
								xy.at(edge.descendant);
						}

						inline void increment      (const Tree::Edge& edge, const mat& y)
						{
							xy.at(edge.ancestor) += xy.at(edge.descendant);
							yy.at(edge.ancestor) += yy.at(edge.descendant);
						}

            inline void clean          (void)
            {
              clean_vector<vec> (xy);
              clean_vector<double> (yy);
            }

						inline double operator()   (const mat& y)
						{
							do_recursion<Tree::Edge, Quadratic, const mat&>(*this, K->tree->begin(), y);
							return yy.at(K->tree->root);
						}
        };

        struct Diagonal
        {
          /* diag(K) */
          private:
            Precision<Marginal<Model>>* const   K;

          public:
            mat                       diag,
                                      xy,
                                      yy;

            Diagonal                  (Precision<Marginal<Model>>& p) :
              K (&p),
              diag (p.model->dim, p.tree->ntips, arma::fill::zeros),
              xy (p.model->dim, p.model->dim),
              yy (p.model->dim, p.model->dim)
            {}

						inline void tip           (const Tree::Rooter& edge)
						{
							xy = yy = K->backbone.xxss.at(edge.position);
						}

						inline void node          (const Tree::Rooter& edge)
						{
							yy -= xy.t() * K->backbone.xxss.at(edge.position) * 
								K->model->edgelength_matrix(edge) * xy;
							xy = K->backbone.xxss.at(edge.position).t() * xy;
						}

						inline void increment     (const Tree::Rooter& edge)
						{
							if (edge.ancestor == K->tree->root)
								diag.col(edge.tip) = yy.diag();
						}

            inline void clean         (void)
            {
            }

						inline mat operator()     (void)
						{
							for(uword i=0; i < K->tree->ntips; ++i)
								do_recursion<Tree::Rooter, Diagonal>(*this, K->tree->begin(i));
							return diag;
						}
        };

        struct Determinant
        {
          /* |K| */
          private:
            Precision<Marginal<Model>>* const     K;

          public:
            std::vector<double>         logdet; 

            Determinant                 (Precision<Marginal<Model>>& p) :
              K (&p),
              logdet (p.tree->nnodes, 0.0)
            {}

						inline void tip             (const Tree::Edge& edge)
						{
							logdet.at(edge.descendant) +=                               
								log(arma::det(K->backbone.xxss.at(edge.position))); 
						}

						inline void node            (const Tree::Edge& edge)
						{
							logdet.at(edge.descendant) +=                               /* -log|I + TX'K^{-1}X| */
								log(arma::det(K->backbone.xxss.at(edge.position))); 
						}

						inline void increment       (const Tree::Edge& edge)
						{
							logdet.at(edge.ancestor) += logdet.at(edge.descendant);
						}

            inline void clean           (void)
            {
              clean_vector<double> (logdet);
            }

						inline double operator()    (void)
						{
							do_recursion<Tree::Edge, Determinant>(*this, K->tree->begin());
							return logdet.at(K->tree->root);
						}
				};

        struct InverseDiagonal
        {
          private:
            Precision<Marginal<Model>>* const   K;
            mat                                 sigma;

          public:
            InverseDiagonal                    (Precision<Marginal<Model>>& p) :
              K (&p),
              sigma (p.model->dim, p.tree->ntips, arma::fill::zeros)
            {}

            mat operator() (void)
            {
              for(uword i=0; i<K->tree->ntips; ++i)
                sigma.col(i) = K->model->tip_variance(K->tree->begin(i));
              return sigma;
            }
        };

        Backbone            backbone;
        Linear              linear;
        Quadratic           quadratic;
        Diagonal            diagonal;
        Determinant         determinant;
        InverseDiagonal     inverse_diagonal;

				Precision (Marginal<Model>& m, Tree& t) : /* cnstr */
					model (&m),
					tree (&t),
					backbone (*this),
					linear (*this),
					quadratic (*this),
					diagonal (*this),
					determinant (*this),
          inverse_diagonal (*this)
  			{
					backbone ();
  			}

				Precision (Marginal<Model>& m, Tree& t, const mat& d) : /* cnstr */
					model (&m),
					tree (&t),
					backbone (*this),
					linear (*this),
					quadratic (*this),
					diagonal (*this),
					determinant (*this),
          inverse_diagonal (*this)
				{
					backbone (d);
				}
    };


 /* Conditional case */ 
  template <>
  template <typename Model>
    struct Precision<Conditional<Model>> 
    {
      /* precision matrix arising from phylogenetic tree. Some notation:
         K is precision matrix
         X is design matrix, mapping observations to traits
         This particular form is for the inverse of a conditional covariance matrix; ie.
         K^{-1} = W^*_{i,i} = W_{i,i} - W_{i,j} W_{j,j}^{-1} W_{j,i}
       */

      private:
        Conditional<Model>* const  model;
        Tree* const                tree;

      public:
        struct Backbone 
        {
          /* X'KX */ 
          private:
            Precision<Conditional<Model>>* const     K;
            uword                                    dim_x,
                                                     dim_z;
            mat                                      Lx, 
                                                     Lz, 
                                                     S, 
                                                     Ix, 
                                                     Iz,
                                                     B2,
                                                     C2,
                                                     B3,
                                                     C3;
          public:
            /* actual recursion */
            std::vector<mat>           xx,      // X'KX
                                       zz,      // Z'MZ
                                       xz,      // X'KDMZ
                                       zxz;     // Z'MD'KDMZ
            /* storage */
            std::vector<mat>           A1,
                                       B5,
                                       D2,
                                       D3,
                                       D5,
                                       BD2,
                                       BD3,
                                       CD2,
                                       CD3,
                                       xx1,
                                       xx3,
                                       zxz2,
                                       zz3, // can't this just be zz4?
                                       zz4,
                                       xz1,
                                       xz2,
                                       xz3;
           

            Backbone                   (Precision<Conditional<Model>>& p) :
              K (&p),
              dim_x (p.model->dim),
              dim_z (p.model->parent->dim),
              xx (p.tree->nnodes, arma::zeros<mat>(dim_x, dim_x)),
              zz (p.tree->nnodes, arma::zeros<mat>(dim_z, dim_z)),
              xz (p.tree->nnodes, arma::zeros<mat>(dim_x, dim_z)), 
              zxz (p.tree->nnodes, arma::zeros<mat>(dim_z, dim_z)),
              A1 (p.tree->nedges, arma::zeros<mat>(dim_x, dim_x)), 
              B5 (p.tree->nedges, arma::zeros<mat>(dim_z, dim_z)), 
              D2 (p.tree->nedges, arma::zeros<mat>(dim_z, dim_z)), 
              D3 (p.tree->nedges, arma::zeros<mat>(dim_z, dim_z)), 
              D5 (p.tree->nedges, arma::zeros<mat>(dim_z, dim_z)), 
              BD2 (p.tree->nedges, arma::zeros<mat>(dim_z, dim_z)),
              BD3 (p.tree->nedges, arma::zeros<mat>(dim_z, dim_x)),
              CD2 (p.tree->nedges, arma::zeros<mat>(dim_x, dim_z)), 
              CD3 (p.tree->nedges, arma::zeros<mat>(dim_x, dim_x)),
              xx1 (p.tree->nedges, arma::zeros<mat>(dim_x, dim_x)), 
              xx3 (p.tree->nedges, arma::zeros<mat>(dim_x, dim_x)), 
              xz1 (p.tree->nedges, arma::zeros<mat>(dim_x, dim_z)), 
              xz2 (p.tree->nedges, arma::zeros<mat>(dim_x, dim_z)), 
              xz3 (p.tree->nedges, arma::zeros<mat>(dim_x, dim_z)),
              zxz2 (p.tree->nedges, arma::zeros<mat>(dim_z, dim_z)),
              zz3 (p.tree->nedges, arma::zeros<mat>(dim_z, dim_z)),
              zz4 (p.tree->nedges, arma::zeros<mat>(dim_z, dim_z)),
              Lx (arma::zeros<mat>(dim_x, dim_x)),
              Lz (arma::zeros<mat>(dim_z, dim_z)),
              S (arma::zeros<mat>(dim_x, dim_z)),
              Ix (arma::eye<mat>(dim_x, dim_x)),
              Iz (arma::eye<mat>(dim_z, dim_z)),
              B2 (arma::zeros<mat>(dim_z, dim_z)),
              B3 (arma::zeros<mat>(dim_z, dim_z)),
              C2 (arma::zeros<mat>(dim_x, dim_z)),
              C3 (arma::zeros<mat>(dim_x, dim_z))
            {}

						inline void tip            (const Tree::Edge& edge)
						{ 
              S = K->model->covariance_matrix(edge);
              zz.at(edge.descendant) = arma::inv(K->model->parent->edgelength_matrix(edge));
							xx.at(edge.descendant) = 
                arma::inv(xx.at(edge.descendant) + K->model->edgelength_matrix(edge) -
                    S * zz.at(edge.descendant) * S.t());
              xz.at(edge.descendant) = xx.at(edge.descendant) * 
                S * zz.at(edge.descendant); 
              zxz.at(edge.descendant) = zz.at(edge.descendant) * 
                S.t() * xz.at(edge.descendant);
						}

						inline void tip            (const Tree::Edge& edge, const mat& d)
						{ 
							xx.at(edge.descendant) = arma::diagmat(d.col(edge.descendant));
              tip (edge);
						}

						inline void node           (const Tree::Edge& edge)
            { 
              S = K->model->covariance_matrix(edge);
              Lz = K->model->parent->edgelength_matrix(edge);
              Lx = K->model->edgelength_matrix(edge);

              /* step 1; update K - KX(I + LX'KX)^{-1}X'K */
              xx1.at(edge.position) = xx.at(edge.descendant);
              xz1.at(edge.position) = xz.at(edge.descendant);

              A1.at(edge.position) = 
                arma::inv(Ix + Lx * xx.at(edge.descendant));
              zxz.at(edge.descendant) -=                                                        // Z'MD'KDMZ - Z'MD'KX(I + LX'KX)^{-1}TX'KDMZ
                xz.at(edge.descendant).t() * A1.at(edge.position) * 
                Lx * xz.at(edge.descendant);
              xz.at(edge.descendant) = A1.at(edge.position).t() * xz.at(edge.descendant);       // (I + LX'KX)^{-T}X'KD'MZ
              xx.at(edge.descendant) *= A1.at(edge.position);                                   // X'KX(I + LX'KX)^{-1}
              
              /* step 2; update K - K(-XS + DMZL')(I + LZ'MZ - S'X'KDMZ + LZ'MD'KDMZ)^{-T} Z'MD'K   */
              xz2.at(edge.position) = xz.at(edge.descendant);
              zxz2.at(edge.position) = zxz.at(edge.descendant);

              B2 = -S.t() * xz.at(edge.descendant) + Lz * zxz.at(edge.descendant);
              C2 = -xx.at(edge.descendant) * S + xz.at(edge.descendant) * Lz; 
              D2.at(edge.position) = arma::inv(Iz + Lz * zz.at(edge.descendant) + B2);
              BD2.at(edge.position) = arma::trans(D2.at(edge.position) * B2);
              CD2.at(edge.position) = C2 * D2.at(edge.position).t();

              xx.at(edge.descendant) -= CD2.at(edge.position) * xz.at(edge.descendant).t();
              mat zx = xz.at(edge.descendant).t();
              xz.at(edge.descendant) -= CD2.at(edge.position) * zxz.at(edge.descendant);
              zxz.at(edge.descendant) = (Iz - BD2.at(edge.position)) * zxz.at(edge.descendant);

              /* step 3; update K - K(-DMZ - XSZ'MZ)(I + LZ'MZ - S'X'KDMZ - S'X'KXSZ'MZ)^{-1} X'K  */
              xx3.at(edge.position) = xx.at(edge.descendant);
              xz3.at(edge.position) = xz.at(edge.descendant);
              zz3.at(edge.position) = zz.at(edge.descendant);

              B3 = -zxz.at(edge.descendant) - 
                (Iz - BD2.at(edge.position)) * zx * S * zz.at(edge.descendant);
              C3 = -xz.at(edge.descendant) - xx.at(edge.descendant) * S * zz.at(edge.descendant);
              D3.at(edge.position) = arma::inv(Iz + Lz * zz.at(edge.descendant) + S.t() * C3);
              CD3.at(edge.position) = C3 * D3.at(edge.position) * S.t();
              BD3.at(edge.position) = B3 * D3.at(edge.position) * S.t();

              xx.at(edge.descendant) = (Ix - CD3.at(edge.position)) * xx.at(edge.descendant);
              zxz.at(edge.descendant) -= BD3.at(edge.position) * xz.at(edge.descendant);
              xz.at(edge.descendant) = (Ix - CD3.at(edge.position)) * xz.at(edge.descendant);

              /* step 4; update D + XSZ' */
              zz4.at(edge.position) = zz.at(edge.descendant);

              zxz.at(edge.descendant) += 
                make_symmetric(zz.at(edge.descendant) * S.t() * xz.at(edge.descendant)) +
                zz.at(edge.descendant) * S.t() * xx.at(edge.descendant) * S * zz.at(edge.descendant);
              xz.at(edge.descendant) +=
                xx.at(edge.descendant) * S * zz.at(edge.descendant);

              /* step 5: update M - MZ(I + LZ'MZ)^{-1}Z'M */
              B5.at(edge.position) = arma::inv(Iz + Lz * zz.at(edge.descendant));
              zz.at(edge.descendant) *= B5.at(edge.position);
              D5.at(edge.position) = Iz - zz.at(edge.descendant) * Lz;
              xz.at(edge.descendant) *= D5.at(edge.position).t();
              zxz.at(edge.descendant) = D5.at(edge.position) * zxz.at(edge.descendant) * D5.at(edge.position).t();
            }

						inline void node           (const Tree::Edge& edge, const mat& d)
						{
							node (edge);
						}

						inline void increment      (const Tree::Edge& edge)
						{
              xz.at(edge.ancestor) += xz.at(edge.descendant); 
              zxz.at(edge.ancestor) += zxz.at(edge.descendant);
              xx.at(edge.ancestor) += xx.at(edge.descendant);
              zz.at(edge.ancestor) += zz.at(edge.descendant);
						}

						inline void increment      (const Tree::Edge& edge, const mat& d)
						{
							increment (edge);
						}

            inline void clean          (void)
            {
              clean_vector<mat> (xx);
              clean_vector<mat> (zz);
              clean_vector<mat> (xz);
              clean_vector<mat> (zxz);
            }

						inline mat operator()      (void)
						{
							do_recursion <Tree::Edge, Backbone> (*this, K->tree->begin());
							return xx.at(K->tree->root);
						}

						inline mat operator()      (const mat& d)
						{
							do_recursion <Tree::Edge, Backbone, const mat&> (*this, K->tree->begin(), d);
							return xx.at(K->tree->root);
						}
        };

        struct Linear
        {
          /* Ky */
          private:
            Precision<Conditional<Model>>* const   K;
            uword                                  dim_x,
                                                   dim_z;
            mat                                    Lx,
                                                   Lz,
                                                   S,
                                                   A,
                                                   B,
                                                   Ix,
                                                   Iz;
            rowvec                                 C,
                                                   D;

          public:
            cube                                   x,       // KX
                                                   z;       // KDMZ
            std::vector<rowvec>                    yx,      // y'KX
                                                   yz;      // y'KDMZ

            Linear                                 (Precision<Conditional<Model>>& p) :
              K (&p),
              dim_x (p.model->dim),
              dim_z (p.model->parent->dim),
              x (dim_x, dim_x, p.tree->ntips),
              z (dim_z, dim_x, p.tree->ntips),
              yx (p.tree->nnodes, arma::zeros<rowvec>(dim_x)),
              yz (p.tree->nnodes, arma::zeros<rowvec>(dim_z)),
              Lx (arma::zeros<mat>(dim_x, dim_x)),
              Lz (arma::zeros<mat>(dim_z, dim_z)),
              S (arma::zeros<mat>(dim_x, dim_z)),
              A (arma::zeros<mat>(dim_x, dim_x)),
              B (arma::zeros<mat>(dim_x, dim_x)),
              C (arma::zeros<rowvec>(dim_z)),
              D (arma::zeros<rowvec>(dim_z)),
              Ix (arma::eye<mat>(dim_x, dim_x)),
              Iz (arma::eye<mat>(dim_z, dim_z))
            {}

						inline void tip                         (const Tree::Edge& edge, mat& y)
						{
              x.slice(edge.descendant) = K->backbone.xx.at(edge.descendant); // ths shoud work right...?
              z.slice(edge.descendant) = K->backbone.xz.at(edge.descendant).t(); //added t() ...
              yx.at(edge.descendant) = arma::trans(x.slice(edge.descendant)*y.col(edge.descendant));
              yz.at(edge.descendant) = arma::trans(z.slice(edge.descendant)*y.col(edge.descendant));
              y.col(edge.descendant) = arma::trans(yx.at(edge.descendant)); // the row/column thing needs to be straightened out
						}

						inline void node                        (const Tree::Edge& edge, mat& y)
						{
              Lx = K->model->edgelength_matrix(edge);
              Lz = K->model->parent->edgelength_matrix(edge);
              S = K->model->covariance_matrix(edge);
              
              /* step 1; update K - KX(I + LX'KX)^{-1}X'K */
              A = K->backbone.A1.at(edge.position) * Lx;
              for (uword i : edge.descent())
              {
                B = A * x.slice(i); 
                z.slice(i) -= K->backbone.xz1.at(edge.position).t() * B;
                y.col(i) -= arma::trans(yx.at(edge.descendant) * B); // wrong order, switch someday
                x.slice(i) -= K->backbone.xx1.at(edge.position) * B;
              }
              
              yx.at(edge.descendant) *= K->backbone.A1.at(edge.position); 
              yz.at(edge.descendant) -= 
                yx.at(edge.descendant) * Lx * K->backbone.xz1.at(edge.position);

              /* step 2; update K - K(-XS + DMZL')(I + LZ'MZ - S'X'KDMZ + LZ'MD'KDMZ)^{-T} Z'MD'K   */
              C = -yx.at(edge.descendant) * S + yz.at(edge.descendant) * Lz;
              D = C * K->backbone.D2.at(edge.position).t();
              for (uword i : edge.descent())
              {
                x.slice(i) -= K->backbone.CD2.at(edge.position) * z.slice(i);
                y.col(i) -= arma::trans(D * z.slice(i)); // backwards, switch someday
                z.slice(i) = (Iz - K->backbone.BD2.at(edge.position)) * z.slice(i);
              }

              yx.at(edge.descendant) -= D * K->backbone.xz2.at(edge.position).t();
              yz.at(edge.descendant) -= D * K->backbone.zxz2.at(edge.position);

              /* step 3; update K - K(-DMZ - XSZ'MZ)(I + LZ'MZ + S'X'KDMZ + S'X'KXSZ'MZ)^{-1} X'K  */
              C = -yz.at(edge.descendant) - 
                yx.at(edge.descendant) * S * K->backbone.zz3.at(edge.position);
              D = C * K->backbone.D3.at(edge.position) * S.t();
              for (uword i : edge.descent())
              {
                y.col(i) -= arma::trans(D * x.slice(i));
                z.slice(i) -= K->backbone.BD3.at(edge.position) * x.slice(i);
                x.slice(i) = (Ix - K->backbone.CD3.at(edge.position)) * x.slice(i);
              }

              yx.at(edge.descendant) -= D * K->backbone.xx3.at(edge.position);
              yz.at(edge.descendant) -= D * K->backbone.xz3.at(edge.position);

              /* step 4 and 5: update 
                    D + XSZ' 
                    M - MZ(I + LZ'MZ)^{-1}Z'M */
              for (uword i : edge.descent())
              {
                z.slice(i) += K->backbone.zz4.at(edge.position) * S.t() * x.slice(i);
                z.slice(i) = K->backbone.B5.at(edge.position).t() * z.slice(i); 
              }

              yz.at(edge.descendant) += yx.at(edge.descendant) * S * K->backbone.zz4.at(edge.position);
              yz.at(edge.descendant) *= K->backbone.D5.at(edge.position).t();
						}

						inline void increment                   (const Tree::Edge& edge, mat& y)
						{
              yx.at(edge.ancestor) += yx.at(edge.descendant);
              yz.at(edge.ancestor) += yz.at(edge.descendant);
						}

            inline void clean          (void)
            {
              clean_vector<rowvec> (yx);
              clean_vector<rowvec> (yz);
            }

						inline mat operator()                   (mat y)
						{
							do_recursion<Tree::Edge, Linear, mat&>(*this, K->tree->begin(), y);
							return y;
						}
        };

        struct Quadratic
        {
          /* y'Ky */
          private:
            Precision<Conditional<Model>>* const    K;
            uword                                   dim_x,
                                                    dim_z;
            rowvec                                  C,
                                                    D;
            mat                                     S,
                                                    L;

          public:
            std::vector<rowvec>           yx,     // y'KX
                                          yz;     // y'KDMZ
            std::vector<double>           yy;     // y'Ky

            Quadratic                  (Precision<Conditional<Model>>& p) :
              K (&p),
              dim_x (p.model->dim),
              dim_z (p.model->parent->dim),
              yx (p.tree->nnodes, arma::zeros<rowvec>(dim_x)),     
              yz (p.tree->nnodes, arma::zeros<rowvec>(dim_z)),
              yy (p.tree->nnodes, 0.0),
              C (arma::zeros<rowvec>(dim_z)),
              D (arma::zeros<rowvec>(dim_z)),
              S (arma::zeros<mat>(dim_x, dim_z)),
              L (arma::zeros<mat>(dim_z, dim_z))
            {}

						inline void tip            (const Tree::Edge& edge, const mat& y)
						{
              yx.at(edge.descendant) = 
                arma::trans(K->backbone.xx.at(edge.descendant) * y.col(edge.descendant));
              yy.at(edge.descendant) = 
                arma::dot(y.col(edge.descendant), yx.at(edge.descendant));
              yz.at(edge.descendant) = yx.at(edge.descendant) * 
                K->model->covariance_matrix(edge) * K->backbone.zz.at(edge.descendant);
						}

						inline void node           (const Tree::Edge& edge, const mat& y)
						{
              S = K->model->covariance_matrix(edge);
              L = K->model->parent->edgelength_matrix(edge);

              /* step 1 */
              yy.at(edge.descendant) -= arma::dot(yx.at(edge.descendant) * 
                K->backbone.A1.at(edge.position) * 
                K->model->edgelength_matrix(edge), yx.at(edge.descendant));
              yx.at(edge.descendant) *= K->backbone.A1.at(edge.position); // intentionally update first
              yz.at(edge.descendant) -= 
                yx.at(edge.descendant) * K->model->edgelength_matrix(edge) *
                K->backbone.xz1.at(edge.position);

              /* step 2 */
              C = -yx.at(edge.descendant) * S + yz.at(edge.descendant) * L;
              D = C * K->backbone.D2.at(edge.position).t();
              yy.at(edge.descendant) -= arma::dot(D, yz.at(edge.descendant));
              vec xy = yx.at(edge.descendant).t() - // necessary as the updates to K are asymmetric
                K->backbone.CD2.at(edge.position) * yz.at(edge.descendant).t();
              yx.at(edge.descendant) -= D * K->backbone.xz2.at(edge.position).t(); 
              yz.at(edge.descendant) -= D * K->backbone.zxz2.at(edge.position);

              /* step 3 */
              C = -yz.at(edge.descendant) - yx.at(edge.descendant) * S * K->backbone.zz3.at(edge.position);
              D = C * K->backbone.D3.at(edge.position) * S.t();
              yy.at(edge.descendant) -= arma::dot(D, xy);
              yx.at(edge.descendant) -= D * K->backbone.xx3.at(edge.position);
              yz.at(edge.descendant) -= D * K->backbone.xz3.at(edge.position);

              /* step 4 and 5 */
              yz.at(edge.descendant) += yx.at(edge.descendant) * S * K->backbone.zz4.at(edge.position);
              yz.at(edge.descendant) *= K->backbone.D5.at(edge.position).t();
						}

						inline void increment      (const Tree::Edge& edge, const mat& y)
						{
              yx.at(edge.ancestor) += yx.at(edge.descendant);
              yz.at(edge.ancestor) += yz.at(edge.descendant);
              yy.at(edge.ancestor) += yy.at(edge.descendant);
						}

            inline void clean          (void)
            {
              clean_vector<rowvec> (yx);
              clean_vector<rowvec> (yz);
              clean_vector<double> (yy);
            }

						inline double operator()   (const mat& y)
						{
							do_recursion<Tree::Edge, Quadratic, const mat&>(*this, K->tree->begin(), y);
							return yy.at(K->tree->root);
						}
        };

        struct Diagonal
        {
          /* diag(K) */
          private:
            Precision<Conditional<Model>>* const    K;
            uword                                   dim_x,
                                                    dim_z;
            mat                                     yx, // y'KX
                                                    yz, // y'KDMZ
                                                    yy; // y'Ky
            mat                                     C,  // temporaries
                                                    D,
                                                    S,
                                                    L;

          public:
            mat                       diag;

            Diagonal                  (Precision<Conditional<Model>>& p) :
              K (&p),
              dim_x (p.model->dim),
              dim_z (p.model->parent->dim),
              diag (p.model->dim, p.tree->ntips, arma::fill::zeros),
              yx (arma::zeros<mat>(dim_x, dim_x)),
              yz (arma::zeros<mat>(dim_x, dim_z)),
              yy (arma::zeros<mat>(dim_x, dim_x)),
              C (arma::zeros<mat>(dim_x, dim_z)),
              D (arma::zeros<mat>(dim_x, dim_z)),
              S (arma::zeros<mat>(dim_x, dim_z)),
              L (arma::zeros<mat>(dim_z, dim_z))
            {}

						inline void tip           (const Tree::Rooter& edge)
						{
              yx = yy = K->backbone.xx.at(edge.descendant);
              yz = K->backbone.xz.at(edge.descendant);
						}

						inline void node          (const Tree::Rooter& edge)
						{
              S = K->model->covariance_matrix(edge);
              L = K->model->parent->edgelength_matrix(edge);

              /* step 1 */
              yy -= yx * K->backbone.A1.at(edge.position) * 
                K->model->edgelength_matrix(edge) * yx.t();
              yx *= K->backbone.A1.at(edge.position);
              yz -= yx * K->model->edgelength_matrix(edge) *
                K->backbone.xz1.at(edge.position);

              /* step 2 */
              C = -yx*S + yz*L;
              D = C * K->backbone.D2.at(edge.position).t();
              yy -= D * yz.t(); 
              mat xy = yx.t() - K->backbone.CD2.at(edge.position) * yz.t();
              yx -= D * K->backbone.xz2.at(edge.position).t(); 
              yz -= D * K->backbone.zxz2.at(edge.position);

              /* step 3 */
              C = -yz - yx * S * K->backbone.zz3.at(edge.position);
              D = C * K->backbone.D3.at(edge.position) * S.t();
              yy -= D * xy;
              yx -= D * K->backbone.xx3.at(edge.position);
              yz -= D * K->backbone.xz3.at(edge.position);

              /* step 4 and 5 */
              yz += yx * S * K->backbone.zz4.at(edge.position);
              yz *= K->backbone.D5.at(edge.position).t();
						}

						inline void increment     (const Tree::Rooter& edge)
						{
							if (edge.ancestor == K->tree->root)
								diag.col(edge.tip) = yy.diag();
						}

            inline void clean          (void)
            {
            }

						inline mat operator()     (void)
						{
							for(uword i=0; i < K->tree->ntips; ++i)
								do_recursion<Tree::Rooter, Diagonal>(*this, K->tree->begin(i));
							return diag;
						}
        };

        struct Determinant
        {
          /* ln |K| */
          private:
            Precision<Conditional<Model>>* const     K;

          public:
            std::vector<double>         logdet; 

            Determinant                 (Precision<Conditional<Model>>& p) :
              K (&p),
              logdet (p.tree->nnodes, 0.0)
            {}

						inline void tip             (const Tree::Edge& edge)
						{
              logdet.at(edge.descendant) +=
                log(arma::det(K->backbone.xx.at(edge.descendant)));
						}

						inline void node            (const Tree::Edge& edge)
						{ 
              logdet.at(edge.descendant) +=
                log(arma::det(K->backbone.A1.at(edge.position))) +
                log(arma::det(K->backbone.D2.at(edge.position))) +
                log(arma::det(K->backbone.D3.at(edge.position))) -
                2.0 * log(arma::det(K->backbone.B5.at(edge.position)));
						}

						inline void increment       (const Tree::Edge& edge)
						{
							logdet.at(edge.ancestor) += logdet.at(edge.descendant);
						}

            inline void clean          (void)
            {
              clean_vector<double> (logdet);
            }

						inline double operator()    (void)
						{
							do_recursion<Tree::Edge, Determinant>(*this, K->tree->begin());
							return logdet.at(K->tree->root);
						}
				};

        struct Offdiagonal
        {
          /* if conditional matrix is K_i - K_ij K_j^{-1} K_{ji}
             then this calculates K_ij z for arbitrary vector z */
          private:
            Precision<Conditional<Model>>* const   K;
            vec                                    S;

          public:
            std::vector<vec>                       xy;
            mat                                    o;

            Offdiagonal                    (Precision<Conditional<Model>>& p) :
              K (&p),
              S (p.model->dim, arma::fill::zeros),
              xy (p.tree->nnodes, arma::zeros<vec>(p.model->parent->dim)),
              o (p.model->dim, p.tree->ntips, arma::fill::zeros)
            {}

            inline void tip (const Tree::Edge& edge, const mat& y)
            {
              xy.at(edge.descendant) = y.col(edge.descendant);
              o.col(edge.descendant) =
                K->model->covariance_matrix(edge) *
                xy.at(edge.descendant);
            }

            inline void node (const Tree::Edge& edge, const mat& y)
            {
              S = K->model->covariance_matrix(edge) *
                xy.at(edge.descendant);
              for(auto i : edge.descent())
                o.col(i) += S;
            }

            inline void increment (const Tree::Edge& edge, const mat& y) 
            {
              xy.at(edge.ancestor) += xy.at(edge.descendant);
            }

            inline void clean (void)
            {
              clean_vector<vec>(xy);
            }

            inline mat operator() (const mat& y)
            {
							do_recursion<Tree::Edge, Offdiagonal, const mat&>(*this, K->tree->begin(), y);
							return o;
            }
        };

        struct InverseDiagonal
        {
          private:
            Precision<Conditional<Model>>* const   K;
            mat                                    sigma;

          public:
            InverseDiagonal                    (Precision<Conditional<Model>>& p) :
              K (&p),
              sigma (p.model->dim, p.tree->ntips, arma::fill::zeros)
            {}

            inline mat operator() (void)
            {
              for(uword i=0; i<K->tree->ntips; ++i)
                sigma.col(i) = K->model->tip_variance(K->tree->begin(i)); //NOT complete TODO add DMD
              return sigma;
            }
        };

        Backbone            backbone;
        Linear              linear;
        Quadratic           quadratic;
        Diagonal            diagonal;
        Determinant         determinant;
        Offdiagonal         offdiagonal;
        InverseDiagonal     inverse_diagonal;

				Precision (Conditional<Model>& m, Tree& t) : /* cnstr */
					model (&m),
					tree (&t),
					backbone (*this),
					linear (*this),
					quadratic (*this),
					diagonal (*this),
					determinant (*this),
          offdiagonal (*this),
          inverse_diagonal (*this)
  			{
					backbone ();
  			}

				Precision (Conditional<Model>& m, Tree& t, const mat& d) : /* cnstr */
					model (&m),
					tree (&t),
					backbone (*this),
					linear (*this),
					quadratic (*this),
					diagonal (*this),
					determinant (*this),
          offdiagonal (*this),
          inverse_diagonal (*this)
				{
					backbone (d);
				}
    };
} /* end namespace epee */

#endif
