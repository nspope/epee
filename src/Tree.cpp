#include <RcppArmadillo.h> 
#include "Tree.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

namespace epee 
{
  using mat = arma::mat;
  using umat = arma::umat;
  using uvec = arma::uvec;
  using vec = arma::vec;
  using uword = arma::uword;
  using subview_col = arma::subview_col<double>;

  /* definitions for base class Edge */
  typename Tree::Edge& 
    Tree::Edge::operator++ (void)
    {
      ++position;
      if (position < tree->nedges)
      {
        descendant = tree->edge_list(position,1);
        ancestor = tree->edge_list(position,0);
        is_tip = descendant < tree->ntips;
      }
      else
      {
        position = tree->nedges;
        descendant = tree->root;
        ancestor = arma::datum::nan;
        is_tip = false;
      }
      return *this;
    }

  typename Tree::Edge 
    Tree::Edge::operator++ (int)
    {
      Edge temp = *this;
      ++*this;
      return temp;
    }

  uword Tree::Edge::operator() (void)
    {
      return position;
    }

  bool Tree::Edge::operator== 
    (const Edge& other)
    {
      return position == other.position;
    }

  bool Tree::Edge::operator!= 
    (const Edge& other)
    {
      return position != other.position;
    }

  bool Tree::Edge::operator< 
    (const Edge& other)
    {
      return position < other.position;
    }

  bool Tree::Edge::operator> 
    (const Edge& other)
    {
      return position > other.position;
    }

  bool Tree::Edge::operator<=
    (const Edge& other)
    {
      return position <= other.position;
    }

  bool Tree::Edge::operator>=
    (const Edge& other)
    {
      return position >= other.position;
    }

  const double& 
    Tree::Edge::length (void) const
    {
      return (tree->lengths(position));
    }

  const double& 
    Tree::Edge::basal_height (void) const
    {
      return (tree->heights(ancestor));
    }

  const double& 
    Tree::Edge::distal_height (void) const
    {
      return (tree->heights(descendant));
    }

  const subview_col
    Tree::Edge::traits (void) const
    {
      if (is_tip)
      {
        return (tree->Y.col(descendant));
      }
      else
      {
        return (tree->Ymiss.col(0)); /* TODO: at some point have tree store ancestral values */
      }
    }

  const uvec&
    Tree::Edge::descent (void) const
    {
      /* if subtending a tip, returns line of descent; if an internal edge, 
         returns all descendants from the distal end of the edge */
      if (is_tip)
        return tree->descent.at(descendant);
      else
        return tree->descendants.at(descendant);
    }

  /* definitions for derived class Rooter : Edge */

  typename Tree::Rooter& 
    Tree::Rooter::operator++ (void)
    {
      ++depth;
      if (depth < tree->descent[tip].n_elem)
      {
        position = tree->descent[tip](depth);
        descendant = tree->edge_list(position,1);
        ancestor = tree->edge_list(position,0);
        is_tip = descendant < tree->ntips;
      }
      else
      {
        depth = tree->descent[tip].n_elem;
        position = tree->nedges;
        descendant = tree->root;
        ancestor = arma::datum::nan;
        is_tip = false;
      }
      return *this;
    }

  typename Tree::Rooter 
    Tree::Rooter::operator++ (int)
    {
      Rooter temp = *this;
      ++*this;
      return temp;
    }

  typename Tree::Edge 
    Tree::Edge::end (void)
    {
      return tree->end(); 
    };

} /* namespace epee */
