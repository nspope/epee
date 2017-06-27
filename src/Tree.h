#ifndef TREE_H
#define TREE_H

#include <RcppArmadillo.h> 

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

namespace epee 
{

  using rowvec = arma::rowvec;
  using vec = arma::vec;
  using mat = arma::mat;
  using uvec = arma::uvec;
  using uword = arma::uword;
  using umat = arma::umat;
  using cube = arma::cube;
  using subview_col = arma::subview_col<double>;
  using map_uv = std::map<uword, uvec>;

  class Tree
  {
    /* contains tree topology, trait data, and an iterator-ish
       class that can be used to access specific edges/nodes.
       Also contains lines of descent for each tip, and set of
       descendants for each internal edge. */
    private:  
      const mat			Ymiss;
      const umat		edge_list; 
      const vec			lengths;
      const vec			heights;
      map_uv        descent;
      map_uv        descendants;

    public:
      const uword		ntraits;
      const uword		ntips;
      const uword		nnodes;
      const uword		nedges;
      const uword   root;
      const mat 		Y;
      double        average_depth;

      Tree (const mat Y_, const umat E_, const vec L_, const vec H_) :
        Y(Y_), 
        edge_list(E_), 
        lengths(L_), 
        heights(H_),
        ntraits(Y_.n_rows), 
        ntips(Y_.n_cols), 
        nnodes(H_.n_elem), 
        nedges(L_.n_elem),
        Ymiss(Y_.n_rows, 1, arma::fill::zeros),
        root(edge_list(E_.n_rows-1, 0)) 
      {
        /* dimension check */
        if ( (ntips-1)*2 != nedges)
          Rcpp::stop("Dimension mismatch: number of tips and number of edges");

        /* construct lines of descent for each tip */
        uword parent,
              tip,
              edge;

        uvec edge_count(ntips, arma::fill::zeros),
             parent_count(nnodes, arma::fill::zeros);

        uvec edge_order = edge_list.col(1);
        edge_order.insert_rows(nedges, 1);
        edge_order(nedges) = root;
        edge_order = sort_index(edge_order);

        /* determine size of vectors */
        for(tip=0;tip<ntips;++tip)
        {
          descent[tip] = uvec(ntips);
          parent = tip;
          while(parent != root)
          {
            edge = edge_order(parent);
            parent = edge_list(edge, 0);
            descent[tip](edge_count(tip)) = edge;
            ++edge_count(tip);
            ++parent_count(parent);
          }
          descent[tip] = descent[tip].head(edge_count(tip));
        }

        /* allocate space */
        for(parent=ntips;parent<nnodes;++parent)
          descendants[parent] = uvec(parent_count(parent));

        /* now use counters for tracking position */
        edge_count.zeros();
        parent_count.zeros();

        /* fill vectors */
        for(tip=0;tip<ntips;++tip)
        {
          parent = tip;
          while(parent != root)
          {
            edge = descent[tip].at(edge_count(tip)); 
            parent = edge_list(edge, 0);
            descendants[parent](parent_count(parent)) = tip;
            ++edge_count(tip);
            ++parent_count(parent);
          }
        }

        /* this is theoretically the complexity */
        average_depth = double(sum(edge_count))/double(ntips);
      }

      vec tip_heights (void)
      {
        return heights(arma::span(0, ntips-1));
      }

      inline double tip_height (const uword i) // single height
      {
        return heights(i);
      }

      /* an iterator-ish nested class for edges */
      class Edge
      {
        protected:
          Tree* const		tree;

        public:
          uword 				position;
          uword         descendant;
          uword         ancestor;
          bool					is_tip;

          Edge (Tree& tree_, uword position_) : 
            tree(&tree_), 
            position(position_), 
            descendant(tree_.edge_list(position_,1)),
            ancestor(tree_.edge_list(position_,0)),
            is_tip(descendant < tree_.ntips)
          {}

          /* operators */
          Edge& 				operator++ (void);
          Edge 					operator++ (int);
          uword					operator() (void);
          bool					operator== (const Edge&);
          bool					operator!= (const Edge&);
          bool					operator> (const Edge&);
          bool					operator< (const Edge&);
          bool					operator>= (const Edge&);
          bool					operator<= (const Edge&);

          /* member functions that return references to tree
             class */
          const double&									length (void) const;
          const double& 								basal_height (void) const;
          const double& 								distal_height (void) const;
          const subview_col   	        traits (void) const;
          const uvec&                   descent (void) const;
          Edge                   				end (void);
      };

      /* container boundaries */
      Edge begin (void)
			{
        return Edge (*this, 0);	
      };

      Edge end (void) 
			{
        return Edge (*this, nedges-1);
      };
			
      /* an iterator-ish nested class for edges that moves vertically on the tree */
      class Rooter : public Edge 
      {

        public:
          const uword   tip;
          uword         depth;

          Rooter (Tree& tree_, uword tip_) : 
            Edge (tree_, tree_.descent[tip_](0)),
            tip(tip_),
            depth(0)
          {}

          /* operators */
          Rooter&				operator++ (void);
          Rooter				operator++ (int);
      };

      /* container boundaries */
      Rooter begin (uword tip) {
        return Rooter (*this, tip);	
      };
  };

} /* end namespace phylopar */

#endif
