#ifndef COVARIANCE_H
#define COVARIANCE_H

#include <RcppArmadillo.h> 
#include "inputs.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

namespace epee 
{

class Phylo
{
  /* relevant tests:
   * does this work with polytomies?
   * does this work with strictly bifurcating trees?
   * internally reorder so as to be postorder?
   */
  public:
    const umat        edge;
          vec         edge_length;
    const uword       edges,
                      nodes,
                      tips,
                      root;
    luvec             descendants;
    double            average_depth;

    Phylo (const umat& edge_list, const vec& lengths, const uword ntip) :
      edge        (edge_list),
      edge_length (lengths),
      edges       (edge.n_rows),
      nodes       (arma::max(edge.col(0))+1),
      tips        (ntip),
      root        (edge(edges-1,0)),
      descendants (nodes)
    {
      /* construct lines of descent for each tip */
      uword parent,
            tip,
            e;

      uvec edge_count(tips, arma::fill::zeros),
           parent_count(nodes, arma::fill::zeros);

      std::vector<uvec> descent (tips);

      uvec edge_order   = edge.col(1);
      edge_order.insert_rows(edges, 1);
      edge_order(edges) = root;
      edge_order        = arma::sort_index(edge_order);

      /* determine size of vectors */
      for(tip=0;tip<tips;++tip)
      {
        descent[tip]    = uvec(tips);
        parent          = tip;
        while(parent != root)
        {
          e             = edge_order(parent);
          parent        = edge (e, 0);
          descent[tip](edge_count(tip)) = e;
          ++edge_count(tip);
          ++parent_count(parent);
        }
        descent[tip]    = descent[tip].head(edge_count(tip));
      }

      /* allocate space */
      for(tip=0;tip<tips;++tip)
        descendants[tip]    = uvec({tip});
      for(parent=tips;parent<nodes;++parent)
        descendants[parent] = uvec(parent_count(parent));

      /* now use counters for tracking position */
      edge_count.zeros();
      parent_count.zeros();

      /* fill vectors */
      for(tip=0;tip<tips;++tip)
      {
        parent = tip;
        while(parent != root)
        {
          e       = descent[tip].at(edge_count(tip)); 
          parent  = edge(e, 0);
          descendants[parent](parent_count(parent)) = tip;
          ++edge_count(tip);
          ++parent_count(parent);
        }
      }

      /* this is theoretically the complexity */
      average_depth = double(arma::sum(edge_count))/double(tips);
    }
};

struct Storage
{
  const uvec  tx, // indices for conditional/marginal traits
              tz; // \/
  const uword mx, // dimensionality
              mz; // \/
  cube  A,        // global storage for recursed quantities
        Ay,       // ||
        B,        // ||
        By,       // ||
        C,        // ||
        D,        // ||
        Dy,       // \/
        U1,       // global storage for update matrices
        U2,       // ||
        U3,       // ||
        U4,       // ||
        U6,       // ||
        Ui1,      // ||
        Ui2,      // ||
        Ui3,      // ||
        Ui4,      // ||
        Ui6,      // ||
        Ud2,      // ||
        Ud3,      // ||
        Ud4,      // ||
        Ud6,      // \/
        Udi2,
        Udi3,
        Udi6,
        G,        // ||
        Gt,
        F,
        H,
        dF,       // global storage for derivatives
        dH,
        dG,
        dA,
        dAy,
        dB,
        dBy,
        dC,
        dD,
        dDy;        // \/
  mat   Lx,  // global storage for covariance parameters
        Lxz, // ||
        Lz,  // ||
        Sx,   // ||
        Sxz,  // ||
        Sz,   // ||
        Eye,      // \/
        Ix,
        Iz,
        yA,
        dyA,
        Mx,
        Mz,
        yx,       // global storage for conditional traits,
        yz,       // for marginal traits,
        d,        // for diagonal modification
        dLx,      // storage for derivatives
        dLxz,
        dLzx,
        dLz,
        dSx,
        dSxz,
        dSzx,
        dSz,
        dY;
  umat  mix,
        miz;

  //Storage (const uvec& type, const uword nodes, const uword tips) :
  Storage (const uvec& type, const mat& y, const uword nodes) :
    tz   (arma::find(type==0)),
    tx   (arma::find(type>0)),
    mx   (tx.n_elem),
    mz   (tz.n_elem),
    A    (mx, mx, nodes),
    Ay   (1, mx, nodes),
    B    (mx, mz, nodes),
    By   (1, mz, nodes),
    C    (mz, mz, nodes),
    D    (mz, mz, nodes),
    Dy   (1, mz, nodes),
    U1   (mx, mx, nodes),
    U2   (mx, mx, nodes),
    U3   (mx, mx, nodes),
    U4   (mx, mx, nodes),
    U6   (mz, mz, nodes),
    Ui1  (mx, mx, nodes),
    Ui2  (mx, mx, nodes),
    Ui3  (mx, mx, nodes),
    Ui4  (mx, mx, nodes),
    Ui6  (mz, mz, nodes),
    Ud2  (mz, mz, nodes),
    Ud3  (mz, mz, nodes),
    Ud4  (mx, mx, nodes),
    Ud6  (mz, mz, nodes),
    Udi2 (mz, mz, nodes),
    Udi3 (mz, mz, nodes),
    Udi6 (mz, mz, nodes),
    G    (mx, mz, y.n_cols),
    Gt   (mz, mx, y.n_cols),
    F    (mx, mx, y.n_cols),
    H    (mx, mz, y.n_cols),
    dF   (mx, mx, y.n_cols),
    dH   (mx, mz, y.n_cols),
    dG   (mx, mz, y.n_cols),
    dA   (mx, mx, nodes),
    dAy  (1, mx, nodes),
    dB   (mx, mz, nodes),
    dBy  (1, mz, nodes),
    dC   (mz, mz, nodes),
    dD   (mz, mz, nodes),
    dDy  (1, mz, nodes),
    Lx   (arma::eye(mx,mx)),
    Lxz  (arma::zeros(mx,mz)),
    Lz   (arma::eye(mz,mz)),
    Sx   (arma::eye(mx,mx)),
    Sxz  (arma::zeros(mx,mz)),
    Sz   (arma::eye(mz,mz)),
    Eye  (arma::eye(mx+mz,mx+mz)),
    Ix   (arma::eye(mx,mx)),
    Iz   (arma::eye(mz,mz)),
    yA   (arma::zeros(mx,1)),
    dyA  (arma::zeros(mx,1)),
    Mx   (arma::zeros(mx, y.n_cols)),
    Mz   (arma::zeros(mz, y.n_cols)),
    yx   (y.rows(tx)),
    yz   (y.rows(tz)),
    d    (mz, y.n_cols, arma::fill::zeros),
    dLx  (mx, mx, arma::fill::zeros),
    dLxz (mx, mz, arma::fill::zeros),
    dLzx (mz, mx, arma::fill::zeros),
    dLz  (mz, mz, arma::fill::zeros),
    dSx  (mx, mx, arma::fill::zeros),
    dSxz (mx, mz, arma::fill::zeros),
    dSzx (mz, mx, arma::fill::zeros),
    dSz  (mz, mz, arma::fill::zeros),
    dY   (mx, y.n_cols),
    mix  (mx, nodes, arma::fill::zeros),
    miz  (mz, nodes, arma::fill::zeros)
  {
  }

  void clear                      (void);
  void clear_derivatives          (void);
  void clear_derivatives_marginal (void);
  void clear_marginal             (void);
  void update                     (const mat&, const mat&, const mat&);
};

class Covariance
{
  friend class Epee;

  private:
    uvec    ze  = {0}, // useful for subsetting matrices
            tip = {0}; // \/
    Phylo   tree;      // phylogeny,
    Storage storage;   // global storage for recursed quantities,
    uword   ce,        // index of child node,
            pe,        // index of parent node,
            e;         // index of edge,
    mat     Lx,        // local copies of covariance matrices 
            Lz,        // ||
            Lxz;       // \/
    uvec    dx,        // conditional traits to update,
            dz;        // marginal traits to update,
    double  le,        // edge length,
            ld,        // temporaries needed for determinant
            ls,        // \/
            dAyy = 0., // derivative multiplier
            dDyy = 0., // ||
            dld  = 0.; // \/
    bool    is_a_tip,  // am I a tip?
            upx,       // update conditional traits?
            upz;       // update marginal traits?

  public:
    double  Ayy      = 0.,
            logdet   = 0.,
            Dyy      = 0.,
            logdet_marginal  = 0.;
    mat     Y,        // linear products and diagonal
            Z,        // linear product, make offdiagonal
            diag,
            dy,
            dr,
            dL,
            dS;     // ||

  Covariance (const umat& edge, 
              const vec&  edge_length, 
              const mat&  y,
              const uvec& type) :
    tree (edge, edge_length, y.n_cols),
    storage (type, y, tree.nodes), 
    Y (storage.mx, tree.tips, arma::fill::zeros),
    Z (storage.mx, tree.tips, arma::fill::zeros),
    diag (storage.mx, tree.tips, arma::fill::zeros),
    dy (storage.mx, tree.tips, arma::fill::zeros),
    dr (storage.mz, tree.tips, arma::fill::zeros),
    dL (storage.mx+storage.mz, storage.mx+storage.mz, arma::fill::zeros),
    dS (storage.mx+storage.mz, storage.mx+storage.mz, arma::fill::zeros)
  {
    missing_data_pattern ();
  }

  /* alternative constructor with design matrix */
  Covariance (const Design& design) :
    tree (design.edge, design.edge_length, design.tips),
    storage (design.type, design.traits, tree.nodes), 
    Y (storage.mx, tree.tips, arma::fill::zeros),
    Z (storage.mx, tree.tips, arma::fill::zeros),
    diag (storage.mx, tree.tips, arma::fill::zeros),
    dy (storage.mx, tree.tips, arma::fill::zeros),
    dr (storage.mz, tree.tips, arma::fill::zeros),
    dL (storage.mx+storage.mz, storage.mx+storage.mz, arma::fill::zeros),
    dS (storage.mx+storage.mz, storage.mx+storage.mz, arma::fill::zeros)
  {
    missing_data_pattern ();
  }

  void missing_data_pattern (void)
  {
    // find missing data pattern
    for(e=0; e<tree.edges; ++e)
    {
      ce = tree.edge(e,1);
      pe = tree.edge(e,0);
      if (ce < tree.tips)
      {
        tip(0) = ce;
        storage.mix(arma::find_finite(storage.yx.cols(tip)), tip) += 1;
        storage.miz(arma::find_finite(storage.yz.cols(tip)), tip) += 1;
      }
      storage.mix.col(pe) += storage.mix.col(ce);
      storage.miz.col(pe) += storage.miz.col(ce);
    }
  }

  /* the marginal part of the covariance matrix */
  void traverse_tree_marginal (void);
  void start_marginal         (void);
  void edge_up_marginal       (void);
  void update_tip_marginal    (void);
  void update_node_marginal   (void);
  void increment_marginal     (void);
  void missing_data           (void);

  /* the conditional part of the covariance matrix */
  void traverse_tree (const mat&, const mat&);
  void start         (const mat&, const mat&);
  void edge_up       (void);
  void update_tip    (void);
  void update_node   (void);
  void increment     (void);
  //TODO: some of the below should be const
  void equationA1    (mat&, mat&, mat&, mat&, mat&, mat&, mat&);
  void equationA2    (mat&, mat&, mat&, mat&, mat&, mat&, mat&, mat&, mat&, mat&);
  void equationA3    (mat&, mat&, mat&, mat&, mat&, mat&, mat&, mat&, mat&, mat&);
  void equationA4    (mat&, mat&, mat&, mat&, mat&, mat&, mat&, mat&, mat&, mat&);
  void equationA5    (mat&, mat&, mat&, mat&, mat&, mat&, mat&);
  void equationA6    (mat&, mat&, mat&, mat&, mat&, mat&, mat&, mat&, mat&);

  /* derivatives for marginal part of covariance matrix */
  void traverse_reverse_marginal      (const mat&, double, double);
  void start_reverse_marginal         (const mat&, const double, const double);
  void edge_down_marginal             (void);
  void derivatives_tip_marginal       (void);
  void derivatives_node_marginal      (void);
  void derivatives_increment_marginal (void);

  /* derivatives for conditional part of covariance matrix */
  void traverse_reverse       (const mat&, double, double);
  void start_reverse          (const mat&, const double, const double);
  void edge_down              (void);
  void derivatives_tip        (void);
  void derivatives_node       (void);
  void derivatives_increment  (void);
  void derivatives_equationA6 (mat&, mat&, mat&, mat&, mat&, mat&, mat&, mat&, mat&, mat&,
                               mat&, mat&, mat&, mat&, mat&, mat&, mat&);
  void derivatives_equationA5 (mat&, mat&, mat&, mat&, mat&, mat&, mat&, mat&, mat&, mat&,
                               mat&, mat&, mat&, mat&);
  void derivatives_equationA4 (mat&, mat&, mat&, mat&, mat&, mat&, mat&, mat&, mat&, mat&,
                               mat&, mat&, mat&, mat&, mat&, mat&, mat&, mat&);
  void derivatives_equationA3 (mat&, mat&, mat&, mat&, mat&, mat&, mat&, mat&, mat&, mat&,
                               mat&, mat&, mat&, mat&, mat&, mat&);
  void derivatives_equationA2 (mat&, mat&, mat&, mat&, mat&, mat&, mat&, mat&, mat&, mat&,
                               mat&, mat&, mat&, mat&, mat&, mat&);
  void derivatives_equationA1 (mat&, mat&, mat&, mat&, mat&, mat&, mat&, mat&, mat&, mat&,
                               mat&, mat&, mat&);

  /* misc */
  void parameters (const mat&, const mat&, const mat&);
};


} /* end namespace epee */
#endif
