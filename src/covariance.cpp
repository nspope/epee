#include "covariance.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

/* new idea:
 *
 * Have covariance store both traits and mean.
 * Mean is a parameter that can be updated.
 * Traits are fixed values.
 * In the *marginal* updates, the mean is subtracted from the trait (errors are calculated)
 * And also the mean for the categorical traits is updated (saved as Z)
 */


namespace epee 
{

void Storage::clear (void)
{
  A.zeros();
  Ay.zeros();
  B.zeros();
  By.zeros();
  C.zeros();
  F.zeros();
  G.zeros();
  U1.zeros();
  U2.zeros();
  U3.zeros();
  U4.zeros();
  U6.zeros();
  Ui1.zeros();
  Ui2.zeros();
  Ui3.zeros();
  Ui4.zeros();
  Ui6.zeros();
  Ud2.zeros();
  Ud3.zeros();
  Ud4.zeros();
  Udi2.zeros();
  Udi3.zeros();
}

void Storage::clear_derivatives (void)
{
  dLx.zeros();
  dLz.zeros();
  dLxz.zeros();
  dLzx.zeros();
  dSx.zeros();
  dSz.zeros();
  dSxz.zeros();
  dSzx.zeros();
  dF.zeros();
  dG.zeros();
  dD.zeros();
  dA.zeros();
  dAy.zeros();
  dB.zeros();
  dBy.zeros();
  dC.zeros();
}

void Storage::clear_derivatives_marginal (void)
{
  dLz.zeros();
  dSz.zeros();
  dLxz.zeros();
  dSxz.zeros();
  dH.zeros();
  dD.zeros();
  dDy.zeros();
}

void Storage::clear_marginal (void)
{
  D.zeros();
  Dy.zeros();
  H.zeros();
  Ud6.zeros();
  Udi6.zeros();
}

void Storage::update (const mat& L, const mat& S, const mat& M)
{
  /* check dimensions TODO */
  Lx   = L(tx,tx);
  Lxz  = L(tx,tz);
  Lz   = L(tz,tz);
  Sx   = S(tx,tx);
  Sxz  = S(tx,tz);
  Sz   = S(tz,tz);
  Mx   = M.rows(tx);
  Mz   = M.rows(tz);
  /* TODO: any need to clean house here? */
}

void Covariance::traverse_tree_marginal (void)
{
  start_marginal ();
  for (e = 0; e < tree.edges; ++e)
  {
    edge_up_marginal ();
    if (is_a_tip)
      update_tip_marginal  ();
    else
      update_node_marginal ();
    increment_marginal ();
  }
}

void Covariance::start_marginal (void)
{
  Dyy       = 0.;
  logdet_marginal  = 0.;
  storage.clear_marginal ();
  Z = storage.Mx;
}

void Covariance::edge_up_marginal (void)
{
  pe       = tree.edge(e,0);
  ce       = tree.edge(e,1);
  le       = tree.edge_length(e);
  is_a_tip = ce < tree.tips;
  missing_data (); 
  Lxz      = storage.Lxz*le;
  Lz       = storage.Lz*le;
}

void Covariance::update_tip_marginal (void)
{
  if (upz)
  {
    mat&    Dy    = storage.Dy.slice(ce);
    mat&    D     = storage.D.slice(ce);
    mat&    H     = storage.H.slice(ce);
    mat&    Sxz   = storage.Sxz;
    mat&    Sz    = storage.Sz;
    vec     yz    = storage.yz.col(ce), //could this be more efficient?
            mu    = storage.Mz.col(ce);

    yz.replace(arma::datum::nan, 0.);
    mu.replace(arma::datum::nan, 0.);

    Lz      += Sz;
    D(dz,dz) = arma::inv(Lz(dz,dz));
    Dy       = (yz - mu).t() * D;
    Dyy     += arma::dot(Dy, yz - mu);
    arma::log_det (ld, ls, D(dz,dz));
    logdet_marginal += ld; 
    if (upx)
    {
      Lxz       += Sxz;
      H          = Lxz * D; 
      Z.col(ce) += Lxz * Dy.t();
    }
  }
}

void Covariance::update_node_marginal (void)
{
  /*
   * (K_{xz} + l_e X_e \Lambda_{xz} Z_e') (K_{z} + l_e Z_e \Lambda_z Z_e') \implies
   *  K_{xz} K_{z}^{-1} - K_{xz} K_z^{-1} Z_e (...) \Lambda_z Z_e' K_z^{-1}
   *  then ...
   *  K_{xz} K_{z}^{-1} + l_e X_e \Lambda_{xz} Z_e' K_{z}^{-1}
   */
  if (upz)
  {
    mat& D = storage.D.slice(ce);
    mat& Dy = storage.Dy.slice(ce);
    mat& Udi6 = storage.Udi6.slice(ce);
    mat& Ud6 = storage.Ud6.slice(ce);
    mat& Iz = storage.Iz;

    Udi6 = Iz + Lz * D;
    Ud6  = arma::inv(Udi6);
    D    = D * Ud6;
    Dyy -= arma::dot(Dy * Ud6 * Lz, Dy);
    Dy   = Dy * Ud6;
    arma::log_det (ld, ls, Ud6);
    logdet_marginal += ld;

    if (upx)
      for (uword ti : tree.descendants.at(ce))
      {
        mat& H     = storage.H.slice(ti); 
        Z.col(ti) += -H * Lz * Dy.t() + Lxz * Dy.t();
        H          = H * Ud6 + Lxz * D; 
      }
  }
}

void Covariance::increment_marginal (void)
{
  if (upz)
  {
    storage.D.slice(pe)  += storage.D.slice(ce);
    storage.Dy.slice(pe) += storage.Dy.slice(ce);
  }
}

void Covariance::missing_data (void)
{
  /* figure out what traits to update */
  dx = arma::find(storage.mix.col(ce));
  dz = arma::find(storage.miz.col(ce));
  upx = dx.n_elem > 0;
  upz = dz.n_elem > 0;
  /* another problem: what if trait covariance is set to 0?                    */
  /* answer: this is fine as long as we don't do the generalized inverse thing */
}

void Covariance::traverse_tree (const mat& y, const mat& d)
{
  start (y, d);
  for (e = 0; e < tree.edges; ++e)
  {
    edge_up ();
    if (is_a_tip)
      update_tip ();
    else
      update_node ();
    increment ();
  }
}

void Covariance::start (const mat& y, const mat& d)
{
  logdet     = 0.;
  Ayy        = 0.;
  storage.yx = y;
  storage.d  = d;
  storage.clear();
  Y.fill(0.);
  diag.fill(0.);
  Y %= storage.yx; // for missing data pattern
  diag %= storage.yx;
}

void Covariance::edge_up (void)
{
  storage.yA.zeros();
  pe       = tree.edge(e,0);
  ce       = tree.edge(e,1);
  tip(0)   = ce;
  le       = tree.edge_length(e);
  is_a_tip = ce < tree.tips;

  missing_data (); 

  /* drop traits that are missing */
  Lx       = storage.Lx*le; 
  Lxz      = storage.Lxz*le;
  Lz       = storage.Lz*le;
}

void Covariance::update_tip (void)
{
  mat&   A     = storage.A.slice(ce);
  mat&   Ay    = storage.Ay.slice(ce);
  mat&   B     = storage.B.slice(ce);
  mat&   By    = storage.By.slice(ce);
  mat&   C     = storage.C.slice(ce);
  mat&   D     = storage.D.slice(ce);
  mat&   Sx    = storage.Sx;
  mat&   Sz    = storage.Sz;
  mat&   Sxz   = storage.Sxz;
  mat&   F     = storage.F.slice(ce);
  mat&   G     = storage.G.slice(ce);
  vec    yx    = storage.yx.col(ce), //could this be more efficient?
         d     = storage.d.col(ce);

  if (upx)
  {
    yx.replace(arma::datum::nan, 0.);
    d.replace(arma::datum::nan, 0.);
    Lx    += Sx + arma::diagmat(d);
    if (upz)
    {
      Lxz += Sxz;
      Lx  -= Lxz * D * Lxz.t();
      A(dx,dx) = arma::inv(Lx(dx,dx));
      B    = A * Lxz * D;
      C    = B.t() * Lxz * D;
      By   = yx.t() * B;
      G    = B;
    }
    else
    {
      A(dx,dx)    = arma::inv(Lx(dx,dx));
    }
    arma::log_det (ld, ls, A(dx,dx));
    logdet       += ld;
    Ay            = yx.t() * A;
    Ayy          += arma::dot(Ay, yx);
    diag.col(ce)  = arma::diagvec(A);
    Y.col(ce)    += Ay.t();
    F             = A;
  }
}

void Covariance::update_node (void)
{
  mat& A  = storage.A.slice(ce);
  mat& Ay = storage.Ay.slice(ce);
  mat& yA = storage.yA;
  mat& B  =storage.B.slice(ce);
  mat& By =storage.By.slice(ce);
  mat& C  =storage.C.slice(ce);
  mat& D  =storage.D.slice(ce);
  mat& U1 =storage.U1.slice(ce);
  mat& Ui1=storage.Ui1.slice(ce);
  mat& U2 =storage.U2.slice(ce);
  mat& Ui2=storage.Ui2.slice(ce);
  mat& Ud2 =storage.Ud2.slice(ce);
  mat& Udi2=storage.Udi2.slice(ce);
  mat& U3 =storage.U3.slice(ce);
  mat& Ui3=storage.Ui3.slice(ce);
  mat& Ud3 =storage.Ud3.slice(ce);
  mat& Udi3=storage.Udi3.slice(ce);
  mat& U4 =storage.U4.slice(ce);
  mat& Ui4=storage.Ui4.slice(ce);
  mat& Ud4 =storage.Ud4.slice(ce);
  mat& U6 =storage.U6.slice(ce);
  mat& Ui6=storage.Ui6.slice(ce);
  mat& Ud6 =storage.Ud6.slice(ce);
  mat& Udi6=storage.Udi6.slice(ce);

  if (upx)
    equationA1(A, Ay, B, By, C, U1, Ui1);
  if (upx && upz)
  {
    equationA2(A, Ay, yA, B, By, C, U2, Ui2, Ud2, Udi2);
    equationA3(A, Ay, yA, B, By, C, U3, Ui3, Ud3, Udi3);
    equationA4(A, Ay, B, By, C, D, U4, Ui4, Ud4, Udi6);
    equationA5(A, Ay, B, By, C, D, Udi6);
  }
  if (upz)
    equationA6(A, Ay, B, By, C, U6, Ui6, Ud6, Udi6);
}

void Covariance::increment ()
{
  if (upx)
  {
    storage.A.slice(pe)  += storage.A.slice(ce);
    storage.Ay.slice(pe) += storage.Ay.slice(ce);
  }
  if (upz && upx)
  {
    storage.B.slice(pe)  += storage.B.slice(ce);
    storage.By.slice(pe) += storage.By.slice(ce);
    storage.C.slice(pe)  += storage.C.slice(ce);
  }
}


void Covariance::equationA1 (mat& A, 
                             mat& Ay, 
                             mat& B, 
                             mat& By, 
                             mat& C, 
                             mat& U1, 
                             mat& Ui1)
{
  mat& Ix = storage.Ix;
  /* 
   * K_{x|z} + l_e X \Lambda_x X' \implies
   * K_{x|z}^{-1} - l_e K_{x|z}^{-1} X (I + l_e \Lambda_x X' K_{x|z}^{-1} X)^{-1} \Lambda_x X' K_{x|z}^{-1}
   */
  Ui1     = Ix + Lx * A;
  U1      = arma::inv(Ui1);
  arma::log_det (ld, ls, U1);
  logdet += ld;
  A       = A * U1;
  Ayy    -= arma::dot(Ay * U1 * Lx, Ay);
  Ay      = Ay * U1;
  if (upz)
  {
    By   -= Ay * Lx * B;
    C    -= B.t() * U1 * Lx * B;
    B     = U1.t() * B;
  }
  for (auto ti : tree.descendants.at(ce))
  {
    mat& F        = storage.F.slice(ti);
    mat& G        = storage.G.slice(ti);
    diag.col(ti) -= arma::diagvec(F.t() * U1 * Lx * F);
    Y.col(ti)    -= arma::trans(Ay * Lx * F);
    G            -= F.t() * Lx * B;
    F             = U1.t() * F;
  }
}

void Covariance::equationA2 (mat& A,
                             mat& Ay,
                             mat& yA,
                             mat& B,
                             mat& By,
                             mat& C,
                             mat& U2,
                             mat& Ui2,
                             mat& Ud2,
                             mat& Udi2)
{
  mat& Ix = storage.Ix;
  mat& Iz = storage.Iz;
  /* 
   * K_{x|z} - l_e X \Lambda_{xz} Z' K_{z}^{-1} K_{zx} \implies 
   *  K_{x|z}^{-1} + 
   *    l_e K_{x|z}^{-1} X (I - l_e \Lambda_{xz} Z' K_{z}^{-1} K_{zx} K_{x|z}^{-1} X)^{-1} \Lambda_{xz} Z' K_{z}^{-1} K_{zx} K_{x|z}^{-1}     
   */
  Ui2     = Ix - Lxz * B.t();
  Udi2    = Iz - B.t() * Lxz;
  U2      = arma::inv(Ui2);
  Ud2     = arma::inv(Udi2);
  arma::log_det (ld, ls, U2);
  logdet += ld;
  A       = A * U2;
  yA      = Ay.t() + A * Lxz * By.t();
  Ay      = Ay * U2;
  Ayy    += arma::dot(Ay * Lxz, By);
  By     += Ay * Lxz * C;
  B      += A * Lxz * C;
  C       = Ud2 * C;
  for (uword ti : tree.descendants.at(ce))
  {
    mat& F        = storage.F.slice(ti);
    mat& G        = storage.G.slice(ti);
    diag.col(ti) += arma::diagvec(F.t() * U2 * Lxz * G.t());
    Y.col(ti)    += arma::trans(Ay * Lxz * G.t());
    G            += F.t() * Lxz * C;
    F            += A * Lxz * (G - F.t() * Lxz * C).t();
  }
}

void Covariance::equationA3 (mat& A,
                             mat& Ay,
                             mat& yA,
                             mat& B,
                             mat& By,
                             mat& C,
                             mat& U3,
                             mat& Ui3,
                             mat& Ud3,
                             mat& Udi3)
{
  mat& Ix = storage.Ix;
  mat& Iz = storage.Iz;
  /* K_{x|z} - l_e K_{xz} K_z^{-1} Z \Lambda_{zx} X' \implies 
   * K_{x|z}^{-1} + l_e K_{x|z}^{-1} K_{xz} K_z^{-1} Z (I - l_e \Lambda_{zx} X' K_{x|z}^{-1} K_{xz} K_z^{-1} Z)^{-1} \Lambda_{zx} X' K_{x|z}^{-1}
   */
  Ui3     = Ix - B * Lxz.t();
  Udi3    = Iz - Lxz.t() * B;
  U3      = arma::inv(Ui3);
  Ud3     = arma::inv(Udi3);
  arma::log_det (ld, ls, U3);
  logdet += ld;
  By      = By * Ud3;
  Ayy    += arma::dot(By * Lxz.t(), yA);
  Ay     += By * Lxz.t() * A;
  A       = U3 * A;
  C       = C * Ud3;
  B       = B * Ud3;
  for (uword ti : tree.descendants.at(ce))
  {
    mat& F        = storage.F.slice(ti);
    mat& G        = storage.G.slice(ti);
    diag.col(ti) += arma::diagvec(G * Ud3 * Lxz.t() * F);
    Y.col(ti)    += arma::trans(By * Lxz.t() * F);
    G             = G * Ud3;
    F             = U3 * F;
  }
}

void Covariance::equationA4 (mat& A,
                             mat& Ay,
                             mat& B,
                             mat& By,
                             mat& C,
                             mat& D,
                             mat& U4,
                             mat& Ui4,
                             mat& Ud4,
                             mat& Udi6)
{
  mat& Ix = storage.Ix;
  /* K_{x|z} - l_e^2 X \Lambda_{xz} Z' K_{z}^{-1} Z \Lambda_{zx} X' \implies 
   *  K_{x|z}^{-1} + l_e^2 K_{x|z}^{-1} X \Lambda_{xz} (I - l_e^2 Z' K_z^{-1} Z \Lambda_{zx} X' K_{x|z}^{-1} X \Lambda_{xz})^{-1} Z' K_z^{-1} Z \Lambda_{zx} X' K_{x|z}^{-1}
   */
  Ud4     = Lxz * D * Udi6 * Lxz.t();
  Ui4     = Ix - Ud4 * A;
  U4      = arma::inv(Ui4);
  arma::log_det (ld, ls, U4);
  logdet += ld;
  Ayy    += arma::dot(Ay * U4 * Ud4, Ay);
  Ay      = Ay * U4;
  A       = A * U4;
  C      += B.t() * U4 * Ud4 * B;
  By     += Ay * Ud4 * B;
  B       = U4.t() * B;
  for (uword ti : tree.descendants.at(ce))
  {
    mat&  F       = storage.F.slice(ti);
    mat&  G       = storage.G.slice(ti);
    diag.col(ti) += arma::diagvec(F.t() * U4 * Ud4 * F);
    Y.col(ti)    += arma::trans(Ay * Ud4 * F);
    G            += F.t() * Ud4 * B + F.t() * U4 * Lxz * D * Udi6;
    F             = U4.t() * F;
  }
}

void Covariance::equationA5 (mat& A,
                             mat& Ay,
                             mat& B,
                             mat& By,
                             mat& C,
                             mat& D,
                             mat& Udi6)
{
  /* K_{xz} + l_e X \Lambda_{xz} Z' */
  Lxz      = Lxz * Udi6.t() * D.t();
  C       += Lxz.t() * B + B.t() * Lxz + Lxz.t() * A * Lxz;
  By      += Ay * Lxz;
  B       += A * Lxz;
}

void Covariance::equationA6 (mat& A,
                             mat& Ay,
                             mat& B,
                             mat& By,
                             mat& C,
                             mat& U6,
                             mat& Ui6,
                             mat& Ud6,
                             mat& Udi6)
{
  /* 
   * K_{z}^{-1} - l_e K_{z}^{-1} Z (I + l_e \Lambda_z Z' K_{z}^{-1} Z)^{-1} \Lambda_z Z' K_{z}^{-1} \implies
   * K_{x|z} + l_e^2 K_{xz} (K_{z}^{-1} - l_e K_{z}^{-1} Z (I + l_e \Lambda_z Z' K_{z}^{-1} Z)^{-1} \Lambda_z Z' K_{z}^{-1}) K_{zx} \implies 
   * K_{x|z}^{-1} - l_e K_{x|z}^{-1} K_{xz} K_{z}^{-1} Z solve(I + l_e \Lambda_z Z' K_{z}^{-1} Z + l_e \Lambda_z Z' K_{z}^{-1} K_{zx} K_{x|z}^{-1} K_{xz} K_{z}^{-1} Z)^{-1} \Lambda_z Z' K_{z}^{-1} K_{zx} K_{x|z}^{-1}
   */
  if (upx)
  {
    Ui6     = Udi6 + Lz * C;
    U6      = arma::inv(Ui6);//Ud6 * arma::inv(Iz + Lz * C * Ud6);
    arma::log_det (ld, ls, U6);
    logdet += ld;
    arma::log_det (ld, ls, Ud6);
    logdet -= ld;
    Ay     -= By * U6 * Lz * B.t();
    Ayy    -= arma::dot(By * U6 * Lz, By);
    By      = By * U6;
    A      -= B * U6 * Lz * B.t();
    B       = B * U6;
    C       = Ud6.t() * C * U6;
    for (uword ti : tree.descendants.at(ce))
    {
      mat&  F       = storage.F.slice(ti);
      mat&  G       = storage.G.slice(ti);
      diag.col(ti) -= arma::diagvec(G * U6 * Lz * G.t());
      Y.col(ti)    -= arma::trans(By * Lz * G.t());
      F            -= B * Lz * G.t();
      G             = G * U6;
    }
  }
}

/* reverse algorithmic differentiation */
void Covariance::traverse_reverse_marginal (const mat& dZ_, double dld_, double dDyy_)
{
  start_reverse_marginal(dZ_, dld_, dDyy_);
  for (e=tree.edges; e>0; --e)
  {
    edge_down_marginal();
    if (is_a_tip)
      derivatives_tip_marginal();
    else
      derivatives_node_marginal();
  }
  derivatives_increment_marginal();
}

void Covariance::start_reverse_marginal (const mat& dZ_, double dld_, double dDyy_)
{
  storage.dY = dZ_;
  dDyy       = dDyy_;
  dld        = dld_;
  storage.clear_derivatives_marginal ();
  dr.zeros();
  dL.zeros();
  dS.zeros();
}

void Covariance::edge_down_marginal (void)
{
  pe       = tree.edge(e-1,0);
  ce       = tree.edge(e-1,1);
  le       = tree.edge_length(e-1);
  is_a_tip = ce < tree.tips;

  missing_data ();
  if (upz)
  {
    Lz  = storage.Lz*le;
    storage.dD.slice(ce)   = storage.dD.slice(pe);
    storage.dDy.slice(ce)  = storage.dDy.slice(pe);
    if (upx)
    {
      Lxz = storage.Lxz*le;
    }
  }
}

void Covariance::derivatives_tip_marginal (void)
{
  /*
   * Di -> D, Dy, Dyy, ld, H, Z
   * Lz,Sz -> Di
   * Lxz,Sxz -> H, Z
   * y -> Dy, Dyy, Z
   */
  mat& D    = storage.D.slice(ce);
  mat& dD   = storage.dD.slice(ce);
  mat& Dy   = storage.Dy.slice(ce);
  mat& dDy  = storage.dDy.slice(ce);
  mat& dH   = storage.dH.slice(ce);
  mat& dLz  = storage.dLz;
  mat& dLxz = storage.dLxz;
  mat& dSz  = storage.dSz;
  mat& dSxz = storage.dSxz;
  mat& Sxz  = storage.Sxz;
  auto dZ   = storage.dY.col(ce);

  if (upz)
  {
    dD       = -D * dD * D +
               -Dy.t() * dDy * D +
               -dDyy * Dy.t() * Dy +
               -dld * D;
    if (upx)
    {
      Lxz        +=  Sxz;
      dD         += -D * Lxz.t() * dH * D +
                    -D * Lxz.t() * dZ * Dy;
      mat tmp     =  dH * D + dZ * Dy;
      dSxz       +=  tmp; 
      dLxz       +=  tmp * le;
      dr.col(ce) +=  D * Lxz.t() * dZ;
    }
    dLz        +=  dD * le; 
    dSz        +=  dD;
    dr.col(ce) +=  D * dDy.t() + 
                   dDyy * 2 * Dy.t();
  }
}

void Covariance::derivatives_node_marginal (void)
{
  /*
   * Ud -> D', Dyy, Dy', H', Z, ld
   * D -> Ud, D', H'
   * H -> H', Z
   * Dy -> Dy', Dyy, Z
   * Lz -> Ud, Dyy, Z
   * Lxz -> H', Z
   */
  mat& Ud6  = storage.Ud6.slice(ce);
  mat& Udi6 = storage.Udi6.slice(ce);
  mat& D    = storage.D.slice(ce);
  mat& dD   = storage.dD.slice(ce);
  mat& Dy   = storage.Dy.slice(ce);
  mat& dDy  = storage.dDy.slice(ce);
  mat& dLz  = storage.dLz;
  mat& dLxz = storage.dLxz;

  if (upz)
  {
    mat dUd   = -D.t() * dD * Ud6.t() +
                 Dy.t() * Dy * Lz * dDyy +
                -Dy.t() * dDy * Ud6.t() +
                -dld * Ud6.t();
    dD        =  dD * Ud6.t();
    dDy       =  dDy * Ud6.t() +
                -dDyy * 2. * Dy * Lz; 
    if (upx)
      for (uword ti : tree.descendants.at(ce))
      {
        mat& H  = storage.H.slice(ti);
        mat& dH = storage.dH.slice(ti);
        auto dZ = storage.dY.col(ti);
        H       = (H - Lxz * D) * Udi6; /* revert */
        dUd    += -Ud6.t() * H.t() * dH * Ud6.t() +
                  -D * dH.t() * Lxz * Ud6.t() +
                   Dy.t() * dZ.t() * H * Lz * Ud6.t() +
                  -Dy.t() * dZ.t() * Lxz * Ud6.t();
        dD     +=  dH.t() * Lxz * Ud6.t();
        dDy    +=  dZ.t() * Lxz * Ud6.t() - dZ.t() * H * Lz * Ud6.t();
        dLz    += -H.t() * dZ * Dy * le;
        dLxz   +=  dZ * Dy * le + dH * D * le;
        dH      =  dH * Ud6.t() - dZ * Dy * Lz;
      }
    dD       +=  Lz * dUd; 
    dLz      +=  dUd * D * Udi6 * le + 
                -dDyy * Dy.t() * Dy * Udi6 * le;
  }
}

void Covariance::derivatives_increment_marginal (void)
{
  dL(storage.tz,storage.tz) = storage.dLz;
  dS(storage.tz,storage.tz) = storage.dSz;
  dL(storage.tx,storage.tz) = storage.dLxz;
  dS(storage.tx,storage.tz) = storage.dSxz;
}

void Covariance::traverse_reverse (const mat& dY_, double dld_, double dAyy_)
{
  start_reverse(dY_, dld_, dAyy_);
  for (e=tree.edges; e>0; --e)
  {
    edge_down();
    if (is_a_tip)
      derivatives_tip();
    else
      derivatives_node();
  }
  derivatives_increment();
}

void Covariance::start_reverse (const mat& dY_, const double dld_, const double dAyy_)
{
  storage.dY = dY_;
  dAyy       = dAyy_;
  dld        = dld_;
  storage.clear_derivatives ();
  dy.zeros();
  dL.zeros();
  dS.zeros();
}

void Covariance::edge_down (void)
{
  pe = tree.edge(e-1,0);
  ce = tree.edge(e-1,1);
  le = tree.edge_length(e-1);
  is_a_tip = ce < tree.tips;

  // clear temporaries
  storage.Gt.zeros();
  storage.dyA.zeros();
  storage.yA.zeros();

  // starting point: derivatives of parent
  missing_data();
  storage.dA.slice(ce) = storage.dA.slice(pe);
  storage.dAy.slice(ce) = storage.dAy.slice(pe);
  storage.dB.slice(ce) = storage.dB.slice(pe);
  storage.dBy.slice(ce) = storage.dBy.slice(pe);
  storage.dC.slice(ce) = storage.dC.slice(pe);
  storage.dD.slice(ce) = storage.dD.slice(pe);

  // get Lx/etc.
  Lx  = storage.Lx*le;
  Lz  = storage.Lz*le;
  Lxz = storage.Lxz*le;
}

void Covariance::derivatives_tip (void)
{
  /* quantities */
  mat& A = storage.A.slice(ce);
  mat& Ay = storage.Ay.slice(ce);
  mat& B = storage.B.slice(ce);
  mat& By = storage.By.slice(ce);
  mat& C = storage.C.slice(ce);
  mat& D = storage.D.slice(ce);
  mat& F = storage.F.slice(ce);
  mat& G = storage.G.slice(ce);
  mat& Sx = storage.Sx;
  mat& Sz = storage.Sz;
  mat& Sxz = storage.Sxz;
  /* derivatives */
  mat& dA = storage.dA.slice(ce);
  mat& dAy = storage.dAy.slice(ce);
  mat& dB = storage.dB.slice(ce);
  mat& dBy = storage.dBy.slice(ce);
  mat& dC = storage.dC.slice(ce);
  mat& dD = storage.dD.slice(ce);
  mat& dF = storage.dF.slice(ce);
  mat& dG = storage.dG.slice(ce);
  mat& dLx = storage.dLx;
  mat& dLz = storage.dLz;
  mat& dLxz = storage.dLxz;
  mat& dLzx = storage.dLzx;
  mat& dSx = storage.dSx;
  mat& dSz = storage.dSz;
  mat& dSxz = storage.dSxz;
  mat& dSzx = storage.dSzx;
  auto dm = dy.col(ce);
  auto dY = storage.dY.col(ce);

  /* */
  dA = -A.t() * dA * A.t() - dld * A.t() - dAyy * Ay.t() * Ay -
        Ay.t() * dAy * A.t() - A.t() * dF * A.t() - Ay.t() * dY.t() * A.t();
  dm = A * dAy.t() + dAyy * 2 * Ay.t() + A * dY;

  if (upz)
  {
    dm += B * dBy.t();
    Lxz += Sxz;
    dA += -(A.t() * dB + Ay.t() * dBy + A.t() * Lxz * D.t() * dC) * D.t() * Lxz.t() * A.t() -
          A.t() * dG * D.t() * Lxz.t() * A.t();

    dD = -D.t() * dD * D.t() +
        D.t() * Lxz.t() * dA * Lxz * D.t() -
        D.t() * Lxz.t() * A.t() * dB * D.t() -
        D.t() * Lxz.t() * Ay.t() * dBy * D.t() -
        D.t() * (dC + dC.t()) * D * Lxz.t() * A * Lxz * D.t() -
        D.t() * Lxz.t() * A.t() * dG * D.t();

    mat tmp = A.t() * dG * D.t();
    dSxz += tmp;
    dLxz += tmp * le;

    dLz += dD * le;
    dSz += dD;

    tmp = -dA * Lxz * D.t() +
      A.t() * dB * D.t() +
      Ay.t() * dBy * D.t() +
      A * Lxz * D * dC * D;
    dSxz += tmp;
    dLxz += tmp*le;

    tmp = -D.t() * Lxz.t() * dA + D.t() * dC * D * Lxz.t() * A;
    dSzx += tmp;
    dLzx += tmp * le;
  }

  dLx += dA * le;
  dSx += dA;
}

void Covariance::derivatives_node (void)
{
  mat& A = storage.A.slice(ce);
  mat& Ay = storage.Ay.slice(ce);
  mat& B = storage.B.slice(ce);
  mat& By = storage.By.slice(ce);
  mat& C = storage.C.slice(ce);
  mat& D = storage.D.slice(ce);
  mat& U1 = storage.U1.slice(ce);
  mat& U2 = storage.U2.slice(ce);
  mat& U3 = storage.U3.slice(ce);
  mat& U4 = storage.U4.slice(ce);
  mat& U6 = storage.U6.slice(ce);
  mat& Ui1 = storage.Ui1.slice(ce);
  mat& Ui2 = storage.Ui2.slice(ce);
  mat& Ui3 = storage.Ui3.slice(ce);
  mat& Ui4 = storage.Ui4.slice(ce);
  mat& Ui6 = storage.Ui6.slice(ce);
  mat& Ud2 = storage.Ud2.slice(ce);
  mat& Ud3 = storage.Ud3.slice(ce);
  mat& Ud4 = storage.Ud4.slice(ce);
  mat& Ud6 = storage.Ud6.slice(ce);
  mat& Udi2 = storage.Udi2.slice(ce);
  mat& Udi3 = storage.Udi3.slice(ce);
  mat& Udi6 = storage.Udi6.slice(ce);
  mat& dA = storage.dA.slice(ce);
  mat& dAy = storage.dAy.slice(ce);
  mat& dyA = storage.dyA;
  mat& dB = storage.dB.slice(ce);
  mat& dBy = storage.dBy.slice(ce);
  mat& dC = storage.dC.slice(ce);
  mat& dD = storage.dD.slice(ce);
  mat& dLx = storage.dLx;
  mat& dLz = storage.dLz;
  mat& dLxz = storage.dLxz;
  mat& dLzx = storage.dLzx;

  if (upz)
  {
    derivatives_equationA6(A, Ay, B, By, C, D, U6, Ui6, Ud6, Udi6, 
                           dLz, dA, dAy, dB, dBy, dC, dD);
  }
  if (upz && upx)
  {
    derivatives_equationA5(A, Ay, B, By, C, D,
                           dLxz, dLzx, dA, dAy, dB, dBy, dC, dD);
    derivatives_equationA4(A, Ay, B, By, C, D, U4, Ui4, Ud4, Ud6,
                           dLxz, dLzx, dA, dAy, dB, dBy, dC, dD);
    derivatives_equationA3(A, Ay, B, By, C, U3, Ui3, Ud3, Udi3,
                           dLzx, dA, dAy, dyA, dB, dBy, dC);
    derivatives_equationA2(A, Ay, B, By, C, U2, Ui2, Ud2, Udi2,
                           dLxz, dA, dAy, dyA, dB, dBy, dC);
  }
  if (upx)
  {
    derivatives_equationA1(A, Ay, B, By, C, U1, Ui1,
                           dLx, dA, dAy, dB, dBy, dC);
  }
}

void Covariance::derivatives_increment (void)
{
  dL(storage.tx,storage.tx) = storage.dLx;
  dS(storage.tx,storage.tx) = storage.dSx;
  dL(storage.tx,storage.tz) = storage.dLxz + storage.dLzx.t();
  dS(storage.tx,storage.tz) = storage.dSxz + storage.dSzx.t();
  dL(storage.tz,storage.tz) = storage.dLz;
  dS(storage.tz,storage.tz) = storage.dSz;
}

void Covariance::derivatives_equationA6 (mat& A,
                                         mat& Ay,
                                         mat& B,
                                         mat& By,
                                         mat& C,
                                         mat& D,
                                         mat& U6,
                                         mat& Ui6,
                                         mat& Ud6,
                                         mat& Udi6,
                                         mat& dLz,
                                         mat& dA,
                                         mat& dAy,
                                         mat& dB,
                                         mat& dBy,
                                         mat& dC,
                                         mat& dD)
{
  mat& Iz   = storage.Iz;
  mat  Ao,//could make these part of storage cuz the size won't vary
       Ayo,
       Bo,
       Byo,
       Co,
       Do,
       dUd,
       dU;

  /* previous states */
  if (upx)
  {
    Co  = Udi6.t() * C * Ui6;
    Bo  = B * Ui6;
    Ao  = A + B * Lz * Bo.t();
    Ayo = Ay + By * Lz * Bo.t();
    Byo = By * Ui6;
  }
  Do    = D * Udi6;

  /* derivatives */
  dUd = -D * dD * Ud6.t();
  if (upx)
  {
    dUd += -C * dC.t() * Ud6.t() + dld * Ud6.t();
    dU   = ((By.t()*dAy + B.t()*dA)*B + dAyy*By.t() * By)*Lz -
            (B.t()*dB + C.t()*dC + By.t()*dBy + dld*Iz)*U6.t();
    dB   = dB*U6.t() - ((dA.t() + dA)*B + dAy.t()*By)*Lz;
    dBy  = dBy*U6.t() - (dAy*B + dAyy*2*By)*Lz;

    for (uword ti : tree.descendants.at(ce))
    {
      mat& dF = storage.dF.slice(ti);
      mat& dG = storage.dG.slice(ti);
      mat& F  = storage.F.slice(ti);
      mat& G  = storage.G.slice(ti);
      auto dY = storage.dY.col(ti);
      mat  Go      = G * Ui6,
           Fo      = F + B * Lz * Go.t();

      dU  += B.t() * dF * G * Lz -
               G.t() * dG * U6.t() +
               By.t() * dY.t() * G * Lz;
      dB  += -dF * G * Lz;
      dBy += -dY.t() * G * Lz;
      dLz += -B.t() * dF * Go * le -
               By.t() * dY.t() * Go * le;
      dG   = dG * U6.t() - 
               dF.t() * B * Lz -
               dY * By * Lz;

      /* revert linear part */
      F    = Fo;
      G    = Go;
    }
    dC     = Ud6 * dC * U6.t() + Lz * dU;
    dLz   += -(By.t() * dAy + B.t() * dA) * Bo * le -
              dAyy * By.t() * Byo * le;
  }
  dD = dD * Ud6.t() + Lz * dUd;
  if (upx) dD += Lz * dU;
  dLz += dUd * Do.t() * le;
  if (upx) dLz += dU * (Do + Co) * le;

  /* revert */
  if (upx)
  {
    A = Ao;
    Ay = Ayo;
    B = Bo;
    By = Byo;
    C = Co;
  }
  D = Do;
}

void Covariance::derivatives_equationA5 (mat& A,
                                         mat& Ay,
                                         mat& B,
                                         mat& By,
                                         mat& C,
                                         mat& D,
                                         mat& dLxz,
                                         mat& dLzx,
                                         mat& dA,
                                         mat& dAy,
                                         mat& dB,
                                         mat& dBy,
                                         mat& dC,
                                         mat& dD)
{
  mat  DLxz = D * Lxz.t();

  /* revert */
  B        += -A * DLxz.t();
  By       += -Ay * DLxz.t();
  C        += -DLxz * (B + A * DLxz.t()) - B.t() * DLxz.t();

  /* derivatives */
  for (uword ti : tree.descendants.at(ce))
  {
    mat& dF = storage.dF.slice(ti);
    mat& dG = storage.dG.slice(ti);
    mat& F  = storage.F.slice(ti);
    mat& G  = storage.G.slice(ti);
    /* revert */
    G      += -F.t() * DLxz.t();

    dF     += DLxz.t() * dG.t();
    dLxz   += F * dG * D * le;
    dD     += Lxz.t() * F * dG;
  }
  dA       += (DLxz.t()*dC + dB)*DLxz;
  dD       += Lxz.t() * (A*dB + Ay.t()*dBy) +
               (dC.t() + dC)*(B.t() + DLxz*A)*Lxz;
  dAy      += dBy*DLxz;
  dLxz     += A * (dB + DLxz.t()*dC)*D * le +
               (B*dC + Ay.t()*dBy)*D * le;
  dLzx     += D*dC*(B.t() + DLxz*A) * le;
  dB       += DLxz.t()*(dC + dC.t());
}

void Covariance::derivatives_equationA4 (mat& A,
                                         mat& Ay,
                                         mat& B,
                                         mat& By,
                                         mat& C,
                                         mat& D,
                                         mat& U4,
                                         mat& Ui4,
                                         mat& Ud4,
                                         mat& Ud6,
                                         mat& dLxz,
                                         mat& dLzx,
                                         mat& dA,
                                         mat& dAy,
                                         mat& dB,
                                         mat& dBy,
                                         mat& dC,
                                         mat& dD)
{
  /* old states */
  mat Bo    = Ui4.t() * B,
      Byo   = By - Ay * Ud4 * Bo,
      Co    = C - B.t() * Ud4 * Bo,
      Ao    = A * Ui4,
      Ayo   = Ay * Ui4;

  mat dU = -Ay.t() * dAy * U4.t() -
            A * dA * U4.t() -
            B * dC * B.t() * Ud4.t() -
            Ay.t() * dBy * B.t() * Ud4.t() -
            B * dB.t() * U4.t() -
            dld * U4.t() -
            dAyy * Ay.t() * Ay * Ud4.t();

  for (uword ti : tree.descendants.at(ce))
  {
    mat& F = storage.F.slice(ti);
    mat& G = storage.G.slice(ti);
    mat& dF = storage.dF.slice(ti);
    mat& dG = storage.dG.slice(ti);
    auto dY = storage.dY.col(ti); 
    dU += -Ay.t() * dY.t() * F.t() * Ud4.t() -
      F * dG * B.t() * Ud4.t() -
      F * dF.t() * U4.t();
    dF = U4 * dF + Ud4 * B * dG.t() + 
          Ud4.t() * Ay.t() * dY.t(); /* try to move this to end */

    /* revert */
    F = Ui4.t() * F;
    G += -F.t() * Ud4 * B;
  }

  dA = dA * U4.t() - Ud4.t() * dU;
  dAy = dAy * U4.t() + dBy * B.t() * Ud4.t() +
        dAyy * 2 * Ay * Ud4;
  dB = U4 * dB + Ud4.t() * B * (dC + dC.t()) +
        Ud4.t() * Ay.t() * dBy;

  dD += -Lxz.t() * dU * Ao.t() * Lxz + 
        Lxz.t() * B * dC * Bo.t() * Lxz +
        Lxz.t() * Ay.t() * dBy * Bo.t() * Lxz +
        dAyy * Lxz.t() * Ay.t() * Ayo * Lxz;
  dLxz += -dU * Ao.t() * Lxz * D * le +
           B * dC * Bo.t() * Lxz * D * le +
           Ay.t() * dBy * Bo.t() * Lxz * D * le +
           dAyy * Ay.t() * Ayo * Lxz * D * le;
  dLzx += -D * Lxz.t() * dU * Ao * le +
          D * Lxz.t() * B * dC * Bo.t() * le +
          D * Lxz.t() * Ay.t() * dBy * Bo.t() * le +
          dAyy * D * Lxz.t() * Ay.t() * Ayo * le;
  for (uword ti : tree.descendants.at(ce))
  {
    auto dY = storage.dY.col(ti);
    mat& F = storage.F.slice(ti); /* can't I revert here istead? */
    mat& dG = storage.dG.slice(ti);
    dAy += dY.t() * F.t() * Ud4.t() * U4.t();
    dB += Ud4.t() * U4.t() * F * dG;
    dD += Lxz.t() * Ay.t() * dY.t() * F.t() * Lxz +
          Lxz.t() * U4.t() * F * dG * Bo.t() * Lxz;
    dLxz += Ay.t() * dY.t() * F.t() * Lxz * D * le +
            U4.t() * F * dG * Bo.t() * Lxz * D * le;
    dLzx += D * Lxz.t() * U4.t() * F * dG * Bo.t() * le +
            D * Lxz.t() * Ay.t() * dY.t() * F.t() * le;
  }

  /* make reversion permanent */
  A = Ao;
  Ay = Ayo;
  B = Bo;
  By = Byo;
  C = Co;
  /* un-revert D */
  D = D * Ud6;
}

void Covariance::derivatives_equationA3 (mat& A,
                                         mat& Ay,
                                         mat& B,
                                         mat& By,
                                         mat& C,
                                         mat& U3,
                                         mat& Ui3,
                                         mat& Ud3,
                                         mat& Udi3,
                                         mat& dLzx,
                                         mat& dA,
                                         mat& dAy,
                                         mat& dyA,
                                         mat& dB,
                                         mat& dBy,
                                         mat& dC)
{
  /* old states */
  mat Bo = B * Udi3,
      Co = C * Udi3,
      Ao = Ui3 * A,
      yA = Ay.t(),
      yAo = Ui3 * Ay.t(),
      Ayo = Ay - By * Lxz.t() * Ao,
      Byo = By * Udi3;
  for (uword ti : tree.descendants.at(ce))
  {
    mat& F = storage.F.slice(ti);
    mat& G = storage.G.slice(ti);
    mat& Gt = storage.Gt.slice(ti);

    /* revert */
    F = Ui3 * F;
    Gt = G.t() - C * Lxz.t() * F;
    G = G * Udi3;
  }

  mat dU = -U3.t() * dA * A.t() - dld * U3.t();
  mat dUd = -By.t() * dBy * Ud3.t() -
            By.t() * dAy * A.t() * Lxz - 
            C.t() * dC * Ud3.t() -
            B.t() * dB * Ud3.t() -
            dAyy * By.t() * yA.t() * Lxz; 

  for (uword ti : tree.descendants.at(ce))
  { 
    auto dY = storage.dY.col(ti);
    mat& dF = storage.dF.slice(ti);
    mat& F = storage.F.slice(ti);
    mat& G = storage.G.slice(ti);
    mat& dG = storage.dG.slice(ti);
    dU += -U3.t() * dF * F.t() * U3.t();
    dUd += -By.t() * dY.t() * F.t() * Lxz * Ud3.t() -
            Ud3.t() * G.t() * dG * Ud3.t();
    dF = U3.t() * dF + Lxz * By.t() * dY.t();
    dG = dG * Ud3.t();
  }
  dBy = dBy * Ud3.t() +
        dAy * A.t() * Lxz + 
        dAyy * yA.t() * Lxz; 
  dB = dB * Ud3.t() - dU * Lxz - Lxz * dUd;
  dyA = dAyy * Lxz * By.t();
  dA = U3.t() * dA + Lxz * By.t() * dAy;
  dC = dC * Ud3.t();
  dLzx += -Bo.t() * dU * le -
            dUd * Bo.t() * le +
            By.t() * dAy * Ao.t() * le +
            dAyy * By.t() * yAo.t() * le;
  for (uword ti : tree.descendants.at(ce))
  {
    mat& F = storage.F.slice(ti);
    auto dY = storage.dY.col(ti);
    dBy += dY.t() * F.t() * Lxz * Ud3.t();
    dLzx += By.t() * dY.t() * F.t() * le;
  }

  /* revert here */
  A = Ao;
  Ay = Ayo;
  B = Bo;
  By = Byo;
  C = Co;
}

void Covariance::derivatives_equationA2 (mat& A,
                                         mat& Ay, 
                                         mat& B,
                                         mat& By,
                                         mat& C,
                                         mat& U2,
                                         mat& Ui2,
                                         mat& Ud2,
                                         mat& Udi2,
                                         mat& dLxz,
                                         mat& dA,
                                         mat& dAy,
                                         mat& dyA,
                                         mat& dB,
                                         mat& dBy,
                                         mat& dC)
{
  /* old states */
  mat  Co   = Udi2 * C,
       Bo   = B - A * Lxz * Co,
       Byo  = By - Ay * Lxz * Co,
       Ayo  = Ay * Ui2,
       Ao   = A * Ui2;
  for (uword ti : tree.descendants.at(ce))
  {
    mat& F = storage.F.slice(ti);
    mat& G = storage.G.slice(ti);
    mat& Gt = storage.Gt.slice(ti);
    Gt = Udi2 * Gt;
    G = Gt.t();
    F += -A * Lxz * Gt;
  }

  mat dU = -A.t() * dA * U2.t() -
            Ay.t() * dAy * U2.t() -
            Ay.t() * dBy * C.t() * Lxz.t() - 
            A.t() * dB * C.t() * Lxz.t() -
            A.t() * dyA * Byo * Lxz.t() * U2.t() -
            dld * U2.t() -
            dAyy * Ay.t() * Byo * Lxz.t() * U2.t();
  mat dUd = -Ud2.t() * dC * C.t();
  for (uword ti : tree.descendants.at(ce))
  {
    auto dY = storage.dY.col(ti);
    mat& G = storage.G.slice(ti);
    mat& F = storage.F.slice(ti);
    mat& dG = storage.dG.slice(ti);
    mat& dF = storage.dF.slice(ti);
    dU += -Ay.t() * dY.t() * G * Lxz.t() * U2.t() -
          A.t() * dF * G * Lxz.t() * U2.t();
    dUd += -Ud2.t() * Lxz.t() * F * dG * C.t();
  }
  dA = dA * U2.t() + dB * C.t() * Lxz.t() + 
        dyA * Byo * Lxz.t() * U2.t();

  dAy = dAy * U2.t() +
        dBy * C.t() * Lxz.t() + 
        dyA.t() + dAyy * Byo * Lxz.t() * U2.t();

  dC = Ud2.t() * dC + Lxz.t() * A.t() * dB +
        Lxz.t() * Ay.t() * dBy;

  dLxz += A.t() * dB * Co.t() * le +
          Ay.t() * dBy * Co.t() * le -
          dU * Bo * le -
          Bo * dUd * le +
          A.t() * dyA * Byo * le +
          dAyy * Ay.t() * Byo * le;

  dB += -dU.t() * Lxz - Lxz * dUd.t();

  dBy += dyA.t() * A * Lxz + dAyy * Ay * Lxz;

  for (uword ti : tree.descendants.at(ce))
  {
    auto dY = storage.dY.col(ti);
    mat& G = storage.G.slice(ti);
    mat& F = storage.F.slice(ti);
    mat& dG = storage.dG.slice(ti);
    mat& dF = storage.dF.slice(ti);
    dAy += dY.t() * G * Lxz.t() * U2.t();
    dA += dF * G * Lxz.t() * U2.t();
    dC += Ud2.t() * Lxz.t() * F * dG;
    dLxz += Ay.t() * dY.t() * G * le +
            F * dG * C.t() * le +
            A.t() * dF * G * le;
    dF += Lxz * C * dG.t();
    dG += dY * Ay * Lxz + 
           (dF - Lxz * C * dG.t()).t() * A * Lxz;
  }
  /* revert to old states */
  A = Ao;
  Ay = Ayo;
  B = Bo;
  By = Byo;
  C = Co;
}

void Covariance::derivatives_equationA1 (mat& A,
                                         mat& Ay,
                                         mat& B,
                                         mat& By,
                                         mat& C,
                                         mat& U1,
                                         mat& Ui1,
                                         mat& dLx,
                                         mat& dA,
                                         mat& dAy,
                                         mat& dB,
                                         mat& dBy,
                                         mat& dC)
{
  mat& Ix   = storage.Ix;
  /* old states */
  mat  Ayo  = Ay * Ui1,
       Ao   = A * Ui1,
       Bo,
       Co,
       Byo;
  if (upz)
  {
    Bo      = Ui1.t() * B;
    Co      = C + B.t() * Lx * Bo;
    Byo     = By + Ay * Lx * Bo;
  }

  /* derivatives */
  mat dU    = -(A.t() * dA + dld * Ix + Ay.t() * dAy) * U1.t() +
                Ay.t() * Ay * Lx * dAyy;
  dA        = dA * U1.t();
  dAy       = dAy * U1.t() - 
                dAyy * 2 * Ay * Lx;
  for (uword ti : tree.descendants.at(ce))
  {
    mat& F  = storage.F.slice(ti);
    mat& dF = storage.dF.slice(ti);
    auto dY = storage.dY.col(ti);
    dU     += Ay.t() * dY.t() * F.t() * Lx -
               F * dF.t() * U1.t();
    dF      = U1 * dF - Lx * Ay.t() * dY.t();
  }
  if (upz)
  {
    dU     += (Ay.t() * dBy + B * dC) * B.t() * Lx.t() -
                B * dB.t() * U1.t();
    for (uword ti : tree.descendants.at(ce))
    {
      mat& F   = storage.F.slice(ti);
      mat& dF  = storage.dF.slice(ti);
      mat& dG  = storage.dG.slice(ti);
      dU      += B * dG.t() * F.t() * Lx;
      dF      += -Lx * B * dG.t();
    }
    dAy    += -dBy * B.t() * Lx.t();
    dB      = U1 * dB - 
                Lx.t() * (Ay.t() * dBy + B * (dC + dC.t()));
    dLx    += -(Ay.t() * dBy + B * dC) * Bo.t() * le;
  }
  dA       += Lx.t() * dU;
  dLx      += dU * Ao.t() * le - dAyy * Ay.t() * Ayo * le;
  for (uword ti : tree.descendants.at(ce))
  {
    mat& F  = storage.F.slice(ti);
    mat& G  = storage.G.slice(ti);
    mat& dF = storage.dF.slice(ti);
    mat& dG = storage.dG.slice(ti);
    auto dY = storage.dY.col(ti);
    dAy    += -dY.t() * F.t() * Lx;
    if (upz)
      dB   += -Lx * F * dG;
    //revert
    F       = Ui1.t() * F;
    if (upz)
      G = G + F.t() * Lx * B;
    dLx    += -Ay.t() * dY.t() * F.t() * le;
    if (upz)
      dLx  += -F * dG * B.t() * le;
  }
}

void Covariance::parameters (const mat& Lambda, const mat& Sigma, const mat& mean)
{
  storage.update(Lambda, Sigma, mean);
  /* TODO: any need to clean house here? */
}

} /* end namespace epee */

