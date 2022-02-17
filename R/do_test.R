sim_polytomy <- function(seed, n, tol)
{
  library(ape)
  set.seed(seed)
  tree <- di2multi(reorder(rcoal(n),"postorder"), tol)
  tree$edge.length <- tree$edge.length/vcv(tree)[1,1] # scale height to 1
  tree
}

test_design <- function(seed=1, n=5, mx=3, mz=2, missing=c(0,0), tol=0)
{
  library(ape)
  set.seed(seed)
  tr    <- sim_polytomy(seed, n, tol)
  typ   <- sample(c(rep(1,mx),rep(0,mz)))
  y     <- matrix(NA,mx+mz,n)
  y[typ==1,] <- sample(c(0,1), n*mx, replace=TRUE)
  y[typ==0,] <- rnorm(n*mz)
  miss       <- matrix(sample(c(TRUE,FALSE), n*(mx+mz), replace=TRUE, prob=c(missing[1],1-missing[1])),mx+mz,n)
  y[miss]    <- NA

  cov <- list()
  beta <- list()
  linear <- list()
  linear2 <- list()
  for(i in 1:(mx+mz))
  {
    ncov <- rpois(1, 1)
    ccc <- matrix(rnorm(n*ncov), n, ncov)
    missx <- matrix(sample(c(TRUE,FALSE), n*(ncov), replace=TRUE, 
                    prob=c(missing[2],1-missing[2])),ncov,n)
    ccc[missx] <- NA
    cov[[i]] <- cbind(1, ccc)
    beta[[i]] <- rnorm(ncov+1)
    linear[[i]] <- cov[[i]] %*% beta[[i]]
    linear[[i]][miss[i,]] <- NA
    linear2[[i]] <- t(ifelse(is.na(linear[[i]]), 0, linear[[i]])) %*% 
                      ifelse(is.na(cov[[i]]),0,cov[[i]])
  }
  linearMat <- t(Reduce(cbind, linear))

  Design = Design$new(y, cov, typ, tr$edge-1, tr$edge.length)
  return(list(
              obj=Design,
              linear=list(epee=Design$linear(beta),r=linearMat),
              adjoint=list(epee=Design$linear2(linearMat),r=linear2)))
}

test_constraints <- function(seed=1, n=5, do_constraints=TRUE)
{
  .random_constraints <- function(n)
  {
    out <- matrix(NA, n, n)
    out[lower.tri(out, diag=TRUE)] <- sample(c(TRUE,FALSE), replace=TRUE, n*(n+1)/2) 
    out
  }

  library(numDeriv)
  set.seed(seed)

  constraints <- .random_constraints(n)

  l <- rnorm(sum(!constraints & !is.na(constraints)))
  dL <- matrix(rnorm(n^2), n, n)

  zeros <- which(constraints & !is.na(constraints), arr.ind=TRUE)-1
  cons <- Constraints$new (zeros, matrix(NA,0,0), n)
  Lambda2 <- cons$apply(l, cons$Lambda)
  dd3 <- cons$apply_d(Lambda2, dL, cons$Lambda)

  dd <- as.vector(dL) %*% jacobian(function(x)
      {
        cons2 <- Constraints$new (zeros, matrix(NA,0,0), n)
        Lambda3 <- cons2$apply(x, cons2$Lambda)
        Lambda3 %*% t(Lambda3)
      }, l)

  list(list(epee=Lambda2%*%t(Lambda2)),
       list(num=dd, 
            epee=dd3), constraints)
}

test_covar <- function(seed=1, n=5, mx=3, mz=2, miss=0, tol=0)
{ 
  library(ape)
  set.seed(seed)
  tr    <- sim_polytomy(seed, n, tol)
  y     <- matrix(rnorm(n*(mx+mz)),mx+mz,n)
  miss  <- matrix(sample(c(TRUE,FALSE), n*(mx+mz), replace=TRUE, prob=c(miss,1-miss)),mx+mz,n)
  y[miss] <- NA
  typ   <- sample(c(rep(1,mx),rep(0,mz)))
  yc    <- y[typ==1,,drop=FALSE] 
  ym    <- y[typ==0,,drop=FALSE] 
  d     <- matrix(rexp(n*mx, 1),mx,n)
  d[is.na(yc)] <- NA
  L     <- as.matrix(rWishart(1, mx+mz, diag(mx+mz))[,,1])
  S     <- as.matrix(rWishart(1, mx+mz, diag(mx+mz))[,,1])
  eta   <- matrix(rnorm((mx+mz)*n), mx+mz, n)
  eta[miss] <- NA
  em    <- eta[typ==0,,drop=FALSE]
  ec    <- eta[typ==1,,drop=FALSE]

  Cova <- Covariance$new(tr$edge-1, tr$edge.length, y, typ)
  timing <- system.time({
    Cova$parameters(L, S, eta)
    Cova$traverse_tree_marginal()
    Cova$traverse_tree(yc, d)})

  dx    <- !is.na(as.vector(t(yc)))
  dz    <- !is.na(as.vector(t(ym)))
  K0    <- kronecker(L[typ==0,typ==0], vcv(tr)) + kronecker(S[typ==0,typ==0], diag(n))
  K1    <- kronecker(L[typ==1,typ==1], vcv(tr)) + kronecker(S[typ==1,typ==1], diag(n))
  K10   <- kronecker(L[typ==1,typ==0], vcv(tr)) + kronecker(S[typ==1,typ==0], diag(n))
  if(mz>0)
    K     <- K1[dx,dx] - K10[dx,dz] %*% solve(K0[dz,dz]) %*% t(K10[dx,dz]) + 
              diag(as.vector(t(d))[dx])
  else
    K     <- K1[dx,dx] + diag(as.vector(t(d))[dx])

  cat("timing: ", timing["elapsed"], "\n")
  list("epee" = 
       list(Ayy=Cova$Ayy, logdet=Cova$logdet, 
            Dyy=Cova$Dyy, logdet_marginal=Cova$logdet_marginal,
            Y=Cova$Y, Z=Cova$Z),
       "inv" =
         list(Ayy=as.vector(t(yc))[dx]%*%solve(K)%*%as.vector(t(yc))[dx],
              logdet=determinant(K)$modulus,
              Dyy=if(mz>0) as.vector(t(ym-em))[dz]%*%
                           solve(K0[dz,dz])%*%
                           as.vector(t(ym-em))[dz] else NULL,
              logdet_marginal=if(mz>0) determinant(K0[dz,dz])$modulus else NULL,
              Y=as.vector(t(yc))[dx]%*%solve(K),
              Z=if(mz>0) as.vector(t(ec))[dx] + K10[dx,dz]%*%solve(K0[dz,dz])%*%
                         as.vector(t(ym-em))[dz] else NULL
              )
       )
}

test_covar_deriv <- function(seed=1, n=5, mx=3, mz=2, miss=0, tol=0, raw=FALSE)
{
  library(ape); library(microbenchmark); library(numDeriv)
  set.seed(seed)
  tr    <- sim_polytomy(seed, n, tol)
  y     <- matrix(rnorm(n*(mx+mz)),mx+mz,n)
  miss  <- matrix(sample(c(TRUE,FALSE), n*(mx+mz), replace=TRUE, prob=c(miss,1-miss)),mx+mz,n)
  y[miss] <- NA
  typ   <- sample(c(rep(1,mx),rep(0,mz)))
  yc    <- y[typ==1,,drop=FALSE] 
  ym    <- y[typ==0,,drop=FALSE] 
  d     <- matrix(rexp(n*mx, 1),mx,n)
  d[is.na(yc)] <- NA
  L     <- as.matrix(rWishart(1, mx+mz, diag(mx+mz))[,,1])
  S     <- as.matrix(rWishart(1, mx+mz, diag(mx+mz))[,,1])

  eta   <- matrix(rnorm((mx+mz)*n), mx+mz, n)
  eta[miss] <- NA
  em    <- eta[typ==0,,drop=FALSE]
  ec    <- eta[typ==1,,drop=FALSE]

  dY    <- matrix(rnorm(mx*n),mx,n)
  dx    <- !is.na(as.vector(t(yc)))
  dz    <- !is.na(as.vector(t(ym)))
  dY[is.na(yc)] = 0;

  Cova <- Covariance$new(tr$edge-1, tr$edge.length, y, typ)
  Cova$parameters(L, S, eta)
  Cova$traverse_tree_marginal()
  yyc <- yc + Cova$Z # conditional mean
  Cova$traverse_tree(yyc, d)
  Cova$traverse_reverse(dY, 1, -1)
  dL <- Cova$dL; dS <- Cova$dS; dy <- Cova$dy
  Cova$traverse_reverse_marginal(dy, 1, -1)
  dL <- dL + Cova$dL; dS <- dS + Cova$dS; dr <- Cova$dr

  sym <- function(X)
  {
    X = X + t(X)
    diag(X) = diag(X)/2
    X[upper.tri(X)] <- NA
    X
  }

  get_mat <- function(L, S, type = c("quad", "det", "linear", "mquad", "mdet", "mlinear"))
  {
    type <- match.arg(type)
    Cova2 <- Covariance$new(tr$edge-1, tr$edge.length, y, typ)
    Cova2$parameters(L, S, eta)
    Cova2$traverse_tree_marginal()
    ymn <- yc + Cova$Z
    Cova2$traverse_tree(ymn, d)

    if (type == "quad")
      return(Cova2$Ayy)
    else if (type == "det")
      return(Cova2$logdet)
    else if (type == "linear")
      return(as.vector(t(Cova2$Y))[dx])
    else if (type == "mquad")
      return (Cova2$Dyy)
    else if (type == "mdet")
      return (Cova2$logdet_marginal)
    else if (type == "mlinear")
      return (as.vector(t(Cova2$Z))[dx])
  }

  get_vec_marg <- function(ym, L, S, type = c("mquad", "mdet", "mlinear"))
  {
    type       <- match.arg(type)
    y[typ==0,] <- ym
    y[miss]    <- NA
    Cova2 <- Covariance$new(tr$edge-1, tr$edge.length, y, typ)
    Cova2$parameters(L, S, eta)
    Cova2$traverse_tree_marginal()

    if (type == "mquad")
      return (Cova2$Dyy)
    else if (type == "mdet")
      return (Cova2$logdet_marginal)
    else if (type == "mlinear")
      return (as.vector(t(Cova2$Z))[dx])
  }

  get_vec <- function(yc, L, S, type = c("quad", "det", "linear"))
  {
    type       <- match.arg(type)
    y[typ==1,] <- yc 
    y[miss]    <- NA
    Cova2 <- Covariance$new(tr$edge-1, tr$edge.length, y, typ)
    Cova2$parameters(L, S, eta)
    Cova2$traverse_tree_marginal()
    ymn <- yc + Cova$Z
    Cova2$traverse_tree(ymn, d)

    if (type == "quad")
      return (Cova2$Ayy)
    else if (type == "det")
      return (Cova2$logdet)
    else if (type == "linear")
      return (as.vector(t(Cova2$Y))[dx])
  }

  dL2 <- matrix(jacobian(function(L) -get_mat(L, S, "quad"), L), mx+mz, mx+mz) +
         matrix(jacobian(function(L) get_mat(L, S, "det"), L), mx+mz, mx+mz) +
         matrix(as.vector(t(dY))[dx] %*% jacobian(function(L) get_mat(L, S, "linear"), L), mx+mz, mx+mz) +
        matrix(jacobian(function(L) -get_mat(L, S, "mquad"), L), mx+mz, mx+mz) +
        matrix(jacobian(function(L) get_mat(L, S, "mdet"), L), mx+mz, mx+mz) +
        matrix(as.vector(t(dy))[dx] %*% jacobian(function(L) get_mat(L, S, "mlinear"), L), mx+mz, mx+mz) 

  dS2 <- matrix(jacobian(function(S) -get_mat(L, S, "quad"), S), mx+mz, mx+mz) +
        matrix(jacobian(function(S) get_mat(L, S, "det"), S), mx+mz, mx+mz) +
        matrix(as.vector(t(dY))[dx] %*% jacobian(function(S) get_mat(L, S, "linear"), S), mx+mz, mx+mz) +
        matrix(jacobian(function(S) -get_mat(L, S, "mquad"), S), mx+mz, mx+mz) +
        matrix(jacobian(function(S) get_mat(L, S, "mdet"), S), mx+mz, mx+mz) +
        matrix(as.vector(t(dy))[dx] %*% jacobian(function(S) get_mat(L, S, "mlinear"), S), mx+mz, mx+mz) 

  dr2 <- 
        matrix(jacobian(function(y) -get_vec_marg(y, L, S, "mquad"), ym),mz,n) +
        matrix(as.vector(t(dy))[dx] %*% jacobian(function(y) get_vec_marg(y, L, S, "mlinear"), ym), mz,n)

  dy2 <-
        matrix(jacobian(function(y) -get_vec(y, L, S, "quad"), yc),mx,n) +
        matrix(as.vector(t(dY))[dx] %*% jacobian(function(y) get_vec(y, L, S, "linear"), yc), mx,n)

  if (raw)
    list("epee" = 
         list(sym(dL), sym(dS), dr, dy),
       "inv" =
         list(sym(dL2), sym(dS2), dr2, dy2),
       "dat" =
         list(ym, yc, typ)
       )
  else
    c("dL"=all(dplyr::near(sym(dL)[lower.tri(dL)], sym(dL2)[lower.tri(dL2)])),
      "dS"=all(dplyr::near(sym(dS)[lower.tri(dS)], sym(dS2)[lower.tri(dS2)])),
      "dr"=all(dplyr::near(dr, ifelse(is.na(dr2), 0, dr2))),
      "dy"=all(dplyr::near(dy, ifelse(is.na(dy2), 0, dy2))))
}

test_epee <- function(seed=1, n=5, mx=3, mz=2, tol=0, miss=0, scala=c(0,0))
{ 
  library(ape)
  set.seed(seed)
  tr    <- sim_polytomy(seed, n, tol)
  y     <- matrix(rnorm(n*(mx+mz)),mx+mz,n)*scala[1]
  miss  <- matrix(sample(c(TRUE,FALSE), n*(mx+mz), replace=TRUE, prob=c(miss,1-miss)),mx+mz,n)
  typ   <- sample(c(rep(1,mx),rep(0,mz)))
  y[typ==1,] <- sample(c(0,1), replace=TRUE, n*mx)
  #y[typ==0,] <- 0#debug
  y[miss] <- NA
  yc    <- y[typ==1,,drop=FALSE] 
  ym    <- y[typ==0,,drop=FALSE] 
  L     <- as.matrix(rWishart(1, mx+mz, diag(mx+mz))[,,1])
  S     <- as.matrix(rWishart(1, mx+mz, diag(mx+mz))[,,1])
  eta   <- matrix(rnorm((mx+mz)*n), mx+mz, n)*scala[2]
  eta[miss] <- NA
  em    <- eta[typ==0,,drop=FALSE]
  ec    <- eta[typ==1,,drop=FALSE]

  Cova <- Covariance$new(tr$edge-1, tr$edge.length, y, typ)
  Cova$parameters(L, S, eta)
  Cova$traverse_tree_marginal()

  Eps   <- Epee$new(Cova)
  Eps$run(Cova$Z, Cova)
  eps_sol <- Eps$Z

#  cova          <- getRefClass("covariance")
#  blah          <- cova$new(L, S, !as.logical(typ), y, tr)
#  eps           <- getRefClass("epee")
#  barf          <- eps$new(blah)
#  barf$run(ec, yc)
#  eps_sol_R <- barf$Z
  eps_sol_R <- NA

  dx    <- !is.na(as.vector(t(yc)))
  dz    <- !is.na(as.vector(t(ym)))
  K0    <- kronecker(L[typ==0,typ==0], vcv(tr)) + kronecker(S[typ==0,typ==0], diag(n))
  K1    <- kronecker(L[typ==1,typ==1], vcv(tr)) + kronecker(S[typ==1,typ==1], diag(n))
  K10   <- kronecker(L[typ==1,typ==0], vcv(tr)) + kronecker(S[typ==1,typ==0], diag(n))
  if(mz>0)
    K     <- K1[dx,dx] - K10[dx,dz] %*% solve(K0[dz,dz]) %*% t(K10[dx,dz])
  else
    K     <- K1[dx,dx]

  mn  <- as.vector(t(ec))[dx]
  if (mz > 0) mn <- mn + K10[dx,dz] %*% solve(K0[dz,dz]) %*% as.vector(t(ym-em))[dz]
  obs <- as.vector(t(yc))[dx]
  lo  <- ifelse(obs, 0, -Inf)
  hi  <- ifelse(obs, Inf, 0)
  int_sol <- log(mvtnorm::pmvnorm(lo, hi, as.vector(mn), sigma = K))

  return(list(epee=eps_sol, epee2=eps_sol_R, int=int_sol))
}

test_epee_deriv <- function(seed=1,n=5,mx=3,mz=2,miss=0,tol=0)
{
  library(ape);library(numDeriv)
  set.seed(seed)
  tr    <- sim_polytomy(seed, n, tol)
  y     <- matrix(rnorm(n*(mx+mz)),mx+mz,n)
  miss  <- matrix(sample(c(TRUE,FALSE), n*(mx+mz), replace=TRUE, prob=c(miss,1-miss)),mx+mz,n)
  typ   <- sample(c(rep(1,mx),rep(0,mz)))
  y[miss] <- NA
  y[typ==1,] <- sample(c(0,1), replace=TRUE, n*mx)
  y[typ==0,] <- 0#debug
  yc    <- y[typ==1,,drop=FALSE] 
  ym    <- y[typ==0,,drop=FALSE] 
  L     <- diag(mx+mz)#as.matrix(rWishart(1, mx+mz, diag(mx+mz))[,,1])
  S     <- diag(mx+mz)#as.matrix(rWishart(1, mx+mz, diag(mx+mz))[,,1])
  eta   <- matrix(0,mx+mz, n)#matrix(rnorm((mx+mz)*n), mx+mz, n)
  eta[miss] <- NA
  em    <- eta[typ==0,,drop=FALSE]
  ec    <- eta[typ==1,,drop=FALSE]

  Cova <- Covariance$new(tr$edge-1, tr$edge.length, y, typ)
  Cova$parameters(L, S, eta)
  Cova$traverse_tree_marginal()

  Eps   <- Epee$new(Cova)
  Eps$run(Cova$Z, Cova)
  dL    <- Eps$dL
  dS    <- Eps$dS
  dY    <- Eps$dm

  get_solution <- function(L, S, Z)
  {
    Cova2 <- Covariance$new(tr$edge-1, tr$edge.length, y, typ)
    Cova2$parameters(L, S, eta)
    Cova2$traverse_tree_marginal()
    Eps2   <- Epee$new(Cova2)
    Eps2$options(1e-4, 100, TRUE)
    Eps2$run(Z, Cova2)
    Eps2$Z[length(Eps2$Z)]
  }

  dL2 <- matrix(jacobian(function(L) get_solution(L, S, Cova$Z), L),mx+mz,mx+mz)
  dS2 <- matrix(jacobian(function(S) get_solution(L, S, Cova$Z), S),mx+mz,mx+mz)
  #dY2 <- matrix(jacobian(function(Z) get_solution(L, S, Z), Cova$Z),mx,n)

  sym <- function(X)
  {
    X = X + t(X)
    diag(X) = diag(X)/2
    X[upper.tri(X)] <- NA
    X
  }

  return(list(
              list(epee=sym(dL), num=sym(dL2)),
              list(epee=sym(dS), num=sym(dS2))
              #list(epee=dY, num=dY2)
              ))

}

test_likelihood <- function(seed=1, n=5, mx=3, mz=2, ncov=0, tol=0, missing=c(0,0),scala=c(1,1), profile=FALSE)
{
  library(ape)
  set.seed(seed)

  # constraints
  .random_constraints <- function(n)
  {
    out <- matrix(NA, n, n)
    out[lower.tri(out)] <- sample(c(TRUE,FALSE), replace=TRUE, n*(n-1)/2) 
    diag(out) <- FALSE
    out
  }
  constraints_L <- .random_constraints(mx+mz)
  constraints_S <- .random_constraints(mx+mz)
  zeros_L <- which(constraints_L & !is.na(constraints_L), arr.ind=TRUE)-1
  zeros_S <- which(constraints_S & !is.na(constraints_S), arr.ind=TRUE)-1
  ll <- rnorm(sum(!constraints_L & !is.na(constraints_L)))
  ss <- rnorm(sum(!constraints_S & !is.na(constraints_S)))
  cons <- Constraints$new(zeros_L, zeros_S, mx+mz)

  tr    <- sim_polytomy(seed, n, tol)
  typ   <- sample(c(rep(1,mx),rep(0,mz)))
  y     <- matrix(NA,mx+mz,n)
  y[typ==1,] <- sample(c(0,1), n*mx, replace=TRUE)
  y[typ==0,] <- rnorm(n*mz)*scala[1]
  miss       <- matrix(sample(c(TRUE,FALSE), n*(mx+mz), replace=TRUE, prob=c(missing[1],1-missing[1])),mx+mz,n)
  y[miss]    <- NA
  yc <- y[typ==1,]
  ym <- y[typ==0,]

  # design
  cov <- list()
  beta <- list()
  for(i in 1:(mx+mz))
  {
    ccc <- matrix(rnorm(n*ncov), n, ncov)
    missx <- matrix(sample(c(TRUE,FALSE), n*(ncov), replace=TRUE, 
                    prob=c(missing[2],1-missing[2])),ncov,n)
    ccc[missx] <- NA
    cov[[i]] <- cbind(1, ccc)
    beta[[i]] <- rnorm(ncov+1)*scala[2]
  }
  desn = Design$new(y, cov, typ, tr$edge-1, tr$edge.length)

  L <- cons$apply(ll, cons$Lambda)
  S <- cons$apply(ss, cons$Sigma)
  L <- L%*%t(L)
  S <- S%*%t(S)
  eta   <- desn$linear(beta)

  myreg <- Regression$new(desn, cons)
  if (profile)
  {
    myreg$profile(profile, ll, ss, beta)
  }
  timing <- system.time(ep_ans <- myreg$loglikelihood(ll, ss, beta))["elapsed"]#note that these all get changed :-/

  #debug
#  Cova <- Covariance$new(tr$edge-1, tr$edge.length, y, typ)
#  Cova$parameters(L, S, eta)
#  Cova$traverse_tree_marginal()
#  eps <- Epee$new(Cova)
#  eps$run(Cova$Z, Cova)
#  eps_ans2 <- eps$Z

  #integration answer
  int_ans <- NA
  em    <- eta[typ==0,]
  ec    <- eta[typ==1,]
  dx    <- !is.na(as.vector(t(yc)))
  dz    <- !is.na(as.vector(t(ym)))
  K0    <- kronecker(L[typ==0,typ==0], vcv(tr)) + kronecker(S[typ==0,typ==0], diag(n))
  K1    <- kronecker(L[typ==1,typ==1], vcv(tr)) + kronecker(S[typ==1,typ==1], diag(n))
  K10   <- kronecker(L[typ==1,typ==0], vcv(tr)) + kronecker(S[typ==1,typ==0], diag(n))
  if(mz>0)
    K     <- K1[dx,dx] - K10[dx,dz] %*% solve(K0[dz,dz]) %*% t(K10[dx,dz])
  else
    K     <- K1[dx,dx]

  mn  <- as.vector(t(ec))[dx]
  if (mz > 0) mn <- mn + K10[dx,dz] %*% solve(K0[dz,dz]) %*% as.vector(t(ym-em))[dz]
  obs <- as.vector(t(yc))[dx]
  lo  <- ifelse(obs, 0, -Inf)
  hi  <- ifelse(obs, Inf, 0)
  int_ans <- 0
  if (mx > 0)
    int_ans <- int_ans + log(mvtnorm::pmvnorm(lo, hi, as.vector(mn), sigma = K))
  if (mz > 0)
    int_ans <- int_ans + mvtnorm::dmvnorm(as.vector(t(ym))[dz], as.vector(t(em))[dz], sigma=K0[dz,dz], log=TRUE)

  cat("timing: ", timing, "\n")
  return(c(ep_ans, int_ans))


}

test_likelihood_deriv <- function(seed=1, n=5, mx=3, mz=2, ncov=0, tol=0, missing=c(0,0),scala=c(1,1))
{
  library(ape)
  set.seed(seed)

  # constraints
  .random_constraints <- function(n)
  {
    out <- matrix(NA, n, n)
    out[lower.tri(out, diag=TRUE)] <- sample(c(TRUE,FALSE), replace=TRUE, n*(n+1)/2) 
    out
  }
  constraints_L <- .random_constraints(mx+mz)
  constraints_S <- .random_constraints(mx+mz)
  zeros_L <- which(constraints_L & !is.na(constraints_L), arr.ind=TRUE)-1
  zeros_S <- which(constraints_S & !is.na(constraints_S), arr.ind=TRUE)-1
  ll <- rnorm(sum(!constraints_L & !is.na(constraints_L)))
  ss <- rnorm(sum(!constraints_S & !is.na(constraints_S)))
  cons <- Constraints$new(zeros_L, zeros_S, mx+mz)

  tr    <- sim_polytomy(seed, n, tol)
  typ   <- sample(c(rep(1,mx),rep(0,mz)))
  y     <- matrix(NA,mx+mz,n)
  y[typ==1,] <- sample(c(0,1), n*mx, replace=TRUE)
  y[typ==0,] <- rnorm(n*mz)*scala[1]
  miss       <- matrix(sample(c(TRUE,FALSE), n*(mx+mz), replace=TRUE, prob=c(missing[1],1-missing[1])),mx+mz,n)
  y[miss]    <- NA
  yc <- y[typ==1,]
  ym <- y[typ==0,]

  # design
  cov <- list()
  beta <- list()
  for(i in 1:(mx+mz))
  {
    ccc <- matrix(rnorm(n*ncov), n, ncov)
    missx <- matrix(sample(c(TRUE,FALSE), n*(ncov), replace=TRUE, 
                    prob=c(missing[2],1-missing[2])),ncov,n)
    ccc[missx] <- NA
    cov[[i]] <- cbind(1, ccc)
    beta[[i]] <- rnorm(ncov+1)*scala[2]
  }
  desn = Design$new(y, cov, typ, tr$edge-1, tr$edge.length)

  L <- cons$apply(ll, cons$Lambda)
  S <- cons$apply(ss, cons$Sigma)
  L <- L%*%t(L)
  S <- S%*%t(S)
  eta   <- desn$linear(beta)

  getLL <- function(desn, cons, ll, ss, beta)
  {
    myreg <- Regression$new(desn, cons)
    myreg$loglikelihood(ll, ss, beta)
  }

  getLLd <- function(desn, cons, ll, ss, beta)
  {
    myreg <- Regression$new(desn, cons)
    nuts <- myreg$loglikelihood(ll, ss, beta)
    return(list(myreg$dL, myreg$dS, myreg$db))
  }

  library(numDeriv)
  ep_ans <- getLLd(desn, cons, ll, ss, beta)
  ch_ans_ll <- jacobian(function(ll) getLL(desn, cons, ll, ss, beta), ll)
  ch_ans_ss <- jacobian(function(ss) getLL(desn, cons, ll, ss, beta), ss)
  ch_ans_be <- list()
  for (i in 1:length(beta))
    ch_ans_be[[i]] <- jacobian(function(be){beta[[i]] <- be; getLL(desn, cons, ll, ss, beta)},
                               beta[[i]])

  return(list(ep_ans, list(ch_ans_ll, ch_ans_ss, ch_ans_be), constraints_L, constraints_S))#ok but signs are off on be
}

sucks <- function(seed=1,n=10,l=0, b=0, nlm=FALSE)
{
  if (nlm)
    optimiz = "nlm"
  else
    optimiz = "Nelder-Mead"
  hi <- random_data(seed=seed,n=n, mz=1, mx=2, l=l, b=b)
  epee(list(cat1~1,cat2~1,cont1~1), data=hi$df, tree=hi$tree, optimizer=optimiz, 
       constraints_within=list(c("cat1","cat2","cont1")),
       constraints_between=list(c("cat1","cat2","cont1")), hessian=TRUE)
#  epee(list(cont1~1), data=hi$df, tree=hi$tree, optimizer=optimiz, hessian=TRUE)
       #constraints_within=list(c("cat1","cat2","cont1")),
       #constraints_between=list(c("cat1","cat2","cont1")), hessian=TRUE)
}
