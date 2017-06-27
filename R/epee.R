ape2epee <- function(N, M, ind, seed=1)
{
	set.seed(seed)
  tree <- rcoal(N)
  traits <- matrix(1, nrow=N, ncol=M)
  tree <- reorder(tree, "postorder")
  y <- matrix(rnorm(N*M))

  if (nrow(traits) != length(tree$tip.label))
    stop("traits must be of same dimension as tree")

  lambda <- solve(rWishart(1, ncol(traits), diag(ncol(traits)))[,,1])
  sigma <- solve(rWishart(1, ncol(traits), diag(ncol(traits)))[,,1])
  
  if (max(ind) > ncol(traits) | min(ind) < 1)
    stop("can only condition on extant traits")

  do_test(lambda, sigma, t(traits), tree$edge-1, tree$edge.length, ape::node.depth.edgelength(tree), ind-1)

	V <- vcv(tree)
	W <- kronecker(V, lambda[ind,ind]) + kronecker(diag(length(tree$tip.label)), sigma[ind,ind])
	K <- solve(W)

  Wf <- kronecker(V, lambda) + kronecker(diag(length(tree$tip.label)), sigma)
  tr <- rep(1:M, N)
  tr_j <- tr %in% ind
  tr_i <- !(tr %in% ind)
  X <- apply(diag(M), 1, rep, N)
  X2 <- apply(diag(length(ind)), 1, rep, N)
  X3 <- apply(diag(M-length(ind)), 1, rep, N)
  Wc = Wf[tr_i,tr_i] - Wf[tr_i,tr_j]%*%solve(Wf[tr_j,tr_j])%*%Wf[tr_j,tr_i]

	return(list(V=V, W=W, K=K, Wc=Wc))
}

epee <- function(trait, type, design, tree)
{
  tree <- ape::reorder.phylo(tree, "postorder")
  E <- tree$edge - 1
  L <- tree$edge.length
  H <- ape::node.depth.edgelength(tree)
  Y <- trait
  dim <- nrow(trait)

  if (ncol(Y) != length(tree$tip.label))
    stop("must have data for all tips in the tree")

  if (dim >= length(tree$tip.label))
    stop("must have more tips than traits")

  if (length(design) != dim || 
      class(design) != "list" || 
      any(lapply(design, nrow) != length(tree$tip.label)))
    stop("must have a design matrix for each trait, with number of rows equalling the number of tips")

  design <- lapply(design, as.matrix)

  return( list(epee = BrownianMotion$new(Y, type, E, L, H, design), trait = trait, design = design, tree = tree, dim = dim) )
}

testcase_continuous <- function(N, K)
{
  library(ape); library(mvtnorm)
  tree <- reorder(rcoal(N), "postorder")
  traits <- matrix(rnorm(N*K),K,N)
  design <- rep(list(cbind(rep(1,N))),K)
  type <- rep(0,K)
  ep <- epee(traits, type, design, tree)
  epans <- ep$epee$loglikelihood(diag(K),diag(K),rep(list(c(1)),K))
  blah <- dmvnorm(as.vector(t(traits)), mean = rep(1,N*K), sigma = kronecker(diag(K),vcv(tree)) + diag(N*K), log=TRUE)
  return(c(epans, blah))
}

testcase_categorical <- function(seed, N, K, tol=1e-4, maxiter=100)
{
  set.seed(seed)
  library(ape); library(mvtnorm)
  tree <- reorder(rcoal(N), "postorder")
  traits <- matrix(rbinom(N*K, size=1, prob=0.5),K,N)
  design <- rep(list(cbind(rep(1,N))),K)
  type <- rep(1,K)
  ep <- epee(traits, type, design, tree)
  epans <- ep$epee$loglikelihood(diag(K),diag(K),rep(list(c(1)),K),tol,maxiter)
  lo <- ifelse(as.vector(t(traits)) == 0, -Inf, 0)
  hi <- ifelse(as.vector(t(traits)) == 0, 0, Inf)
  blah <- log(pmvnorm(lo, hi, mean = rep(1,N*K), sigma = kronecker(diag(K),vcv(tree)) + diag(N*K)))
  return(c(epans, blah))
}

testcase_mixed <- function(seed, N, K, i, tol=1e-4, maxiter=100)
{
  set.seed(seed)
  library(ape); library(mvtnorm)
  tree <- reorder(rcoal(N), "postorder")

  ncon <- length(i)
  ncat <- K - ncon

  # at this point redo traits & type
  categorical <- matrix(rbinom(ncat * N, size=1, prob=0.5), ncat, N)
  continuous <- matrix(rnorm(N*ncon), ncon, N)
  #continuous <- matrix(1, ncon, N)
  traits <- rbind(continuous, categorical)

  type <- c(rep(0,ncon), rep(1,ncat))

  design <- rep(list(cbind(rep(1,N))),K)
  ep <- epee(traits, type, design, tree)
  epans <- ep$epee$loglikelihood(0.2 + diag(K),diag(K),rep(list(c(1)),K),tol,maxiter)

  # approximation
  fM <- kronecker(0.2 + diag(K),vcv(tree)) + diag(N*K)
  tr <- rep(1:K, each=N)
  mM <- fM[tr %in% i, tr %in% i]
  cM <- fM[!(tr %in% i), !(tr %in% i)] - fM[!(tr %in% i), tr %in% i] %*% solve(mM) %*% fM[tr %in% i, !(tr %in% i)]
  cmean <- 1 + fM[!(tr %in% i), tr %in% i] %*% solve(mM) %*% (as.vector(t(continuous)) - 1)

  blah1 <- dmvnorm(as.vector(t(continuous)), mean = rep(1,N*ncon), sigma = mM, log=TRUE)
  lo <- ifelse(as.vector(t(categorical)) == 0, -Inf, 0)
  hi <- ifelse(as.vector(t(categorical)) == 0, 0, Inf)
  blah2 <- log(pmvnorm(lo, hi, mean = as.vector(cmean), sigma = cM))
  return(c(epans, blah1+blah2))
}
