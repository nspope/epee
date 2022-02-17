random_data <- function(seed=1, n=10, mx=3, mz=2, l=0, b=0)
{
  library(ape)
  set.seed(seed)
  tree=ape::rcoal(n)
  tree$edge.length=tree$edge.length/ape::vcv(tree)[1,1]#scale to height 1
  V <- cov2cor(ape::vcv(tree))
  df <- t(MASS::mvrnorm(mx+mz, rep(b, n), exp(l)^2*V + diag(n))) 
  df[,1:mx] <- ifelse(df[,1:mx] > 0, 1, 0)
  colnames(df) <- c(if(mx) paste0("cat", 1:mx) else NULL, if(mz) paste0("cont", 1:mz) else NULL)
  df <- as.data.frame(df)
  for(i in 1:mx)
    df[[i]] <- as.factor(df[[i]])
  list(df=as.data.frame(df), tree=tree)
}

parse_formulae <- function(f, data)
{
  name <- as.character(f)[2]
  if (!(name %in% colnames(data)))
    stop("one of variable names not found in data")
  type <- is.factor(data[[name]]) || is.character(data[[name]])
  if (type)
  {
    uv <- unique(data[[name]])
    if (length(uv[!is.na(uv)])>2)
      stop("currently can only use binary categorical variables")
  }
  attr(f, "type") <- type
  f <- list(f)
  names(f) = name
  f
}

likelihood2d <- function(data, tree, Sigma=diag(nrow(data)), Lambda=diag(nrow(data)), mu=rep(0,nrow(data)))
{
  # for comparison against hiscott2d package
  dim      <- nrow(data)
  ntips    <- length(tree$tip.label)
  tree     <- ape::reorder.phylo(tree, "postorder")

  traits   <- data
  design   <- rep(list(matrix(1,ntips,1)), dim)
  type     <- as.numeric(apply(data, 1, function(x) all(x < .Machine$double.eps || 1-x < .Machine$double.eps)))
  desn     <- Design$new(traits, design, type, tree$edge-1, tree$edge.length)

  cons     <- Constraints$new(matrix(NA,0,2), matrix(NA,0,2), dim) # no constraints

  llfac    <- Regression$new(desn, cons)

  # covariance must be parameterized as *lower triangular* log-cholesky factors
  Sigma    <- t(chol(Sigma))
  diag(Sigma) <- log(diag(Sigma))
  Lambda   <- t(chol(Lambda))
  diag(Lambda) <- log(diag(Lambda))

  sigmaf   <- Sigma[lower.tri(Sigma,diag=TRUE)]
  lambdaf  <- Lambda[lower.tri(Lambda,diag=TRUE)]

  llfac$loglikelihood(lambdaf, sigmaf, as.list(mu))
}

# note on identifiability: the following are equivalent ...
# mvtnorm::pmvnorm(c(0,0),c(Inf,Inf),c(1,1),sigma=matrix(c(1,0.2,0.2,1),2,2))
# k = 2
# mvtnorm::pmvnorm(c(0,0),c(Inf,Inf),c(1,1)*k,sigma=matrix(c(1,0.2,0.2,1)*k^2,2,2))
# mvtnorm::pmvnorm(c(0,0),c(Inf,Inf),c(1,k),sigma=matrix(c(1,0.2,0.2*k,k^2),2,2))

# so ... how to make this identifiable? basically, need the categorical traits done in terms of covariance matrix

# formula + data creates design matrix
epee <- function(formulas, data, tree, taxa = NULL, 
                 constraints_within=list(), constraints_between=list(), 
                 optimizer=c("nlm","Nelder-Mead"),
                 ...) # arguments passed to optimizer
{
  library(ape)
  formulas <- sapply(formulas, parse_formulae, data=data) # stop if categorical > 2 traits
  dim      <- length(formulas)
  tips     <- length(tree$tip.label)
  tree     <- ape::reorder.phylo(tree, "postorder")
  tree$edge.length <- tree$edge.length/ape::vcv(tree)[1,1]

  # check that taxon labels are valid
  if (is.null(taxa))
  {
    warning("taxon names not provided: assuming data is ordered as per tree$tip.labels")
    taxa   <- tree$tip.label
  } else if (length(taxa) != tips)
    stop("must have a row in data for each tip on tree")

  # reorder data to match tip ordering on tree
  taxa_map <- match(taxa, tree$tip.label)
  if (any(is.na(taxa_map)))
    stop("one of supplied taxon names is not found in tip labels of tree")
  data     <- data[taxa_map,]

  # reorder traits so as to put categorical traits first
  formulas <- formulas[order(sapply(formulas, attr, which="type"), decreasing=TRUE)]

  # assemble traits, design matrices
  traits   <- matrix(NA, dim, tips)
  type     <- rep(NA, dim)
  design   <- list()
  for (i in 1:dim)
  {
    type[i]     <- attr(formulas[[i]], "type")
    traits[i,]  <- as.numeric(data[[names(formulas)[i]]])
    if (type[i]) # categorical traits -> {0, 1}
      traits[i,] <- traits[i,] - 1
    design[[i]] <- model.matrix(formulas[[i]], data=data)
    if (type[i])
      constraints_within <- c(constraints_within, names(formulas)[i])
  }
  desn     <- Design$new(traits, design, type, tree$edge-1, tree$edge.length)

  # assemble constraints
  build_constraints <- function(co)
  {
    C <- matrix(NA, 0, 2)
    if (is.list(co))
      for (i in co) 
      {
        if (length(i)==1)
          i <- c(i, i)
        mn <- match(i, names(formulas))-1
        if (any(is.na(mn)))
          stop("variable name supplied in constraints does not match response variables")
        C <- rbind(C, t(combn(mn, 2)))
      }
    else
      stop("constraints must be lists")
    C 
  }
  cL       <- build_constraints(constraints_between)
  cS       <- build_constraints(constraints_within)
  cons     <- Constraints$new(cL, cS, dim)

  # loglikelihood object
  llfac    <- Regression$new(desn, cons)

  # starting values?
  # for now, using beta=0 and Lambda, Sigma = I
  nbeta    <- sapply(design, ncol)
  nLambda  <- dim*(dim+1)/2 - nrow(cons$Lambda)
  nSigma   <- dim*(dim+1)/2 - nrow(cons$Sigma)
  start    <- rep(0, sum(nbeta)+nLambda+nSigma)

  # optimize
  # nlm is slow, but easy to supply a gradient to ... 
  # switch in the future
  opti     <- function(p, return_grad=FALSE)
  {
    vp  <- nSigma + nLambda
    if (nSigma)
      Sin <- matrix(p[1:nSigma], nSigma, 1)
    else
      Sin <- matrix(NA, 0, 0)
    if (nLambda)
      Lin <- matrix(p[(nSigma + 1):vp], nLambda, 1)
    else
      Lin <- matrix(NA, 0, 0)

    beta <- list()
    for (i in 1:dim)
    {
      np        <- vp + sum(nbeta[1:i])
      beta[[i]] <- p[(np-nbeta[i]+1):np]
    }
    ll <- -llfac$loglikelihood(Lin,
                               Sin,
                               beta)
    #debug
#    K<-vcv(tree)*exp(Lin)^2 + diag(tips)
#    m<-rep(beta[[1]],tips)
#    lo<-ifelse(traits[1,], 0, -Inf)
#    hi<-ifelse(traits[1,], Inf, 0)
#    ll2<-log(mvtnorm::pmvnorm(lo,hi,m,sigma=K))
#    print(c(ll,ll2,Lin,beta[[1]]))
    gra <- c(llfac$dS, llfac$dL, unlist(llfac$db))
    attr(ll, "gradient") <- -gra
    if (return_grad)
      return(gra)
    else
      return(ll)
  }

  if (optimizer=="nlm")
    ans      <- nlm(opti, start, typsize = rep(0.25, sum(nbeta)+nLambda+nSigma), hessian=TRUE)
  else
  {
    ans           <- optim(start, opti, ...)
    ans$estimate  <- ans$par
    ans$minimum   <- ans$value
    ans$code      <- ans$convergence
  }

  # format output
  hessian   <- solve(ans$hessian)
  if (nSigma)
  {
    S   <- matrix(ans$estimate[1:nSigma], nSigma, 1)
    Sse <- matrix(diag(hessian)[1:nSigma], nSigma, 1)
  }
  else
  {
    S   <- matrix(NA, 0, 0)
    Sse <- matrix(NA, 0, 0)
  }
  if (nLambda)
  {
    L <- matrix(ans$estimate[(nSigma + 1):(nSigma+nLambda)], nLambda, 1)
    Lse <- matrix(diag(hessian)[(nSigma+1):(nSigma+nLambda)], nLambda, 1)
  }
  else
  {
    L <- matrix(NA, 0, 0)
    Lse <- matrix(NA, 0, 0)
  }
  SS <- cons$apply(S, cons$Sigma)
  LL <- cons$apply(L, cons$Lambda)
  be <- ans$estimate[(nSigma+nLambda+1):(nSigma+nLambda+sum(nbeta))]

  Sigma   <- SS %*% t(SS)
  Sest    <- cbind("Estimate"=c(S), "StdErr"=c(sqrt(Sse)))
  attr(Sigma, "unconstrained") <- Sest
  colnames(Sigma) <- rownames(Sigma) <- names(formulas)
  Lambda                             <- LL %*% t(LL)
  Lest    <- cbind("Estimate"=c(L), "StdErr"=c(sqrt(Lse)))
  attr(Lambda, "unconstrained")      <- Lest
  colnames(Lambda) <- rownames(Lambda) <- names(formulas)
  beta    <- cbind("Estimate"=be, 
                   "StdErr"=sqrt(diag(hessian)[(nSigma+nLambda+1):(nSigma+nLambda+sum(nbeta))]))
  rownames(beta) <- names(formulas)
  loglik <- -ans$minimum
  attr(loglik, "df") <- sum(nLambda+nSigma+sum(nbeta))
  list(
       Lambda=Lambda,
       Sigma=Sigma,
       beta=beta, 
       loglik=loglik,
       convergence=ans$code,
       vcov=round(hessian,4))
}


