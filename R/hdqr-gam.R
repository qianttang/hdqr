#' @export

hdqr.gam <- function(x, y, tau, nlambda=100, lambda.factor=ifelse(nobs < nvars, 0.01, 1e-04), 
    lambda=NULL, lam2=0, hval=.125, pf=rep(1, nvars), pf2=rep(1, nvars), 
    exclude, dfmax=nvars + 1, pmax=min(dfmax * 1.2, nvars), standardize=TRUE, 
    eps=1e-08, maxit=1e+06, sigma=0.05, ddd=5L, is.ns=TRUE) {
  ####################################################################
  this.call = match.call()
  y = drop(y)
  x = as.matrix(x)
  np = dim(x)
  nobs = as.integer(np[1])
  nvars = as.integer(np[2])
  model.matrix = array(NA, c(nobs, ddd * nvars))
  nknots = matrix(0, ddd - 1, nvars)
  bdknots = matrix(0, 2, nvars)
  if(!is.ns) nknots = matrix(0, ddd - 3, nvars)
  for (j in 1:nvars) {
    if(is.ns) {
      tmp = ns(x[, j], df=ddd)  # natural spline
      nknots[, j] = attr(tmp, 'knots')
    } else {
      tmp = bs(x[, j], df=ddd)  # B-spline
      if(ddd != 3) nknots[, j] = attr(tmp, 'knots')
    }
    model.matrix[, j + nvars * (seq(ddd) - 1)] = tmp
    bdknots[, j] = attr(tmp, 'Boundary.knots')
  }
  x = model.matrix 
  np = dim(x)
  nobs = as.integer(np[1])
  nvars = as.integer(np[2])
  vnames = colnames(x)
  if (is.null(vnames)) 
    vnames = paste0("V", seq(nvars))
  if (length(y) != nobs) 
    stop("x and y have different number of observations")
  ####################################################################
  #parameter setup
  alpha = NULL # user can change this to tune elastic-net based on alpha.  
  if (!is.null(alpha)) {
    alpha = as.double(alpha)
    if (alpha <= 0 || alpha > 1) 
      stop("alpha: 0 < alpha <= 1")
    if (!is.null(lam2)) 
      warning("lam2 has no effect.")
    lam2 = -1.0
  } else {
    if (!is.null(lam2)) {
      if (lam2 < 0) stop("lam2 is non-negative.")
      alpha = -1.0
    } else {
      alpha = 1.0 # default lasso
    }
  }
  maxit = as.integer(maxit)
  isd = as.integer(standardize)
  eps = as.double(eps)
  dfmax = as.integer(dfmax)
  pmax = as.integer(pmax)
  if (!missing(exclude)) {
    jd = match(exclude, seq(nvars), 0)
    if (!all(jd > 0)) 
      stop("Some excluded variables are out of range.")
    jd = as.integer(c(length(jd), jd))
  } else jd = as.integer(0)
  ####################################################################
  #lambda setup
  nlam = as.integer(nlambda)
  if (is.null(lambda)) {
    if (lambda.factor >= 1) 
      stop("lambda.factor should be less than 1.")
    flmin = as.double(lambda.factor)
    ulam = double(1)
  } else {
    #flmin=1 if user define lambda
    flmin = as.double(1)
    if (any(lambda < 0)) 
      stop("The values of lambda should be non-negative.")
    ulam = as.double(rev(sort(lambda)))
    nlam = as.integer(length(lambda))
  }
  pfncol = NCOL(pf)
  pf = matrix(as.double(pf), ncol=pfncol)
  if (NROW(pf) != nvars) {
    stop("The size of L1 penalty factor must be the same with the number of input variables.")
  } else {
    if (pfncol != nlambda & NCOL(pf) != 1)
    stop("pf should be a matrix with 1 or nlambda columns.")
  }
  if (length(pf2) != nvars) 
    stop("The size of L2 penalty factor must be the same with the number of input variables.")
  pf2 = as.double(pf2)

  ####################################################################
  fit = .Fortran("lqr_hd", alpha, lam2, hval, nobs, nvars, 
    as.double(x), as.double(y), as.double(tau), jd, pfncol, pf, pf2, dfmax, 
    pmax, nlam, flmin, ulam, eps, isd, maxit, 
    nalam=integer(1), b0=double(nlam), 
    beta=double(pmax * nlam), ibeta=integer(pmax), 
    nbeta=integer(nlam), alam=double(nlam), npass=integer(nlam), jerr=integer(1),
    sigma=as.double(sigma), PACKAGE="hdqr")
  outlist = getoutput(fit, maxit, pmax, nvars, vnames)
  fit = c(outlist, list(npasses = fit$npass, jerr = fit$jerr))
  if (is.null(lambda)) 
    fit$lambda = lamfix(fit$lambda)
  fit$call = this.call
  fit$nknots  = nknots
  fit$bdknots = bdknots
  fit$ddd = ddd
  fit$is.ns = is.ns
  ####################################################################
  class(fit) = c("hdqr.gam")
  fit
} 

