#' Solve the linear quantile regression. The solution path is computed
#' at a grid of values of tuning parameter \code{lambda}.
#'
#' @param x Matrix of predictors, of dimension (nobs * nvars);
#'   each row is an observation.
#' @param y Response variable. The length is \eqn{n}.
#' @param tau The quantile level \eqn{\tau}. The value must be
#'   in (0,1). Default is 0.5.
#' @param nlambda The number of \code{lambda} values (default is 100).
#' @param lambda.factor The factor for getting the minimal value
#'   in the \code{lambda} sequence, where
#'   \code{min(lambda)} = \code{lambda.factor} * \code{max(lambda)}
#'   and \code{max(lambda)} is the smallest value of \code{lambda}
#'   for which all coefficients (except the intercept when it is present)
#'   are penalized to zero. The default depends on the relationship
#'   between \eqn{n} (the number of rows in the design matrix) and
#'   \eqn{p} (the number of predictors). If \eqn{n < p}, it defaults to
#'   \code{0.05}. If \eqn{n > p}, the default is \code{0.001},
#'   closer to zero.  A very small value of \code{lambda.factor} will
#'   lead to a saturated fit. The argument takes no effect if there is a
#'   user-supplied \code{lambda} sequence.
#' @param lambda A user-supplied \code{lambda} sequence. Typically,
#'   by leaving this option unspecified, users can have the program
#'   compute its own \code{lambda} sequence based on \code{nlambda}
#'   and \code{lambda.factor}. It is better to supply, if necessary,
#'   a decreasing sequence of \code{lambda} values than a single
#'   (small) value. The program will ensure that the user-supplied
#'   \code{lambda} sequence is sorted in decreasing order before
#'   fitting the model to take advanage of the warm-start technique.
#' @param lam2 Regularization parameter \code{lambda2} for the
#'   quadratic penalty of the coefficients. Unlike \code{lambda},
#'   only one value of \code{lambda2} is used for each fitting process.
#' @param exclude Indices of variables to be excluded from the model.
#'   Default is none. Equivalent to an infinite penalty factor.
#' @param dfmax The maximum number of variables allowed in the model.
#'   Useful for very large \eqn{p} when a partial path is desired.
#'   Default is \eqn{p+1}.
#' @param pmax The maximum number of coefficients allowed ever
#'   to be nonzero along the solution path. For example, once
#'   \eqn{\beta} enters the model, no matter how many times it
#'   exits or re-enters the model through the path, it will be
#'   counted only once. Default is \code{min(dfmax*1.2, p)}.
#' @param standardize Logical flag for variable standardization,
#'   prior to fitting the model sequence. The coefficients are
#'   always returned to the original scale. Default is \code{TRUE}.
#' @param hval The smoothing index for \code{method='huber'}. Default is 0.125^2.
#' @param eps Stopping criterion.
#' @param maxit Maximum number of iterates.
#' @param sigma Penalty parameter appearing in the quadratic term
#'   of the augmented Lagrangian function. Must be positive.
#' @param is_exact Exact or approximated solutions. Default is \code{FALSE}.
#'
#' @details
#' Note that the objective function in the penalized quantile
#' regression is
#'   \deqn{1'\rho_{\tau}(y-X\beta-b_0))/N + \lambda_1\cdot|pf_1\circ\beta|_1 +
#'     0.5*\lambda_2\cdot|\sqrt{pf_2}\circ\beta|^2,}{
#'     1'\rho[\tau](y-X\beta))/N + \lambda[1]*|pf1•\beta|[1]
#'     + 0.5*\lambda_{2}*|\sqrtpf2•\beta|^2,}
#'   where \eqn{\rho_{\tau}}{\rho_{\tau}} the quantile or check loss
#'   and the penalty is a combination of weighted L1 and L2 terms and
#'   \eqn{\circ}{•} denotes the Hadmamard product.
#'
#' For faster computation, if the algorithm is not converging or
#' running slow, consider increasing \code{eps}, increasing
#' \code{sigma}, decreasing \code{nlambda}, or increasing
#' \code{lambda.factor} before increasing \code{maxit}.
#'
#' @return
#' An object with S3 class \code{hdqr} consisting of
#'   \item{call}{the call that produced this object}
#'   \item{b0}{intercept sequence of length \code{length(lambda)}}
#'   \item{beta}{a \code{p*length(lambda)} matrix of coefficients,
#'               stored as a sparse matrix (\code{dgCMatrix} class,
#'               the standard class for sparse numeric matrices in
#'               the \code{Matrix} package.). To convert it into
#'               normal type matrix, use \code{as.matrix()}.}
#'   \item{lambda}{the actual sequence of \code{lambda} values used}
#'   \item{df}{the number of nonzero coefficients for each value
#'             of \code{lambda}.}
#'   \item{npasses}{the number of iterations for every lambda value}
#'   \item{jerr}{error flag, for warnings and errors, 0 if no error.}
#'
#' @keywords quantile regression
#' @useDynLib hdqr, .registration=TRUE
#' @export

hdqr.class <- function(x, y, tau, nlambda=100, lambda.factor=ifelse(nobs < nvars, 0.01, 1e-04), 
    lambda=NULL, lam2=0, hval=.125, pf=rep(1, nvars), pf2=rep(1, nvars), 
    exclude, dfmax=nvars + 1, pmax=min(dfmax * 1.2, nvars), standardize=TRUE, 
    eps=1e-08, maxit=1e+06, sigma=0.05) {
  ####################################################################
  this.call = match.call()
  y = drop(y)
  x = as.matrix(x)
  np = dim(x)
  nobs = as.integer(np[1])
  nvars = as.integer(np[2])
  ntau = as.integer(length(tau))
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
  beta <- array(0, c(ntau, pmax, nlam))
  ####################################################################
  fit = .Fortran("qr_class", alpha, lam2, hval, nobs, nvars, 
    as.double(x), as.double(y), as.double(tau), ntau, jd, pfncol, pf, pf2, dfmax, 
    pmax, nlam, flmin, ulam, eps, isd, maxit, 
    nalam=integer(1), b0=double(ntau * nlam), 
    beta=as.double(beta), ibeta=integer(pmax), 
    nbeta=integer(ntau * nlam), alam=double(ntau * nlam), npass=integer(ntau * nlam), 
    jerr=integer(1), sigma=as.double(sigma), PACKAGE="hdqr")
  # outlist = getoutput(fit, maxit, pmax, nvars, vnames) ! rewrite getoutput
  # fit = c(outlist, list(npasses = fit$npass, jerr = fit$jerr))
  # fit$call = this.call
  # ####################################################################
  # class(fit) = c("hdqr.class")
  beta = array(fit$beta[seq(ntau * pmax * nlam)], c(ntau, pmax, nlam))
  b0 = matrix(fit$b0[seq(ntau * nlam)], ntau, nlam)
  fit=list(beta=beta, b0=b0, tau=tau)
  # if (is.null(lambda)) 
  fit$lambda = ulam
  class(fit) = c("hdqr.class")
  fit

} 