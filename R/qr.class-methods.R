#' @export

coef.hdqr.class <- function(object, s=NULL, t=NULL,
    type=c("coefficients", "nonzero"), ...) {
  type = match.arg(type)
  inx = which(abs(object$tau - t) < 1e-5)
  b0 = t(as.matrix(object$b0[inx,]))
  rownames(b0) = "(Intercept)"
  nbeta = rbind2(b0, object$beta[inx,,])
  if (!is.null(s)) {
    vnames = dimnames(nbeta)[[1]]
    dimnames(nbeta) = list(NULL, NULL)
    lambda = object$lambda
    lamlist = lambda.interp(lambda, s)
    nbeta = nbeta[,lamlist$left,drop=FALSE] %*% 
      Diagonal(x=lamlist$frac) +
      nbeta[,lamlist$right,drop=FALSE] %*% 
      Diagonal(x=1-lamlist$frac)
    dimnames(nbeta) = list(vnames, paste(seq(along=s)))
  }
  if (type == "coefficients") 
    return(nbeta)
  if (type == "nonzero") 
    return(nonzero(nbeta[-1, , drop=FALSE], bystep=TRUE))
} 
 

#' @export
predict.hdqr.class <- function(object, newx, s=NULL, t=NULL, ...) {
  inx = which(abs(object$tau - t) < 1e-5)
  b0 = t(as.matrix(object$b0[inx,]))
  rownames(b0) = "(Intercept)"
  nbeta = rbind2(b0, object$beta[inx,,])
  if (!is.null(s)) {
    vnames = dimnames(nbeta)[[1]]
    dimnames(nbeta) = list(NULL, NULL)
    lambda = object$lambda
    lamlist = lambda.interp(lambda, s)
    nbeta = nbeta[ , lamlist$left, drop=FALSE] %*% 
            Diagonal(x=lamlist$frac) +
            nbeta[ , lamlist$right, drop=FALSE] %*% 
            Diagonal(x=1-lamlist$frac)
    dimnames(nbeta) = list(vnames, paste(seq(along=s)))
  }
  nfit = as.matrix(as.matrix(cbind2(1, newx)) %*% nbeta)
  nfit
} 

#' @export
coef.cv.hdqr.class <- function(object,
                            s = c("lambda.1se", "lambda.min"), t,
                            ...) {
  if (is.numeric(s)) {
    lambda <- s
  } else if (is.character(s)) {
    s <- match.arg(s)
    lambda <- object[[s]]
  } else stop("Invalid form for s")
  coef(object$hdqr.class.fit, s = lambda, t=t, ...)
}

#' @export
predict.cv.hdqr.class <- function(object, newx,
                               s = c("lambda.1se", "lambda.min"), t, 
                               ...) {
  if (is.numeric(s)) {
    lambda <- s
  } else if (is.character(s)) {
    s <- match.arg(s)
    lambda <- object[[s]]
  } else stop("Invalid form for s")
    predict(object$hdqr.class.fit, newx, s = lambda, t=t, ...)
}