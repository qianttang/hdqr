#' @export

coef.hdqr.gam <- function(object, s=NULL, 
    type=c("coefficients", "nonzero"), ...) {
  type = match.arg(type)
  b0 = t(as.matrix(object$b0))
  rownames(b0) = "(Intercept)"
  nbeta = rbind2(b0, object$beta)
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
predict.hdqr.gam <- function(object, newx, s=NULL, ...) {
  b0 = t(as.matrix(object$b0))
  rownames(b0) = "(Intercept)"
  nbeta = rbind2(b0, object$beta)
  ddd = as.integer(object$ddd)
  is.ns = object$is.ns
  np = dim(newx)
  nobs = as.integer(np[1])
  nvars = as.integer(np[2])
  model.matrix = array(NA, c(nobs, ddd * nvars))
  if(is.ns){
    for (j in seq.int(nvars)) {
      model.matrix[, j + nvars* (seq(ddd) - 1)] = ns(newx[, j], df=ddd, 
        knots=object$nknots[,j], Boundary.knots=object$bdknots[,j])
    }
  } else {
      for (j in seq.int(nvars)) {
        if( ddd !=3 ){
          model.matrix[, j + nvars* (seq(ddd) - 1)] = bs(newx[, j], df=ddd, 
          knots=object$nknots[,j], Boundary.knots=object$bdknots[,j])
        } else{
          model.matrix[, j + nvars* (seq(ddd) - 1)] = bs(newx[, j], df=ddd, 
          Boundary.knots=object$bdknots[,j])
        }
      }
  }

  newx = model.matrix 
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

