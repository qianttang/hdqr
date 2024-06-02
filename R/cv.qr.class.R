#' @export
#'
cv.tau <- function(x, y, tau, lambda = NULL, nfolds=5L,
                       foldid, delta=.125, lam2=0.01, q, ...){
  ############################################################################
  ## data setup
  y <- drop(y)
  x <- as.matrix(x)
  x.row <- as.integer(NROW(x))
  ntau <- length(tau)
  if (length(y) != x.row)
    stop("x and y have different number of observations.")
  hdqr.object = hdqr.class(x, y, lambda=lambda, tau=tau, ...)
  ###Now fit the nfold models and store them
  if (missing(foldid)) {
    foldid <- sample(rep(seq(nfolds), length = x.row))
  } else nfolds <- max(foldid)
  if (nfolds < 3L)
    stop("nfolds must be at least 3; nfolds = 5 recommended")
  lambda <- sort(lambda, decreasing=TRUE)
  outlist <- as.list(seq(nfolds))
  for (i in seq(nfolds)) {
    which <- foldid == i
    outlist[[i]] <- hdqr.class(x=x[!which, , drop=FALSE], y=y[!which],
      lambda=lambda, tau=tau, ...)
    # if (outlist[[i]]$jerr != 0)
    #   stop(paste("Error occurs when fitting the", i, "th folder."))
  }
  cvmat <- f1.tau.path(outlist, x, y, tau, ntau, lambda, foldid, x.row, delta, lam2, q, ...)
  ## wrap up output
  inxmax <- which(cvmat==max(cvmat, na.rm=TRUE), arr.ind=TRUE)
  lambda.max <- lambda[inxmax[1,2]]
  tau.max <- tau[inxmax[1,1]]
  out <- list(lambda=lambda, tau=tau, cvmat=cvmat, lambda.max=lambda.max, 
    tau.max=tau.max, hdqr.class.fit=hdqr.object)
  obj <- c(out)
  class(obj) <- "cv.hdqr.class"
  obj
}

f1.tau.path <- function(outlist, x, y, tau, ntau, lambda, foldid, x.row,
                       delta, lam2, q, ...){
  nfolds <- max(foldid)
  cvmat <- matrix(NA, ntau, length(lambda))
  nlams <- double(nfolds)
  for (t in 1:ntau) {
    predmat <- matrix(NA, x.row, length(lambda))
    for (i in seq(nfolds)) {
      whichfold <- foldid == i
      fitobj <- outlist[[i]]
      preds <- predict(fitobj, x[whichfold, , drop = FALSE], t=tau[t])
      nlami <- length(lambda)
      # nlami <- length(fitobj$lambda)
      predmat[whichfold, seq(nlami)] <- preds
      nlams[i] <- nlami
    }
    cvmat[t, ] <- apply(predmat, 2, f1, y=y, q=q)
  }
  cvmat
}

f1 <- function(pred, y, q){
  pred <- matrix(pred, length(pred), 1)
  # compared with q, if val >= q, 1, else -1
  pred <- apply(pred, 1, phi_x, q=q)
  res <- confusionMatrix(as.factor(pred),as.factor(y))
  f1 = ifelse(is.na(res$byClass[7]), 0, res$byClass[7])
  return(f1)
}

phi_x <- function(z,qval){
  ifelse(z>=qval, 1, 0)
}