stan = function(train, validation=NULL, test=NULL) {
  # standardize data set
  # Args:
  #   train:      original training set
  #   validation: original validation set
  #   test:       original test set
  # Returns:
  #   train:      standardized training set
  #   validation: standardized validation set
  #   test:       standardized test set
  train_m     = colMeans(train$X)
  std_train_x = t(apply(train$X, 1, function(x) x - train_m))
  train_sd    = apply(std_train_x, 2,
                function(x) sqrt(x %*% x / length(x)))
  train_sd[train_sd==0] = 1
  train$X = t(apply(std_train_x, 1, function(x) x / train_sd))
  if (!is.null(validation)) validation$X =
    scale(validation$X, center=train_m, scale=train_sd)
  if (!is.null(test)) test$X =
    scale(test$X, center=train_m, scale=train_sd)
  rm.att = function(x) {
    attributes(x) = attributes(x)[c(1,2)]
    x
  }
  train$X = rm.att(train$X)
  validation$X = rm.att(validation$X)
  test$X = rm.att(test$X)
  # returns:
  list(train=train, validation=validation, test=test)
}

Data1Gen <- function(n, d, seed) {
  set.seed(seed)
  x = matrix(rnorm(n * d), n, d)
  beta0 = rnorm(d)
  y = c(x %*% beta0) + rnorm(n, sd = 0.1)
  return(list(x=x,y=y))
}

  # Sigma <- matrix(rho, p, p)
  # diag(Sigma) <- 1
  # eig <- eigen(Sigma)
  # Sigma.sq <- eig$vectors %*% tcrossprod(diag(sqrt(eig$values)),
  #                             eig$vectors)

# Data3Gen <- function(i, num, dim, rho=0.5, Sigma.sq){
#   set.seed(i)
#   n = num
#   p = dim
#   X <- matrix(rnorm(n*p, 0, 1), n, p) %*% Sigma.sq
#   beta <- rep(0, p)
#   beta_sub <- c(2, 0, 1.5, 0, 0.8, 0, 0, 1, 0, 1.75, 0, 0, 0.75, 0,0,0.3)
#   beta[1:16] <- beta_sub
#   # E <- rnorm(n, 0 , 2)
#   # E <- 0.9 * rnorm(n, 0 , 1) + 0.1 * rnorm(n, 0 , 25)
#   # E <- rnorm(n, 0, runif(1, 1, 5))
#   # E <- rlaplace(n, location = 0, scale = 1)
#   # E <- sqrt(2) * rt(n, 4)
#   E <- rcauchy(n)
#   Y = X %*% beta + E
#   Y <- as.matrix(Y,n,1)
#   return(list(x=X,y=Y))
# }

Data3Gen <- function(i, num, dim, rho=0.1, Sigma.sq){
  set.seed(i)
  n = num
  p = dim
  # Sigma <- matrix(rho, p, p)
  # diag(Sigma) <- 1
  # eig <- eigen(Sigma)
  # Sigma.sq <- eig$vectors %*% tcrossprod(diag(sqrt(eig$values)),
  #                             eig$vectors)
  X <- matrix(rnorm(n*p, 0, 1), n, p) %*% Sigma.sq
  cvec <- matrix(c(1:p), p, 1)
  Beta <- apply(cvec, 2, function(x) (-1)^x*exp(-(2*x-1)/20))
  E <- rnorm(n, 0, 1)
  Y = X %*% Beta + 3*E
  Y <- as.matrix(Y,n,1)
  return(list(x=X,y=Y))
}

Data2Gen <- function(i, num, dim, rho=0.5, Sigma.sq){
  set.seed(i)
  n = num
  p1 = dim
  X <- matrix(rnorm(n*7, 0, 1), n, 7) %*% Sigma.sq
  X <- cbind(X, matrix(rnorm(n*(p1-7), 0, 1), n, (p1-7)))
  beta <- rep(0, p1)
  beta_sub <- c(2, 1.5, 0.8, 1, 1.75, 0.75, 0.3)
  beta[1:7] <- beta_sub
  # E <- rnorm(n, 0 , 2)
  # E <- 0.9 * rnorm(n, 0 , 1) + 0.1 * rnorm(n, 0 , 25)
  # E <- rnorm(n, 0, runif(1, 1, 5))
  # E <- rlaplace(n, location = 0, scale = 1)
  # E <- sqrt(2) * rt(n, 4)
  E <- rcauchy(n)
  Y = X %*% beta + E
  Y <- as.matrix(Y,n,1)
  return(list(x=X,y=Y))
}


objfn_qr = function(beta, b0=0, x, y, tau, lam1, lam2)
{
    pred = b0 + x %*% beta
    loss = check_loss(y-pred, tau)
    reg = lam1 * sum(abs(beta)) + lam2 * sum(beta^2) / 2 
    mean(loss) + reg
}

huber_loss = function(x, delta)
{
    xabs = abs(x)
    ifelse(xabs <= delta, xabs^2 / 2, delta * (xabs - delta / 2))
}

objfn_huber = function(beta, b0, x, y, delta, lam1, lam2)
{
    pred = x %*% beta + b0
    loss = huber_loss(y - pred, delta)
    reg = lam1 * sum(abs(beta)) + lam2 * sum(beta^2) / 2
    mean(loss) + reg
}

cvx_qr <- function(beta, x, y, lam1, lam2, tau) {
  # p = dim(beta)[1]
  pred <- x %*% beta[-1] + beta[1]
  loss1 <- mean(check_loss(y - pred, tau=tau))
  loss2 <- lam2/2 * sum(beta[-1]^2) + lam1 * p_norm(beta[-1], 1)
  loss1 + loss2
}

