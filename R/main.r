#' Generate a dataset from Kang & Schafer simulation
#'
#' @param n Sample size
#'
#' @return A data.frame with n rows
KSdatagen <- function(n) {
  kBeta <- c(210, 27.4, 13.7, 13.7, 13.7)
  kTheta <- c(-1, 0.5, -0.25, -0.1)
  Z <- matrix(rnorm(n * 4), nrow = n, ncol = 4)
  X <- matrix(nrow = n, ncol = 4)
  X[, 1] <- exp(Z[, 1] / 2)
  X[, 2] <- Z[, 2] / (1 + exp(Z[, 1])) + 10
  X[, 3] <- (Z[, 1] * Z[, 3] / 25 + 0.6) ^ 3
  X[, 4] <- (Z[, 2] + Z[, 4] + 20) ^ 2
  y <- rnorm(n, mean = cbind(1, Z) %*% kBeta)
  true.prop <- 1 / (1 + exp(-Z %*% kTheta))
  T <- rbinom(n, 1, true.prop)
  dat <- list(Y = y, T = T, Z = Z, X = X)
  return(dat)
}

.trim <- function(x) pmax(x, 1e-10)

.compute.quantile <- function(Q, w, q, r) {
  F <- function(y)sapply(y, function(x)mean(rowSums((Q <= x) * w)))
  inv <- function(qq){
    froot <- function(x) (F(x) - qq)
    uniroot(froot, r, extendInt = 'yes')$root
  }
  return(sapply(q, function(qq)inv(qq)))
}

#' Outcome distribution estimator (G-computation)
#'
#' @param y Vector with outcome values.
#' @param t Vector with binary treatment indicator.
#' @param Q Conditional outcome distribution estimate. This should come in the for of an n x p matrix, where each column represents a conditional quantile. (See Kang & Schafer example in tmle() function)
#' @param g Propensity score.
#' @param q Quantile to be computed (e.g., q = 0.5 for the median.)
#'
#' @return A point estimate

od <- function(y, t, Q, g, q){
  n <- length(y)
  w  <- matrix(1/dim(Q)[2], ncol = dim(Q)[2], nrow = n)
  chiq <- .compute.quantile(Q, w, q, range(y))
  return(chiq)

}

#' Targeted minimum loss based estimator (TMLE)
#'
#' @param y Vector with outcome values.
#' @param t Vector with binary treatment indicator.
#' @param Q Conditional outcome distribution estimate. This should come in the for of an n x p matrix, where each column represents a conditional quantile. (See Kang & Schafer example)
#' @param g Propensity score.
#' @param q Quantile to be computed (e.g., q = 0.5 for the median.)
#' @return A point estimate
#' @export
#' @examples
#' n.quant <- 500
#' formT <- T ~ X1+X2+X3+X4
#' formY <- Y ~ X1+X2+X3+X4
#' data <- KSdatagen(1000)
#' X <- data$X
#' Y <- data$Y
#' T <- data$T
#' fitT <- glm(formT, data = data.frame(T=T, X), family = binomial)
#' fitY1 <- lm(formY, data = data.frame(Y=Y, T=T, X), subset = T == 1)
#' fitY0 <- lm(formY, data = data.frame(Y=Y, T=T, X), subset = T == 0)
#' median1 <- predict(fitY1, newdata = data.frame(T=1, X))
#' median0 <- predict(fitY0, newdata = data.frame(T=0, X))
#' Q1 <- sapply(seq(1/n.quant, 1 - 1/n.quant, 1/n.quant),
#'             function(q)qnorm(q, mean = median1, sd = summary(fitY1)$sigma))
#' Q0 <- sapply(seq(1/n.quant, 1 - 1/n.quant, 1/n.quant),
#'             function(q)qnorm(q, mean = median0, sd = summary(fitY0)$sigma))
#' g1 <- trim(predict(fitT, type = 'response'))
#' ame <- tmle(Y, T, Q1, g1, 0.5) - tmle(Y, 1 - T, Q0, 1 - g1, 0.5)

tmle <- function(y, t, Q, g, q){

  n <- length(y)
  D  <- function(y, w, chiq){
    1 / g * ((y <= chiq) - rowSums((Q <= chiq) * w))
  }
  w  <- matrix(1/dim(Q)[2], ncol = dim(Q)[2], nrow = n)
  h <- t
  chiq <- .compute.quantile(Q, w, q, range(y))
  Do <- D(y, w, chiq)
  Dq <- D(Q, w, chiq)
  iter     <- 1
  crit     <- TRUE
  max.iter <- 20
  while(crit && iter <= max.iter){
    est.eq <- function(eps){
      out <- - mean(h * (Do - rowSums(Dq * exp(eps * Dq) * w) /
                           rowSums(exp(eps * Dq) * w)))
      return(out)
    }
    loglik <- function(eps){
      out <- - mean(h * (eps * Do - log(rowSums(exp(eps * Dq) * w))))
      return(out)
    }
    eps <- optim(par = 0, loglik, gr = est.eq, method = 'BFGS')$par
    w <- exp(eps * Dq) * w / rowSums(exp(eps * Dq) * w)
    chiq <- .compute.quantile(Q, w, q, range(y))
    Do <- D(y, w, chiq)
    Dq <- D(Q, w, chiq)
    iter <- iter + 1
    crit <- abs(eps)  > 1e-4 / n^0.6
  }
  return(chiq)
}

#' Firpo estimator (Firpo)
#'
#' @param y Vector with outcome values.
#' @param t Vector with binary treatment indicator.
#' @param Q Conditional outcome distribution estimate. This should come in the for of an n x p matrix, where each column represents a conditional quantile. (See Kang & Schafer example)
#' @param g Propensity score.
#' @param q Quantile to be computed (e.g., q = 0.5 for the median.)
#' @return A point estimate

firpo <- function(y, t, Q, g, q){
  library(quantreg)
  h <- t / g
  chiq <- coef(rq(y ~ 1, weights = h, tau = q))
  names(chiq) <- NULL
  return(chiq)
}

#' Inverse probability weighted estimator (IPW)
#'
#' @param y Vector with outcome values.
#' @param t Vector with binary treatment indicator.
#' @param Q Conditional outcome distribution estimate. This should come in the for of an n x p matrix, where each column represents a conditional quantile. (See Kang & Schafer example)
#' @param g Propensity score.
#' @param q Quantile to be computed (e.g., q = 0.5 for the median.)
#' @return A point estimate

ipw <- function(y, t, Q, g, q){
  n  <- length(y)
  w  <- matrix(1/dim(Q)[2], ncol = dim(Q)[2], nrow = n)
  h <- t / g
  D  <- Vectorize(function(chiq)  mean(h * (y <= chiq) - q))
  chiq <- uniroot(D, c(-1000, 1000), extendInt = 'yes')$root
  return(chiq)
}

#' Augmented inverse probability weighted estimator (AIPW)
#'
#' @param y Vector with outcome values.
#' @param t Vector with binary treatment indicator.
#' @param Q Conditional outcome distribution estimate. This should come in the for of an n x p matrix, where each column represents a conditional quantile. (See Kang & Schafer example)
#' @param g Propensity score.
#' @param q Quantile to be computed (e.g., q = 0.5 for the median.)
#' @return A point estimate

aipw <- function(y, t, Q, g, q){
  n <- length(y)
  w <- matrix(1/dim(Q)[2], ncol = dim(Q)[2], nrow = n)
  h <- t / g
  D  <- Vectorize(function(chiq){
    mean(h * ((y <= chiq) - rowSums((Q <= chiq) * w))) +
      mean(rowSums((Q <= chiq) * w) - q)
  })
  chiq <- uniroot(D, c(-1000, 1000), extendInt = 'yes')$root
  return(chiq)
}

