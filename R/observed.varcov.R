#' @importFrom stats D
#'
#' @name observed.varcov
#' @aliases observed.varcov
#'
#' @title Observed Fisher Information
#'
#' @description \code{observed.varcov} calculates the inverse of the observed Fisher Information. Analytical second-order partial log-density derivatives are used in the calculations.
#'
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#'
#' @param logdensity An expression with the log of the probability density function.
#' @param X A numeric vector with the observations.
#' @param parms A character vector with the parameter name(s) specified in the logdensity expression.
#' @param mle A numeric vector with the parameter estimate(s).
#'
#' @return \code{observed.varcov} returns a list with two components (i) \bold{mle}: the inputted maximum likelihood estimate(s) and (ii) \bold{varcov}: the observed variance-covariance evaluated at the inputted mle argument.
#'
#' @return If the observed information is singular an error message is returned.
#
#' @seealso \code{\link[stats]{deriv}}, \code{\link[stats]{D}}, \code{\link[mle.tools]{expected.varcov}}.
#'
#' @details The second-order partial log-density derivatives are calculated via \code{D} function.
#'
#' @examples
#' {library(mle.tools); library(fitdistrplus); set.seed(1)};
#'
#' ##Normal distribution
#' lpdf <- quote(-log(sigma) - 0.5 / sigma ^ 2 * (x - mu) ^ 2)
#'
#' x <- rnorm(n = 100, mean = 0.0, sd = 1.0)
#'
#' observed.varcov(logdensity = lpdf, X = x, parms = c("mu", "sigma"),
#'  mle = c(mean(x), sd(x)))
#'
#' ################################################################################
#'
#' ## Weibull distribution
#' lpdf <- quote(log(shape) - shape * log(scale) + shape * log(x) - (x / scale) ^ shape)
#'
#' x <- rweibull(n = 100, shape = 1.5, scale = 2.0)
#'
#' fit <- fitdist(data = x, distr = 'weibull')
#' fit$vcov
#'
#' observed.varcov(logdensity = lpdf, X = x, parms = c("shape", "scale"), mle = fit$estimate)
#'
#' ################################################################################
#'
#' ## Exponetial distribution
#' lpdf <- quote(log(rate) - rate * x)
#'
#' x <- rexp(n = 100, rate = 0.5)
#'
#' fit <- fitdist(data = x, distr = 'exp')
#' fit$vcov
#'
#' observed.varcov(logdensity = lpdf, X = x, parms = c("rate"), mle = fit$estimate)
#'
#' ################################################################################
#'
#' ## Gamma distribution
#' lpdf <- quote(-shape * log(scale) - lgamma(shape) + shape * log(x) -
#'  x / scale)
#'
#' x <- rgamma(n = 100, shape = 1.5, scale = 2.0)
#'
#' fit <- fitdist(data = x, distr = 'gamma', start = list(shape = 1.5, scale =  2.0))
#' fit$vcov
#'
#' observed.varcov(logdensity = lpdf, X = x, parms = c("shape", "scale"), mle = fit$estimate)
#'
#' ################################################################################
#'
#' ## Beta distribution
#' lpdf <- quote(lgamma(shape1 + shape2) - lgamma(shape1) - lgamma(shape2) +
#'   shape1 * log(x) + shape2 * log(1 - x))
#'
#' x <- rbeta(n = 100, shape1 = 2.0, shape2 = 2.0)
#'
#' fit <- fitdist(data = x, distr = 'beta', start = list(shape1 = 2.0, shape2 =  2.0))
#' fit$vcov
#'
#' observed.varcov(logdensity = lpdf, X = x, parms = c("shape1", "shape2"), mle = fit$estimate)
#'
#' @rdname observed.varcov
#' @export
observed.varcov <- function(logdensity, X, parms, mle)
{
  {p <- length(parms); l <- length(mle); n <- length(X)};
  if(p != l) stop("The arguments 'parms' and 'mle' must be have the same size")

  obs <- matrix(NA_real_, ncol = p, nrow = p)
  {colnames(obs) <- parms; rownames(obs) <- parms; names(mle) <- parms};

  {H <- function(x){}; ll <- H};

  for(i in 1:p)
  {
    assign(parms[i], mle[i])
  }

  #body(ll) <- bquote(.(logdensity))
  #{n <- length(X); a <- -2 * sum(ll(x = X), na.rm = T); b <- a + 2 * p; d <- b + 2 * p * (p + 1) / (n - p - 1); e <- a + p * log(n); f = a + 2 * log(log(n)) * p}
  #{fitstats <- data.frame(Neg2LogLike = a, AIC = b, AICc = d, BIC = e, HQC = f); rownames(fitstats) <- ''}

  first <- sapply(1:p, function(i) D(logdensity, parms[i]))

  for(i in 1:p)
  {
    for(j in i:p)
    {
      second <- D(first[[i]], parms[j])
      if(length(grep('x', second)) == 0)  second <- bquote(.(n) * .(second))
      body(H) <- bquote(.(second))
      obs[i, j] <- sum(H(x = X), na.rm = T)
      if(j > i) obs[j, i] <- obs[i, j]
    }
  }
  if(any(eigen(-obs)$values < 0)) stop("The final Hessian matrix has at least one negative eigenvalue")
  else return(list(mle = mle, varcov = solve(-obs)))  ##likelihood.based.statistics= fitstats))
}

