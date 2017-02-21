#' @importFrom stats D
#' @importFrom stats integrate
#'
#' @name coxsnell.bc
#' @aliases coxsnell.bc
#'
#' @title Bias-Corrected Maximum Likelihood Estimate(s)
#'
#' @description \code{coxsnell.bc} calculates the bias-corrected maximum likelihood estimate(s) using the bias formula introduced by Cox and Snell (1968).
#'
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#'
#' @param density An expression with the probability density function.
#' @param logdensity An expression with the logarithm of the probability density function.
#' @param n A numeric scalar with the sample size.
#' @param parms A character vector with the parameter name(s) specified in the density and logdensity expressions.
#' @param mle A numeric vector with the parameter estimate(s).
#' @param lower The lower integration limit (lower = ``-Inf'' is the default).
#' @param upper The upper integration limit (upper = ``Inf'' is the default).
#' @param ... Additional arguments passed to \code{integrate} function.
#'
#' @return \code{coxsnell.bc} returns a list with five components (i) \bold{mle}: the inputted maximum likelihood estimate(s), (ii) \bold{varcov}: the expected variance-covariance evaluated at the inputted mle argument, (iii) \bold{mle.bc}: the bias-corrected maximum likelihood estimate(s), (iv) \bold{varcov.bc}: the expected variance-covariance evaluated at the bias-corrected maximum likelihood estimate(s) and (v) \bold{bias}: the bias estimate(s).
#'
#' @return If the numerical integration fails and/or the expected information is singular an error message is returned.
#
#' @seealso \code{\link[stats]{deriv}}, \code{\link[stats]{D}}, \code{\link[mle.tools]{expected.varcov}}, \code{\link[stats]{integrate}}, \code{\link[mle.tools]{observed.varcov}}.
#'
#' @details The first, second and third-order partial log-density derivatives are analytically calculated via \code{D} function. The expected values of the partial log-density derivatives are calculated via \code{integrate} function.
#'
#' @examples
#' {library(mle.tools); library(fitdistrplus); set.seed(1)};
#'
#' ## Normal distribution
#' pdf <- quote(1 / (sqrt(2 * pi) * sigma) * exp(-0.5 / sigma ^ 2 * (x - mu) ^ 2))
#' lpdf <- quote(- log(sigma) - 0.5 / sigma ^ 2 * (x - mu) ^ 2)
#'
#' x <- rnorm(n = 100, mean = 0.0, sd = 1.0)
#' {mu.hat <- mean(x); sigma.hat = sqrt((length(x) - 1) * var(x) / length(x))}
#'
#' coxsnell.bc(density = pdf, logdensity = lpdf, n = length(x), parms = c("mu", "sigma"),
#'  mle = c(mu.hat, sigma.hat), lower = '-Inf', upper = 'Inf')
#'
#' ################################################################################
#'
#' ## Weibull distribution
#' pdf <- quote(shape / scale ^ shape * x ^ (shape - 1) * exp(-(x / scale) ^ shape))
#' lpdf <- quote(log(shape) - shape * log(scale) + shape * log(x) -
#'  (x / scale) ^ shape)
#'
#' x <- rweibull(n = 100, shape = 1.5, scale = 2.0)
#'
#' fit <- fitdist(data = x, distr = 'weibull')
#' fit$vcov
#'
#' coxsnell.bc(density = pdf, logdensity = lpdf, n = length(x), parms = c("shape", "scale"),
#'  mle = fit$estimate, lower = 0)
#'
#' ################################################################################
#'
#' ## Exponentiated Weibull distribution
#' pdf <- quote(alpha * shape / scale ^ shape * x ^ (shape - 1) * exp(-(x / scale) ^ shape) *
#'  (1 - exp(-(x / scale) ^ shape)) ^ (alpha - 1))
#' lpdf <- quote(log(alpha) + log(shape) - shape * log(scale) + shape * log(x) -
#'  (x / scale) ^ shape + (alpha - 1) * log((1 - exp(-(x / scale) ^ shape))))
#'
#' coxsnell.bc(density = pdf, logdensity = lpdf, n = 100, parms = c("shape", "scale", "alpha"),
#'  mle = c(1.5, 2.0, 1.0), lower = 0)
#'
#' ################################################################################
#'
#' ## Exponetial distribution
#' pdf <- quote(rate * exp(-rate * x))
#' lpdf <- quote(log(rate) - rate * x)
#'
#' x <- rexp(n = 100, rate = 0.5)
#'
#' fit <- fitdist(data = x, distr = 'exp')
#' fit$vcov
#'
#' coxsnell.bc(density = pdf, logdensity = lpdf, n = length(x), parms = c("rate"),
#'  mle = fit$estimate, lower = 0)
#'
#' ################################################################################
#'
#' ## Gamma distribution
#' pdf <- quote(1 /(scale ^ shape * gamma(shape)) * x ^ (shape - 1) * exp(-x / scale))
#' lpdf <- quote(-shape * log(scale) - lgamma(shape) + shape * log(x) -
#'  x / scale)
#'
#' x <- rgamma(n = 100, shape = 1.5, scale = 2.0)
#'
#' fit <- fitdist(data = x, distr = 'gamma', start = list(shape = 1.5, scale =  2.0))
#' fit$vcov
#'
#' coxsnell.bc(density = pdf, logdensity = lpdf, n = length(x), parms = c("shape", "scale"),
#'  mle = fit$estimate, lower = 0)
#'
#' ################################################################################
#'
#' ## Beta distribution
#' pdf <- quote(gamma(shape1 + shape2) / (gamma(shape1) * gamma(shape2)) * x ^ (shape1 - 1) *
#'  (1 - x) ^ (shape2 - 1))
#' lpdf <- quote(lgamma(shape1 + shape2) - lgamma(shape1) - lgamma(shape2) +
#'  shape1 * log(x) + shape2 * log(1 - x))
#'
#' x <- rbeta(n = 100, shape1 = 2.0, shape2 = 2.0)
#'
#' fit <- fitdist(data = x, distr = 'beta', start = list(shape1 = 2.0, shape2 =  2.0))
#' fit$vcov
#'
#' coxsnell.bc(density = pdf, logdensity = lpdf, n = length(x), parms = c("shape1", "shape2"),
#' mle = fit$estimate, lower = 0, upper = 1)
#'
#' @rdname coxsnell.bc
#' @export

coxsnell.bc <- function(density, logdensity, n, parms, mle, lower = '-Inf', upper = 'Inf', ...)
{
  {p <- length(parms); l <- length(mle)};
  if(p != l) stop("The arguments 'parms' and 'mle' must be have the same size")

  kappa_ij <- matrix(NA_real_, ncol = p, nrow = p)
  kappa_ijl <- array(NA_real_, dim = c(p, p, p))
  kappa_ij_l <- kappa_ijl

  {colnames(kappa_ij) <- parms; rownames(kappa_ij) <- parms};

  {integrand <- function(x){}; snd <- integrand}

  for(i in 1:p)
  {
    assign(parms[i], mle[i])
  }

  first <- sapply(1:p, function(i) D(logdensity, parms[i]))
  for(i in 1:p)
  {
    for(j in 1:p)
    {
      second <- D(first[[i]], parms[j])
      body(integrand) <- bquote(.(second) * .(density))
      aux <- tryCatch(integrate(integrand, lower, upper, stop.on.error = FALSE)[c("message", "value")], error = function(e) list(message = "fails"))
      if(aux$message != 'OK') stop('The integrate function failed')
      kappa_ij[i, j] <- -n * aux$value
      for(l in 1:p)
      {
	      third <- D(second, parms[l])
	      body(integrand) <- bquote(.(third) * .(density))
	      aux <- tryCatch(integrate(integrand, lower, upper, stop.on.error = FALSE)[c("message", "value")], error = function(e) list(message = "fails"))
	      if(aux$message != 'OK') stop('The integrate function failed')
	      kappa_ijl[i, j, l] <- n * aux$value

	      body(integrand) <- bquote(.(second) * .(first[[l]]) * .(density))
	      aux <- tryCatch(integrate(integrand, lower, upper, stop.on.error = FALSE)[c("message", "value")], error = function(e) list(message = "fails"))
	      if(aux$message != 'OK') stop('The integrate function failed')
	      kappa_ij_l[i, j, l] <- n * aux$value
      }
    }
  }
  if(any(eigen(kappa_ij)$values < 0)) stop("The final Hessian matrix has at least one negative eigenvalue")
  inv_kappa_ij <- solve(kappa_ij)

  bc <- vector(length = p, mode = 'numeric')
  names(bc) <- parms
  for(s in 1:p)
  {
    bc[s] <- 0
    for(i in 1:p)
    {
      for(j in 1:p)
      {
	      for(l in 1:p)
	      {
	        bc[s] <- bc[s] + inv_kappa_ij[s, i] * inv_kappa_ij[j, l] * (0.5 * kappa_ijl[i, j, l] + kappa_ij_l[i, j, l])
	      }
      }
    }
  }
  {mle.bc <- mle - bc; varcov.bc <- expected.varcov(density, logdensity, n, parms, mle.bc, lower, upper, ...)$varcov}
  {names(mle) <- parms; names(mle.bc) <- parms};
  return(list(mle = mle, varcov = inv_kappa_ij, mle.bc = mle.bc, varcov.bc = varcov.bc, bias = bc))
}
