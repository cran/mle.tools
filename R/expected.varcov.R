#' @importFrom stats D
#' @importFrom stats integrate
#'
#' @name expected.varcov
#' @aliases expected.varcov
#'
#' @title Expected Fisher Information
#'
#' @description \code{expected.varcov} calculates the inverse of the expected Fisher information. Analytical second-order partial log-density derivatives and numerical integration are used in the calculations.
#'
#' @author Josmar Mazucheli \email{jmazucheli@gmail.com}
#'
#' @param density An expression with the probability density function.
#' @param logdensity An expression with the log of the probability density function.
#' @param n A numeric scalar with the sample size.
#' @param parms A character vector with the parameter name(s) specified in the density and logdensity expressions.
#' @param mle A numeric vector with the parameter estimate(s).
#' @param lower The lower integration limit (lower = ``-Inf'' is the default).
#' @param upper The upper integration limit (upper = ``Inf'' is the default).
#' @param ... Additional arguments passed to \code{integrate} function.
#'
#' @return \code{expected.varcov} returns a list with two components (i) \bold{mle}: the inputted maximum likelihood estimate(s) and (ii) \bold{varcov}: the expected variance-covariance evaluated at the inputted mle argument.
#'
#' @return If the numerical integration fails and/or the expected information is singular an error message is returned.
#
#' @seealso \code{\link[stats]{deriv}}, \code{\link[stats]{D}}, \code{\link[stats]{integrate}}, \code{\link[mle.tools]{expected.varcov}}.
#'
#' @details The second-order partial log-density derivatives and its expected values are calculated via \code{D} and \code{integrate} functions, respectively.
#'
#' @examples
#' {library(mle.tools); library(fitdistrplus); set.seed(1)};
#'
#' ## Normal distribution
#' pdf <- quote(1 / (sqrt(2 * pi) * sigma) * exp(-0.5 / sigma ^ 2 * (x - mu) ^ 2))
#' lpdf <- quote(-log(sigma) - 0.5 / sigma ^ 2 * (x - mu) ^ 2)
#'
#' x <- rnorm(n = 100, mean = 0.0, sd = 1.0)
#'
#' expected.varcov(density = pdf, logdensity = lpdf, n = length(x), parms = c("mu", "sigma"),
#'  mle = c(mean(x), sd(x)), lower = '-Inf', upper = 'Inf')
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
#' expected.varcov(density = pdf, logdensity = lpdf, n = length(x), parms = c("shape", "scale"),
#'  mle = fit$estimate, lower = 0)
#'
#' ################################################################################
#'
#' ## Expoentiated Weibull distribution
#' pdf <- quote(alpha * shape / scale ^ shape * x ^ (shape - 1) * exp(-(x / scale) ^ shape) *
#'  (1 - exp(-(x / scale) ^ shape)) ^ (alpha - 1))
#' lpdf <- quote(log(alpha) + log(shape) - shape * log(scale) + shape * log(x) -
#'  (x / scale) ^ shape + (alpha - 1) * log((1 - exp(-(x / scale) ^ shape))))
#'
#' expected.varcov(density = pdf, logdensity = lpdf, n = 100, parms = c("shape", "scale", "alpha"),
#'  mle = c(1.5, 2.0, 1.0), lower = 0)
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
#' expected.varcov(density = pdf, logdensity = lpdf, n = length(x), parms = c("rate"),
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
#' expected.varcov(density = pdf, logdensity = lpdf, n = length(x), parms = c("shape", "scale"),
#'  mle = fit$estimate, lower = 0)
#'
#' ################################################################################
#'
#' ## Beta distribution
#' pdf <- quote(gamma(shape1 + shape2) / (gamma(shape1) * gamma(shape2)) * x ^ (shape1 - 1) *
#' (1 - x) ^ (shape2 - 1))
#' lpdf <- quote(lgamma(shape1 + shape2) - lgamma(shape1) - lgamma(shape2) +
#'  shape1 * log(x) + shape2 * log(1 - x))
#'
#' x <- rbeta(n = 100, shape1 = 2.0, shape2 = 2.0)
#'
#' fit <- fitdist(data = x, distr = 'beta', start = list(shape1 = 2.0, shape2 =  2.0))
#' fit$vcov
#'
#' expected.varcov(density = pdf, logdensity = lpdf, n = length(x), parms = c("shape1", "shape2"),
#'  mle = fit$estimate, lower = 0, upper = 1)
#'
#' @rdname expected.varcov
#' @export

expected.varcov <- function(density, logdensity, n, parms, mle, lower = '-Inf', upper = 'Inf', ...)
{
  {p <- length(parms); l <- length(mle)};
  if(p != l) stop("The arguments 'parms' and 'mle' must be have the same size")

  H <- matrix(NA_real_, ncol = p, nrow = p)
  {colnames(H) <- parms; rownames(H) <- parms; names(mle) = parms};

  integrand <- function(x){}

  for(i in 1:p)
  {
    assign(parms[i], mle[i])
  }

  first <- sapply(1:p, function(i) D(logdensity, parms[i]))

  for(i in 1:p)
  {
    for(j in i:p)
    {
      second <- D(first[[i]], parms[j])
      body(integrand) 	<- bquote(.(second) * .(density))
      aux <- tryCatch(integrate(integrand, lower, upper, stop.on.error = FALSE)[c("message", "value")], error = function(e) list(message = "fails"))
      if(aux$message != 'OK') stop('The integrate function failed')
      H[i, j] <- -n * aux$value
      if(j > i) H[j, i] <- H[i, j]
    }
  }
  if(any(eigen(H)$values < 0)) stop("The final Hessian matrix has at least one negative eigenvalue")
  else return(list(mle = mle, varcov = solve(H)))
}

