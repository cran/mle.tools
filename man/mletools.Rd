\name{mle.tools-package}
\alias{mle.tools-package}
\docType{package}
\title{Overview of the ``mle.tools'' Package}

\description{The current version of the \pkg{mle.tools} package has implemented three functions which are of great interest in maximum likelihood estimation. These functions calculates the expected /observed Fisher information and the bias-corrected maximum likelihood estimate(s) using the bias formula introduced by Cox and Snell (1968). They can be applied to any probability density function whose terms are available in the derivatives table of \code{D} function (see ``deriv.c'' source code for further details). Integrals, when required, are computed numerically via \code{integrate} function. Below are some mathematical details of how the returned values are calculated.

Let \eqn{X_{1},\ldots ,X_{n}} be \emph{i.i.d.} random variables with
probability density functions \eqn{f(x_{i}\mid \bold{\theta })} depending on a \eqn{p}-dimensional parameter vector \eqn{\bold{\theta } = (\theta_1,\ldots,\theta_p)}. The \emph{(j,k)-th} element of the observed, \eqn{H_{jk}}, and expected, \eqn{I_{jk}}, Fisher information are calculated, respectively,  as
\deqn{H_{jk} =\left. {-\sum\limits_{i=1}^{n}\frac{%
\partial ^{2}}{\partial \theta _{j}\partial \theta _{k}}\log f\left(
x_{i}\mid {\bold{\theta} }\right) }\right\vert _{\bold{\theta }=\widehat{\bold{%
\theta }}}}
and
\deqn{I_{jk}=-n\times E\left( \frac{\partial ^{2}}{\partial \theta _{j}\partial
\theta _{k}}\log f\left( x\mid \bold{\theta }\right) \right) =\left. -n\times
\int\limits_{\mathcal{X} }\frac{\partial ^{2}}{\partial \theta _{j}\partial
\theta _{k}}\log f\left( x\mid \bold{\theta }\right) \times f\left(
x\mid \bold{\theta }\right) dx\right\vert _{\bold{\theta }=\widehat{\bold{%
\theta }}}}

where \eqn{(j,k=1,\ldots,p)}, \eqn{\bold{\widehat{\theta}}} is the maximum likelihood estimate of \eqn{\bold{\theta}} and \eqn{\mathcal{X}} denotes the support of the random variable \eqn{X}.

The \code{observed.varcov} function returns the inputted maximum likelihood estimate(s) and the inverse of \eqn{\bold{H}} while the \code{expected.varcov} function returns the inputted maximum likelihood estimate(s) and the inverse of \eqn{\bold{I}}. If \eqn{\bold{H}} and/or \eqn{\bold{I}} are singular an error message is returned.

Furthermore, the bias corrected maximum likelihood estimate of \eqn{\theta_s}   (\eqn{s=1,\ldots,p)}, denoted by \eqn{\widetilde{\theta_s}}, is
calculated as \eqn{\widetilde{\theta_s} = \widehat{\theta} - \widehat{Bias}(\widehat{\theta}_s)}, where \eqn{\widehat{\theta}_s} is the maximum likelihood estimate of \eqn{{\theta}_s} and
\deqn{{\widehat{Bias}\left( {\widehat{\theta }}_{s}\right) =}\left. {%
\sum\limits_{j=1}^{p}\sum\limits_{k=1}^{p}\sum\limits_{l=1}^{p}\kappa
^{sj}\kappa ^{kl}\left[ 0.5\kappa _{{jkl}}+\kappa _{{jk,l}}\right] }%
\right\vert _{\bold{\theta }=\widehat{\bold{\theta }}}} where \eqn{\kappa ^{jk}} is the \emph{(j,k)-th} element of the inverse of the expected Fisher information, \eqn{{\kappa_{jkl}=} n\times E\left( \frac{\partial ^{3}}{\partial \theta _{j}\partial {{\theta}}_{k}{\theta }_{l}}\log f\left( x\mid \bold{\theta }\right) \right)} and
\eqn{\kappa_{jk,l}= n \times E\left( \frac{\partial ^{2}}{\partial \theta _{j}\partial \theta_{k}}\log f\left( x\mid\bold{\theta }\right) \times \frac{\partial }{{\theta }_{l}}\log f\left( x\mid\bold{\theta }\right) \right) }.

The bias-corrected maximum likelihood estimate(s) and some other quantities are calculated via \code{coxsnell.bc} function. If the numerical integration fails
and/or \eqn{\bold{I}} is singular an error message is returned.

It is noteworthy that for a series of probability distributions it is possible, after extensive algebra, to obtain the analytical expressions for \eqn{Bias({\widehat{\theta}_s)}}. In Stosic and Cordeiro (2009) are the analytic expressions for 22 two-parameter continuous probability distributions. They also present the \emph{Maple} and \emph{Mathematica} scripts used to obtain all analytic expressions (see Cordeiro and Cribari-Neto 2014 for further details).
}

\author{
Josmar Mazucheli \email{jmazucheli@gmail.com}
}

\references{

Azzalini, A. (1996). \emph{Statistical Inference: Based on the Likelihood}. London: Chapman and Hall.

Cordeiro, G. M. and Cribari-Neto, F., (2014). An introduction to Bartlett correction and bias reduction. SpringerBriefs in Statistics, New-York.

Cordeiro, G. M. and McCullagh, P., (1991). Bias correction in generalized linear models. \emph{Journal of the Royal Statistical Society, Series B}, \bold{53}, 3, 629--643.

Cox, D. R. and Hinkley, D. V. (1974). \emph{Theoretical Statistics}. London: Chapman and Hall.

Cox, D. R. and Snell, E. J., (1968). A general definition of residuals (with discussion). \emph{Journal of the Royal Statistical Society, Series B}, \bold{30}, 2, 24--275.

Efron, B. and Hinkley, D. V. (1978). Assessing the accuracy of the maximum likelihood estimator: Observed versus expected Fisher information. \emph{Biometrika}, \bold{65}, 3, 457--482.

Pawitan, Y. (2001). \emph{In All Likelihood: Statistical Modelling and Inference Using Likelihood}. Oxford: Oxford University Press.

Stosic, B. D. and Cordeiro, G. M., (2009). Using Maple and Mathematica to derive bias corrections for two parameter distributions. \emph{Journal of Statistical Computation and Simulation}, \bold{79}, 6, 751--767.

}
