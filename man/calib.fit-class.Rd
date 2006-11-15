\name{calib.fit-class}
\docType{class}
\alias{calib.fit-class}
\alias{coefficients,calib.fit-method}
\alias{fitted,calib.fit-method}
\alias{loq,calib.fit-method}
\alias{mdc,calib.fit-method}
\alias{plot,calib.fit,missing-method}
\alias{rdl,calib.fit-method}
\alias{residuals,calib.fit-method}

\title{Class "calib.fit" ~~~ }
\description{	 ~~ A concise (1-5 lines) description of what the class is.  ~~}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("calib.fit", ...)}.
	 ~~ describe objects here ~~ 
}
\section{Slots}{
	 \describe{
    \item{\code{coefficients}:}{Estimates of the coefficients}
    \item{\code{se.coefficients}:}{Estimates of the standard errors for the coefficients}
    \item{\code{sigma}:}{Standard deviation of model}
    \item{\code{cov.unscaled}:}{Unscaled variance-covariance matrix of the coefficients}
    \item{\code{theta}:}{The estimated value of POM parameter theta}
    \item{\code{df.residual}:}{The residual degree of freedom}
    \item{\code{fitted.values}:}{The estimated fitted values}
    \item{\code{residuals}:}{The values of the residuals}
    \item{\code{vfe.method}:}{The type of variance function used to estimate theta}
    \item{\code{kused}:}{The number of iteration to convergence of the GLS fit}
    \item{\code{status}:}{Indicated whether the algorithm converged}
    \item{\code{x}:}{The x values}
    \item{\code{y}:}{The y values}
    \item{\code{parm}:}{Indicated whether a log parameterization of the model used}
    \item{\code{m}:}{The number of repeated measures used in the model}
    \item{\code{cv}:}{Coefficient of variation used}
    \item{\code{mdc}:}{Minimum detectable concentration}
    \item{\code{rdl}:}{Reliable detection limit}
    \item{\code{loq}:}{Limit of quantitation}
    \item{\code{cf}:}{Measure of assay performance}
    \item{\code{gradient}:}{The gradient matrix based off of the final parameter estimates}
    \item{\code{lof.test}:}{An object of class lof.test}
    \item{\code{var.model}:}{The type of variance model used}
    \item{\code{conf.level}:}{Confidence level used}
    \item{\code{mmod}:}{In the linear regression fit whether a quadrattic or linear
  			model was fit}
    \item{\code{SST}:}{Sum of squares total. Linear regression fit only}
    \item{\code{SSE}:}{Sum of squares error. Linear regression fit only}
    \item{\code{SSR}:}{Sum of squares regression. Linear regression fit only}
    \item{\code{Rsq}:}{R squared. Linear regression fit only}
    \item{\code{type}:}{Object of class \code{"character"} ~~ }
    \item{\code{rdlwarn}:}{The general model type used, fpl, thpl or lin}
  }
}
\section{Methods}{
  \describe{
    \item{coefficients}{Accessor function for coefficients}
    \item{fitted}{Accessor function for the fitted values}
    \item{loq}{Accessor function for the limit of quanitization}
    \item{mdc}{Accessor function for the minimum detectable concentration}
    \item{plot}{Plot methods for object of class calib.fit}
    \item{rdl}{Accessor function the reliable detection limit}
    \item{residuals}{Accessor function for residuals}
	 }
}
\references{ ~put references to the literature/web site here ~ }
\author{Perry Haaland, Daniel Samarov, Elaine McVey}

\seealso{\code{\link{calib.fit}}}}

\examples{
showClass("calib.fit")
}
\keyword{classes}
