\name{calib.fit-class}
\docType{class}
\alias{calib.fit-class}
\alias{coef,calib.fit-method}
\alias{coefficients,calib.fit-method}
\alias{fitted,calib.fit-method}
\alias{loq,calib.fit-method}
\alias{mdc,calib.fit-method}
\alias{plot,calib.fit,missing-method}
\alias{print,calib.fit-method}
\alias{rdl,calib.fit-method}
\alias{resid,calib.fit-method}
\alias{residuals,calib.fit-method}
\alias{show,calib.fit-method}
\alias{summary,calib.fit-method}

\title{General Standard Curve Fitting}
\description{calib.fit utilizes two main model types, linear and logistic regression
  for the purpose of standard curve fitting. It also incorporates
  several alternatives within each model type to allow for more
  flexible and reliable fitting.}
%\section{Objects from the Class}{
%Objects can be created by calls of the form \code{new("calib.fit", ...)}.
%	 ~~ describe objects here ~~ 
%}
\section{Slots}{
	 \describe{
    \item{\code{coefficients}:}{Estimates of the coefficients}
    \item{\code{se.coefficients}:}{Estimates of the standard errors for the coefficients}
    \item{\code{sigma}:}{Standard deviation of model}
    \item{\code{cov.unscaled}:}{Unscaled variance-covariance matrix of the coefficients}
    \item{\code{pom}:}{Whether or not POM model was used, TRUE or FALSE.}
    \item{\code{theta}:}{The estimated value of POM parameter theta}
    \item{\code{df.residual}:}{The residual degree of freedom}
    \item{\code{fitted.values}:}{The estimated fitted values}
    \item{\code{residuals}:}{The values of the residuals}
    \item{\code{method}:}{Which algorithm was used to do the optimization (i.e. ML, RML, EM, etc.)}
    \item{\code{kused}:}{The number of iteration to convergence of the GLS fit}
    \item{\code{status}:}{Indicated whether the algorithm converged}
    \item{\code{x}:}{The x values}
    \item{\code{y}:}{The y values}
    \item{\code{logParm}:}{Indicated whether a log parameterization of the model used}
    \item{\code{m}:}{The number of repeated measures used in the model}
    \item{\code{cv}:}{Coefficient of variation used}
    \item{\code{mdc}:}{Minimum detectable concentration}
    \item{\code{rdl}:}{Reliable detection limit}
    \item{\code{loq}:}{Limit of quantitation}
    \item{\code{gradient}:}{The gradient matrix based off of the final parameter estimates}
    \item{\code{lof.test}:}{An slot returing information from a lack of fit (LOF) test.}
    \item{\code{var.model}:}{The type of variance model used}
    \item{\code{conf.level}:}{Confidence level used}
    \item{\code{type}:}{Whether, fpl, tpl, thpl, lin or quad was used. }
    \item{\code{rdlwarn}:}{The general model type used, fpl, thpl or lin}
  }
}
\section{Methods}{
  \describe{
    \item{coef}{Accessor function for coefficients}
    \item{coefficients}{Accessor function for coefficients}
    \item{fitted}{Accessor function for the fitted values}
    \item{loq}{Accessor function for the limit of quanitization}
    \item{mdc}{Accessor function for the minimum detectable concentration}
    \item{plot}{Plot methods for object of class calib.fit}
    \item{print}{Print results from calibration curve fit.}
    \item{rdl}{Accessor function the reliable detection limit}
    \item{resid}{Accessor function for residuals}
    \item{residuals}{Accessor function for residuals}
    \item{show}{Same as print}
    \item{summary}{Same as print}
	 }
}
\author{Perry Haaland, Elaine McVey, Daniel Samarov}

\seealso{
\code{\link{calib-class}, \link{calib.fit}, \link{calib}, \link{plot}, \link{resid}, \link{residuals}, 
\link{show}, \link{summary}, \link{print}, \link{fitted}, \link{coefficients},
\link{coef}}
}
\examples{
showClass("calib.fit")
}
\keyword{classes}
