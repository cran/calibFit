\name{calib.fit}
\alias{calib.fit}
\title{General Standard Curve Fitting}
\description{
 calib.fit utilizes two main model types, linear and logistic,
  for the purpose of standard curve fitting. It also incorporates
  several alternatives within each model type to allow for more
  flexible and reliable fitting.}
\usage{
calib.fit(x, y, b1start, b2start, b3start, b4start, calcDiagnostics = TRUE, 
m, cv = 0.2, conf = 0.95, mx = 50, lof.calc = T, lowLim = 0.001, 
type = c("log.fpl.pom", "fpl.pom", "log.fpl", "fpl", "log.tpl.pom", "tpl.pom", 
"log.tpl", "tpl", "quad.pom", "log.thpl.pom", "thpl.pom", "log.thpl", 
"thpl", "lin.pom"))
}

\arguments{
  \item{x}{Independent variable (for example dose)}
  \item{y}{Dependent variable (for example response)}
  \item{b1start}{Starting value for \code{b1} in either the \code{fpl}, \code{tpl} or \code{thpl} models.}
  \item{b2start}{Starting value for \code{b2} in either the \code{fpl}, \code{tpl} or \code{thpl} models.}
  \item{b3start}{Starting value for \code{b3} in either the \code{fpl}, \code{tpl} or \code{thpl} models.}
  \item{b4start}{Starting value for \code{b4} in either the \code{fpl} or \code{tpl} models.}
  \item{calcDiagnostics}{Should diagnostics (i.e. mdc, rdl and loq) be calculated. Default to \code{TRUE}.}
  \item{m}{Number of repeated measurements}
  \item{cv}{The acceptable coefficient of variation. The limits of quantitation 
  are calculated with this constraint}
  \item{conf}{The confidence level used for the determining the prediction interval}
  \item{mx}{The maximum number of iterations used in the non-linear least-squares fit}
  \item{lof.calc}{Should the lack of fit statistics be calculated. The 
  			default is TRUE.}
  \item{lowLim}{If there are \code{x} values equal to zero what (small) positive value be used to
  approximate it. Defaults is 0.001.}
  \item{type}{Can take the values \code{log.fpl.pom}, the log parameterized four
  			parameter logistic regression (\code{fpl}) fit with power of the mean (POM),
  			\code{fpl.pom}, \code{fpl} fit with POM, \code{log.fpl}, log parameterized \code{fpl}, \code{fpl},
  			standard \code{fpl}, \code{thpl.pom}, three parameter logistic (\code{thpl}) regression fit
  			with POM, \code{thpl} without POM \code{thpl}, log parameterized POM \code{thpl},
  			\code{log.thpl.pom}, log parameterized \code{thpl}, \code{log.thpl}, \code{lin.quad.pom}, 
  			linear regression with 
  			a quadrattic term fit by POM, \code{lin.pom}, linear regression fit by POM.}
}

\value{
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
\references{Davidian and Haaland 1990}
\author{Perry Haaland, Elaine McVey, Daniel Samarov}

\seealso{
\code{\link{calib-class}, \link{calib.fit}, \link{calib}, \link{plot}, \link{resid}, \link{residuals}, 
\link{show}, \link{summary}, \link{print}, \link{fitted}, \link{coefficients},
\link{coef}}
}
\examples{
data(HPLC)
attach(HPLC)
model <- calib.fit(Concentration, Response)
}

\keyword{models}

