\name{calib}
\alias{calib}
\title{Calibration function}
\description{
Computes the calibration statistics
}
\usage{
calib(x, y0, conf = 0.9, dilution = 1, samp.names, m = x@m, truth, 
times, samp.units = "", dose.units = "", dose.name = "", 
maxit = 1000, toler = 1e-05, rname = "response", extrap = FALSE, xname = x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Output from calib.fit.}
  \item{y0}{Points to be calibrated.}
  \item{conf}{A vector of mean response values to predict the unknown x.}
  \item{dilution}{Dilution factor.}
  \item{samp.names}{Names of the unknowns.}
  \item{m}{Number of replicates.}
  \item{truth}{Optional argument to provide true concentrations if they are known}
  \item{times}{...}
  \item{samp.units}{Names of the unknowns}
  \item{dose.units}{Units of dose}
  \item{dose.name}{Name of dose}
  \item{maxit}{Maximum number of iterations to use in optimization}
  \item{toler}{Tolerance for optimization step}
  \item{rname}{This is the name of the reponse variable}
  \item{extrap}{Option to extrapoloate out of range values}
  \item{xname}{Names of concentrations}
}

\value{
   \item{Estimated.x}{Predicted values of x (for example concentration)}
  \item{PredStdErr}{The predicted standard errors of the estimated x's}
  \item{inver.low}{The estimate of the lower confidence limit for the 
  	predicted x's using inverse estimation}
  \item{inver.up}{The estimate of the upper confidence limit for the
  	predicted x's using inverse estimation}
  \item{wald.low}{The estimate of the lower confidence limit for the
  	predicted x's using Wald estimation}
  \item{wald.up}{The estimate of the upper confidence limit for the
  	predicted x's using Wald estimation}
  \item{avg.response}{y0 values}
}
\references{Davidian and Haaland 1990}
\author{Perry Haaland, Elaine McVey, Daniel Samarov}

\seealso{\code{\link{calib-class}, \link{calib.fit}, \link{calib.fit-class}, \link{plot}}}
\examples{
data(HPLC)
attach(HPLC)
model <- calib.fit(Concentration, Response)
calib(model, Concentration)
}
\keyword{models}
