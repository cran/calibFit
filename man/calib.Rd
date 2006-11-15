\name{calib}
\alias{calib}

\title{Calibration function}
\description{
Computes the calibration statistics
}
\usage{
calib(calib.fit.out, y0, conf = 0.9, dilution = 1, 
samp.names, truth, times, samp.units = "", dose.units = "", 
dose.name = "", maxit = 1000, toler = 1e-05, 
rname = "response", extrap = F, xname = x)
}

\arguments{
  \item{calib.fit.out}{Output from calib.fit}
  \item{y0}{ ~~Describe \code{y0} here~~ }
  \item{conf}{A vector of mean response values to predict the unknown x}
  \item{dilution}{Dilution factor}
  \item{samp.names}{Names of the unknowns}
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
\references{Haaland and Davidian 1992}
\author{Perry Haaland, Daniel Samarov, Elaine McVey}

\examples{
data(ELISA)
data(experimentData)
attach(ELISA)
attach(experimentData)

model <- calib.fit(Concentration, Response)
cal <- calib(model, OD)

}
\keyword{models}
