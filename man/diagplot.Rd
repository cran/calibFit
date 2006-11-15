\name{diagplot}
\alias{diagplot}

\title{Diagnostic plots for standard curve models in calib}
\description{
Diagnostic plot for the POM fit
}
\usage{
diagplot(pom.out, ptype = "all", ..., cex = 1, oma = c(1, 1, 3, 1))
}

\arguments{
  \item{pom.out}{Output from calib.fit with var.model="POM"}
  \item{ptype}{Different diagnostic plot types to display. Takes values
  "all" for diagnostic plots, "stud" for plot of predicted values versus 
  studentized residuals, 
  "raw" for plot of predicted values versus residuals, "abs" for plot of
  log predicted values versus log absolute residuals , 
  "cube" for a plot of log predicted values versus cube root squared 
  student residuals}
  \item{cex}{Graphing parameter}
  \item{oma}{Graphing parameter}
  \item{\dots}{ ~~Describe \code{\dots} here~~ }

}
\details{
Note that diagplot is accessed through the plot function. See examples
below
}

\author{Perry Haaland, Daniel Samarov, Elaine McVey}

\seealso{\code{\link{calib}}}}
\examples{
data(HPLC)
attach(HPLC)
model <- calib.fit(Concentration, Response)
plot(model, type = "diagplot")
}
\keyword{hplot}
