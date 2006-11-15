\name{precprof}
\alias{precprof}

\title{Plotting function}
\description{
This creates a plot of the dose vs. coefficient of variation (cv) 
output from an object of class calib.fit.
}
\usage{
precprof(calib.fit.out, m = calib.fit.out@m, 
cv = calib.fit.out@cv, vlen = 500, mit = 1000, 
tolder = 0.001, ylim, xlab = "dose", ylab = "CV", 
txt.size = 1)
}

\arguments{
  \item{calib.fit.out}{An object outputed from calib.fit}
  \item{m}{Number of replicates}
  \item{cv}{Acceptable coefficient of variation}
  \item{ylim}{y limit}
  \item{xlab}{x label}
  \item{ylab}{y label}
  \item{txt.size}{Text size}
  \item{tolder}{...}
  \item{vlen}{...}
  \item{mit}{Max iterations}
}
\details{
At the top of the graph the lowest cv acceptable, the working range of the
assay, and the number of reps is printed.
Additionally note that precprof is accessed as a generic function through
the plot function. 
}

\references{Davidian and Haaland 1992}
\author{Perry Haaland, Dnaiel Samarov, Elaine McVey}

\seealso{\code{\link{calib.fit}}}}
\examples{
data(ELISA)
attach(ELISA)
model <- calib.fit(Concentration, Response)

plot(model, type="precprof")

}
\keyword{hplot}
