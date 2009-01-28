\name{mdc}
\alias{mdc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Minimum Detectable Concentration}
\description{
Returns the minimum detectable concentration based on the model fitted by \code{calib.fit}
}
\usage{
mdc(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{calib.fit}}
  \item{\dots}{NA}
}
\details{
Calculates the minimum detectable concentrations
}
\value{

  \item{value}{\code{mdc} value}
  
}

\author{Perry Haaland, Elaine McVey, Daniel Samarov}

\examples{
data(HPLC)
attach(HPLC)
model <- calib.fit(Concentration, Response)
mdc(model)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}

