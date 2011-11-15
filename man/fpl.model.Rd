\name{fpl.model}
\alias{fpl.model}
\title{Four Parameter Logistic Model}
\description{
Internal function, calculates the response values for given set of \code{x} values and the
associated \code{gradient}. 
}
\usage{
fpl.model(x, b1, b2, b3, b4, w = 1, logParm = TRUE)
}

\arguments{
  \item{x}{Independent variables.}
  \item{b1}{Value for \code{b1}.}
  \item{b2}{Value for \code{b2}.}
  \item{b3}{Value for \code{b3}.}
  \item{b4}{Value for \code{b4}.}
  \item{w}{Weight vector.}
  \item{logParm}{Should log parameterization be used.}
}

\author{Perry Haaland, Elaine McVey, Daniel Samarov}
\seealso{\code{\link{calib.fit}}}
\keyword{models}
