\name{thpl.model}
\alias{thpl.model}
\title{Three Parameter Logistic Model}
\description{
Internal function, calculates the response values for given set of \code{x} values and the
associated \code{gradient}. 
}
\usage{
thpl.model(x, b1, b2, b3, w = 1, logParm = TRUE)
}

\arguments{
  \item{x}{Independent variables.}
  \item{b1}{Value for \code{b1}.}
  \item{b2}{Value for \code{b2}.}
  \item{b3}{Value for \code{b3}.}
  \item{w}{Weight vector.}
  \item{logParm}{Should log parameterization be used.}
}

\author{Perry Haaland, Elaine McVey, Daniel Samarov}
\seealso{\code{\link{calib.fit}}}
\keyword{models}
