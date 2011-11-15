\name{lin.model}
\alias{lin.model}
\title{Linear or Quadratic model.}
\description{
Internal function, calculates the response values for given set of \code{x} values and the
associated \code{gradient}. 
}
\usage{
lin.model(x, beta, w = 1, type)
}

\arguments{
  \item{x}{Independent variables.}
  \item{beta}{If \code{type = "lin"}, \code{beta} is a single value, if \code{type = "quad"}, \code{beta} is a two values.}
  \item{w}{Weight vector.}
  \item{type}{\code{lin} or \code{quad}.}
}
\author{Perry Haaland, Elaine McVey, Daniel Samarov}
\seealso{\code{\link{calib.fit}}}
\keyword{models}

