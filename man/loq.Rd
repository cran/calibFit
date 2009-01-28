\name{loq}
\alias{loq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Limit of Quantization}
\description{
Computes \code{loq} value.
}
\usage{
loq(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{calib.fit}}
  \item{\dots}{NA}
}
\value{
 
  \item{value}{Limit of quantization}
  
}

\author{Perry Haaland, Elaine McVey, Daniel Samarov}

\examples{
data(ELISA)
attach(ELISA)
model <- calib.fit(Concentration, Response)
loq(model)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}

