\name{rdl}
\alias{rdl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Reliable Detection Limit}
\description{
Methods for function \code{rdl}
}
\usage{
rdl(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{calib.fit}}
  \item{\dots}{NA}
}
\value{
 
  \item{value}{Reliable detection limit}
  
}

\author{Perry Haaland, Elaine McVey, Daniel Samarov}

\examples{
data(ELISA)
attach(ELISA)
model <- calib.fit(Concentration, Response)
rdl(model)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}

