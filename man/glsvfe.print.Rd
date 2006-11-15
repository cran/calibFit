\name{glsvfe.print}
\alias{glsvfe.print}
\title{Printing options for calib}
\description{
This function is called when running calib.fit and prints out 
the results from the gls and vfe from each iteration
}
\usage{
glsvfe.print(k, b, se.b, sigma, theta, digits, mod = 4, bios = T)
}

\arguments{
  \item{k}{The kth iteration}
  \item{b}{The vector of current estimates ofr b1,b2,b3,b4}
  \item{se.b}{The standard error of the estimates}
  \item{sigma}{The residual standard error}
  \item{theta}{Current estimate of the POM parameter}
  \item{digits}{Used to determilne the number of significant digits. It is set to 3
  in calib.fit}
  \item{mod}{4 if fpl is being used, 3 for thpl}
  \item{bios}{}
}

\author{Perry Haaland, Matthew Mitchell, Daniel Samarov}

\examples{
# Example
}
\keyword{print}

