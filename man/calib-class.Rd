\name{Calib-class}
\docType{class}
\alias{Calib-class}
\alias{plot,Calib,missing-method}

\title{Class "Calib" ~~~ }
\description{	 ~~ A concise (1-5 lines) description of what the class is.  ~~}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Calib", ...)}.
	 ~~ describe objects here ~~ 
}
\section{Slots}{
	 \describe{
    \item{\code{Estimated.x}:}{Predicted values of x}
    \item{\code{PredStdErr}:}{Predicted standard errors}
    \item{\code{inver.low}:}{Estimated lower confidence intervals using inverse method}
    \item{\code{inver.up}:}{Estimated upper confidence intervals using inverse method}
    \item{\code{wald.low}:}{Esimated lower confidence interval using Wald intervals}
    \item{\code{wald.up}:}{Estimated upper confidence interval using Wald intervals}
    \item{\code{waldl.low}:}{Internal use}
    \item{\code{waldl.up}:}{Internal use}
    \item{\code{avg.response}:}{y0 values}
    \item{\code{dilution}:}{Dilution factors}
    \item{\code{oor}:}{Values which fall out of range of assay}
    \item{\code{lmdc}:}{Internal use}
    \item{\code{row.names}:}{Internal use}
    \item{\code{labels}:}{Internal use}
    \item{\code{max.x}:}{Internal use}
    \item{\code{extrap}:}{Should values outside the concentration range be estimated}
    \item{\code{repeq}:}{Internal use}
    \item{\code{rname}:}{Internal use}
    \item{\code{conf.level}:}{Values to be used for estimating confidence intervals of
    esimated concentrations}
    \item{\code{type}:}{Model type used (fpl, thpl or lin)}
    \item{\code{mdc}:}{Minimum detectable concentration values}
    \item{\code{truth}:}{If true concentration values are provided 
    in the calib function}
    \item{\code{times}:}{Internal use}
    \item{\code{samp.names}:}{Names of samples}
    \item{\code{lof.test}:}{Lack of fit object}
    \item{\code{calib.fit}:}{calib.fit object}
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "Calib", y = "missing")}: ... }
	 }
}

\author{Perry Haaland, Daniel Samarov, Elaine McVey}

\seealso{\code{\link{calib}}}}
%\seealso{\code{\link{calib.fit}}}}

\examples{
showClass("Calib")
}
\keyword{classes}
