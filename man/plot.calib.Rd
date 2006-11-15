\name{plot.calib}
\alias{plot.calib}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plotting function for calibration}
\description{
Generates a plot calibrated values with various display options
}
\usage{
plot.calib(calib.out, int.meth = "w", 
ylim, xlim = c(min(0, calib.out@labels$plot.ind), 
1.1 * max(calib.out@labels$plot.ind)), ylab, xlab, 
join = calib.out@labels$cont.pl, mdcline = T, txt.size = 1, 
addcal.mat, ..., oma = c(2, 2, 2, 2), cex = 1, 
extrap = calib.out@extrap, main, las = 2, cex.axis = 0.7)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{calib.out}{Object of class calib.fit}
  \item{int.meth}{Indicates which intervals are shown on the plot, Wald, "w",
  or inverse, "i". Default is Wald}
  \item{ylim}{y limit}
  \item{xlim}{x limit}
  \item{ylab}{y label}
  \item{xlab}{x label}
  \item{join}{...}
  \item{mdcline}{Whether a line indicating where the MDC should be plotted. Default is
  TRUE}
  \item{txt.size}{Similar to cex}
  \item{addcal.mat}{...}
  \item{\dots}{...}
  \item{oma}{Graphing parameter}
  \item{cex}{Graphing parameter}
  \item{extrap}{Should values outside the range of concentrations be extrapolated}
  \item{main}{Title of plot}
  \item{las}{Graphing paramter}
  \item{cex.axis}{Graphing parameter}
}
\details{
Note that this function is not directly accessible by the user but is accessed through
the plot functions. See examples below for details.
}

\author{Perry Haaland, Daniel Samarov, Elaine McVey}

\seealso{\code{\link{calib}}}}
\examples{
data(ELISA)
attach(ELISA)
model <- calib.fit(Concentration, Response)
cal <- calib(model, Response)
plot(cal, main = "Plot of Calibrated Values")
}
\keyword{hplot}

