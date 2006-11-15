\name{plot.lin.pom}
\alias{plot.lin.pom}

\title{Plotting function for calib.fit function}
\description{
Provides plotting method for the various lin fits.   
}
\usage{
plot.lin.pom(lin.out, conf = 0.9, m = lin.out@m, main = "", 
xlab = "", ylab = "", xlim = c(0, max(lin.out@x)), 
ylim = c(0, 1.1 * max(lin.out@y, na.rm = T)), print.txt = T, 
start.txt, increm = 0.03, digits = c(2, 2, 4, 2, 2), txt.size = 0.8, 
lof.print = F, logplot = F, pred.lim = F, sub = "", ...)
}

\arguments{
  \item{lin.out}{Object of class calib.fit}
  \item{conf}{Confidence limits to plot}
  \item{m}{Number of replicates}
  \item{main}{Title of graph}
  \item{xlab}{x label}
  \item{ylab}{y label}
  \item{xlim}{x limit}
  \item{ylim}{y limit}
  \item{print.txt}{...}
  \item{start.txt}{...}
  \item{increm}{Helps determine the position of the text}
  \item{digits}{Graphing parameter determining text size on the plot}
  \item{txt.size}{Controls text size proportion (similar to cex par
  argument)}
  \item{lof.print}{Should the lack of fit statistics be printed}
  \item{logplot}{Determines if x variable to printed on the log scale}
  \item{pred.lim}{Should prediction limits for confidence limits be
  plotted. Default is FALSE.}
  \item{sub}{Subtitle}
  \item{\dots}{...}
}
\details{
Note that this function is not directly accessible by the user but is accessed
through plot. See example below. 
}

\author{Perry Haaland, Daniel Samarov, Elaine McVey}

\seealso{\code{\link{calib.fit}}}}
\examples{
data(HPLC)
attach(HPLC)
model <- calib.fit(Concentration, Response, type = "lin.pom")
plot(model)
}
\keyword{hplot}
