\name{plot.thpl.pom}
\alias{plot.thpl.pom}

\title{Plotting functions for the calib.fit function}
\description{
Provides plotting method for the various thpl fits.
}
\usage{
plot.thpl.pom(thpl.out, conf = 0.9, m = thpl.out@m, main = "", 
xlab = "", ylab = "", xlim = c(0, max(thpl.out@x)), 
ylim = c(0, 1.1 * max(thpl.out@y, na.rm = T)), logplot = F, 
print.txt = T, click.txt = F, start.txt, increm, dig = c(4, 3), 
txt.size = 1, heading = T, lof.print = F, bios = F, pred.lim = F, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{thpl.out}{Object of class thpl.fit}
  \item{conf}{Confidence limits to plot}
  \item{m}{Number of replicates}
  \item{main}{Title of graph}
  \item{xlab}{x label}
  \item{ylab}{y label}
  \item{xlim}{x limit}
  \item{ylim}{y limit}
  \item{logplot}{Determines if x variable to printed on the log scale}
  \item{print.txt}{Should text be printed}
  \item{click.txt}{...}
  \item{start.txt}{...}
  \item{increm}{Helps determine the position of the text}
  \item{dig}{Graphing parameter determining text size of the plot}
  \item{txt.size}{Controls text size proportion (similar to cex par argument)}
  \item{heading}{...}
  \item{lof.print}{Should the lack of fit statistics be printed}
  \item{bios}{If TRUE then instead of having the parameterization beta1,...,
  beta3 to describe the thpl model, Finf, F0 and Kd 
  printed instead}
  \item{pred.lim}{Should prediction limits for confidence limits be plotted.
  Default is FALSE}
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
model <- calib.fit(Concentration, Response, type = "thpl.pom")
plot(model)
}
\keyword{hplot}

