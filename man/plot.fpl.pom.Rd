\name{plot.fpl.pom}
\alias{plot.fpl.pom}

\title{Plotting function for calib.fit function}
\description{
Provides plotting method for the various fpl fits. 
}
\usage{
plot.fpl.pom(fpl.out, conf = 0.9, m = fpl.out@m, main = "", 
xlab = "", ylab = "", xlim = c(0, max(fpl.out@x)), 
ylim = c(0, 1.1 * max(fpl.out@y, na.rm = T)), logplot = F, 
print.txt = T, click.txt = F, start.txt, increm, 
dig = c(4, 3), txt.size = 1, heading = T, lof.print = F, 
pred.lim = F, ...)
}

\arguments{
  \item{fpl.out}{Object of class calib.fit}
  \item{conf}{Confidence limits to plot}
  \item{m}{Number of replicates}
  \item{main}{Title of graph}
  \item{xlab}{x label}
  \item{ylab}{y label}
  \item{xlim}{x limit}
  \item{ylim}{y limit}
  \item{logplot}{Determines if x variable to be printed on the log scale}
  \item{print.txt}{Should text be printed on the plot}
  \item{click.txt}{...}
  \item{start.txt}{...}
  \item{increm}{Helps determine the position of the text}
  \item{dig}{Graphing parameter determining text size on the plot}
  \item{txt.size}{Controls text size proportion (similar to cex par argument)}
  \item{heading}{...}
  \item{lof.print}{Should the lack of fit statistics be printed}
  \item{pred.lim}{Should prediction limits for confidence intervals be plotted.
  Default is FALSE.}
  \item{\dots}{...}
}
\details{
Note that this function is not directly accessible by the user but is accessed
through plot. See example below. 
}

\author{Perry Haaland, Daniel Samarov, Elaine McVey}

\seealso{\code{\link{calib.fit}}}}
\examples{
data(ELISA)
attach(ELISA)
model <- calib.fit(Concentration, Response, type = "fpl.pom")
plot(model)
}

\keyword{hplot}
