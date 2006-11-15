\name{calib.fit}
\alias{calib.fit}

\title{General Standard Curve Fitting}
\description{
 calib.fit utilizes two main model types, linear and logistic,
  for the purpose of standard curve fitting. It also incorporates
  several alternatives within each model type to allow for more
  flexible and reliable fitting.
}
\usage{
calib.fit(x, y, m = 2, cv = 0.2, conf = 0.95, b3start, 
b4start, mx = 50, ap.calc = T, lof.calc = T, theta, 
method = "squared", kmax = 10, theta.init = 1, digits = 3, 
type = c("log.fpl.pom", "fpl.pom", "log.fpl", "fpl", 
"lin.quad.pom", "log.thpl.pom", "thpl.pom", "log.thpl", 
"thpl", "lin.pom"), intercept = T)
}

\arguments{
\item{x}{Independent variable (for example dose)}
  \item{y}{Dependent variable (for example response)}
  \item{m}{Number of repeated measurements}
  \item{cv}{The acceptable coefficient of variation. The limits of quantitation
  			are calculated with this constraint}
  \item{conf}{The confidence level used for the determining the prediction interval}
  \item{b3start}{The starting value for b3}
  \item{b4start}{The starting value for b4}
  \item{mx}{The maximum number of iterations used in the non-linear least-squares
  			fit}
  \item{ap.calc}{Should assay performance measures (MDC, RDL and LOQ) 
  			be calculated. The default is TRUE.}
  \item{lof.calc}{Should the lack of fit statistics be calculated. The 
  			default is TRUE.}
  \item{theta}{If a value is specified for theta then calib.fit will not 
  			estimate theta and will only perform the iterative reweighting 
  			of the response.}
  \item{method}{The loss function used in the optimization process when estimating
  			theta. The values for method can be "absolute", "squared", "ols",
  			or "reml". The default is "squared".}
  \item{kmax}{The maximum number of iterations used to estimate theta. The
  			default is 10.}
  \item{theta.init}{The initial starting value for theta.}
  \item{digits}{A graphing parameter specifying the number of digits which will 
  			be displayed when an object of class calib.fit is plotted.}
  \item{type}{Can take the values "log.fpl.pom", the log parameterized four
  			parameter logistic regression (fpl) fit with power of the mean (POM),
  			"fpl.pom", fpl fit with POM, "log.fpl", log parameterized fpl, "fpl",
  			standard fpl, "thpl.pom", three parameter logistic (thpl) regression fit
  			with POM, thpl without POM "thpl", log parameterized POM thpl,
  			"log.thpl.pom", log parameterized thpl, "log.thpl", "lin.quad.pom", 
  			linear regression with 
  			a quadrattic term fit by POM, "lin.pom", linear regression fit by POM.}
  \item{intercept}{Whether or not to include an intercept in a model of type "lin.quad.pom"
  			or "lin.pom". Default is TRUE.}
}

\value{
  \item{coefficients}{Estimates of the coefficients}
  \item{se.coefficients}{Estimates of the standard errors for the coefficients}
  \item{cov.unscaled}{Unscaled variance-covariance matrix of the coefficients}
  \item{theta} {The estimated value of POM parameter theta}
  \item{df.residuals}{The residual degree of freedom}
  \item{fitted.values}{The estimated fitted values}
  \item{residuals}{The values of the residuals}
  \item{vfe.method}{The type of variance function used to estimate theta}
  \item{kused}{The number of iteration to convergence of the GLS fit}
  \item{status}{Indicated whether the algorithm converged}
  \item{x}{The x values}
  \item{y} {The y values}
  \item{parm}{Indicated whether a log parameterization of the model used}
  \item{m}{The number of repeated measures used in the model}
  \item{cv}{Coefficient of variation used}
  \item{mdc}{Minimum detectable concentration}
  \item{rdl}{Reliable detection limit}
  \item{loq}{Limit of quantitation}
  \item{cf}{Measure of assay performance}
  \item{gradient}{The gradient matrix based off of the final parameter estimates}
  \item{Fstat}{The F statistic for the goodness of fit test}
  \item{p.value}{The p-value for the goodness of fit test}
  \item{lofss}{Lack of fit sum of squares}
  \item{df.lof}{Lack of fit degrees of freedom}
  \item{pure.error}{Lack of fit error deviation from the mean}
  \item{df.pure.error}{Degrees of freedom for pure error}
  \item{sse}{Lack of fit error deviation from the fitted values}
  \item{sse}{Degrees of freedom for sse}
  \item{var.model}{The type of variance model used}
  \item{conf.level}{Confidence level used}
  \item{mmod}{In the linear regression fit whether a quadrattic or linear
  			model was fit.}
  \item{SST}{Sum of squares total. Linear regression fit only.}
  \item{SSE}{Sum of squares error. Linear regression fit only.}
  \item{SSR}{Sum of squares regression. Linear regression fit only.}
  \item{Rsq}{R squared. Linear regression fit only.}
  \item{type}{The general model type used, fpl, thpl or lin.}
}
\references{Davidian and Haaland 1992}
\author{Perry Haaland, Daniel Samarov, Elaine McVey}

\examples{
data(HPLC)
attach(HPLC)
model <- calib.fit(Concentration, Response)
}

\keyword{models}

