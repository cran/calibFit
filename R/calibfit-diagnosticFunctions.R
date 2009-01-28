## calibfit-diagnostic-functions.R
## 
## Functions that are called by other calibfit functions. 
## These are not accessible to the user.
##
## Author: samarov
###############################################################################

## 11-17-08: Have a method for mdc, do not need this anymore. DVS
#mdc.calc <- function(out, conf = out@conf.level, m = out@m){
#	if((out@type == "lin") | (out@type == "quad")){
#		mdc.out <- mdc.lin.pom(out, m = m, conf = conf)
#	}
#	if(out@type == "thpl"){
#		mdc.out <- mdc.thpl.pom(out, m = m, conf = conf)
#	}
#	if(out@type == "fpl"){
#		mdc.out <- mdc.fpl.pom(out, m = m, conf = conf)
#	}
#	if(out@type == "tpl"){
#		mdc.out <- mdc.tpl.pom(out, m = m, conf=conf)
#	}
##	if(!(calib.fit.out@type%in%c("lin","fpl","tpl","thpl")))
##		stop("Model type selected is not one of fpl, thpl, tpl or lin")
#	return(mdc.out)
#}

mdc.fpl.pom <- function(out, conf = out@conf.level, m = out@m)
## 07-01-08: Does not appear to me that these are in use any more. DVS
#	, mit = 1000, toler = 0.0001)
{
	######################################################################
	## This function computes the Mimimum Detectable Concentration      ##
	##  for the FPL  POM model.                                         ##
	##                                                                  ##
	## input variables:                                                 ##
	##   fpl.out: the output from the fitting function                  ##
	##   m: number of reps (number of y-values for each x)              ##
	##   xu: upper x value                                              ##
	##   conf: confidence level (in point-wise intervals)               ##
	##                                                                  ##
	## old comments have ###                                            ##
	## new version created by MM on 4/14/03                             ##
	######################################################################
	
	## reads in the necessary parameters
	cov.un <- out@cov.unscaled
	b <- out@coefficients
	s <- out@sigma
#	xu <- max(out@x)
	
	## mdc is computed by finding the x that corresponds to the upper
	## prediction limit of y when x=0 (for increasing curve - b1 < b2)
	## or the lower prediction limit of y when x=0 (decreasing curve b1 > b2)
	## The predicted value of y when x=0 is b1. So we find the prediction interval for b1.
	
	f0b <- b[1]
	var.f0b <-  s^2 * (((f0b^2)^out@theta)/m + (cov.un[1, 1]))
	
	tcrit <- qt( (1 + conf)/2, out@df.residual)
	lim0 <- f0b + (sign(b[2] - b[1])) * tcrit * sqrt(var.f0b)
	mdcval <- fpl.inverse(out, lim0)
	
	## 07-01-08: I'm not really sure why mdc returns a data frame. I don't think
	## that this is necessary. Having mdc return mdcval from above. DVS
#	mdc.df <- data.frame(mdc = mdcval, m = m)  
#	return(mdc.df)  
	return(mdcval)
}

mdc.thpl.pom <- function(out, conf = out@conf.level, m = out@m){
	######################################################################
	## This function computes the Mimimum Detectable Concentration      ##
	##  for the THPL  POM model.                                        ##
	##                                                                  ##
	## input variables:                                                 ##
	##   thpl.out: the output from the fitting function                 ##
	##   m: number of reps (number of y-values for each x)              ##
	##   xu: upper x value                                              ##
	##   conf: confidence level (in point-wise intervals)               ##
	##   mit, toler: used if dos=F                                      ##
	##                                                                  ##
	## written by Matthew Mitchell on 5/20/03                           ##
	######################################################################
	
	## reads in the necessary parameters
	cov.un <- out@cov.unscaled
	b <- out@coefficients
	s <- out@sigma
	
	## mdc is computed by finding the x that corresponds to the upper
	## prediction limit of y when x=0 (for increasing curve - b1 < b2)
	## or the lower prediction limit of y when x=0 (decreasing curve b1 > b2)
	## The predicted value of y when x=0 is b1. So we find the prediction interval for b1.
	
	f0b <- b[1]
	var.f0b <-  s^2 * (((f0b^2)^out@theta)/m + (cov.un[1, 1]))
	
	tcrit <- qt( (1 + conf)/2, out@df.residual)
	lim0 <- f0b + (sign(b[2] - b[1])) * tcrit * sqrt(var.f0b)
	mdcval <- thpl.inverse(out,lim0)
	
	## 07-01-08: I'm not really sure why mdc returns a data frame. I don't think
	## that this is necessary. Having mdc return mdcval from above. DVS
#	mdc.df <- data.frame(mdc = mdcval, m = m)
#	return(mdc.df)    
	return(mdcval)
}



mdc.tpl.pom <- function(out, conf = out@conf.level, m = out@m){
	######################################################################
	## This function computes the Mimimum Detectable Concentration      ##
	##  for the FPL  POM model.                                         ##
	##                                                                  ##
	## input variables:                                                 ##
	##   fpl.out: the output from the fitting function                  ##
	##   m: number of reps (number of y-values for each x)              ##
	##   xu: upper x value                                              ##
	##   conf: confidence level (in point-wise intervals)               ##
	##   mit, toler: used if dos=F                                      ##
	##   dos: is T or F, used to change the data type                   ##
	##                                                                  ##
	## new version created by MM on 4/14/03                             ##
	######################################################################
	
	## reads in the necessary parameters
	cov.un <- out@cov.unscaled
	b <- out@coefficients
	s <- out@sigma
	
	## mdc is computed by finding the x that corresponds to the upper
	## prediction limit of y when x=0 (for increasing curve - b1 < b2)
	## or the lower prediction limit of y when x=0 (decreasing curve b1 > b2)
	## The predicted value of y when x=0 is b1. So we find the prediction 
	## interval for b1.
	
	f0b <- b[1]
	## 08-08-08: See comment under rdl.tpl.pom. DVS
#	var.f0b <-  s^2 * (((f0b^2)^out@theta) / m)
	var.f0b <- s^2 * (((f0b^2)^out@theta)/m)
	
	tcrit <- qt((1 + conf) / 2, out@df.residual)
	lim0 <- f0b + (sign(b[2] - b[1])) * tcrit * sqrt(var.f0b)
	mdcval <- tpl.inverse(out, lim0)
	if(is.na(mdcval))
		mdcval <- NULL
	
	## 07-01-08: I'm not really sure why mdc returns a data frame. I don't think
	## that this is necessary. Having mdc return mdcval from above. DVS
#	mdc.df <- data.frame(mdc = mdcval, m = m)
#	return(mdc.df)
	return(mdcval)	
}

mdc.lin.pom <- function(out, conf = out@conf.level, m = out@m)
{
	######################################################################
	## This function computes the Mimimum Detectable Concentration      ##
	##  for the Linear POM model.                                       ##
	##                                                                  ##
	## input variables:                                                 ##
	##   lin.out: the output from the fitting function                  ##
	##   m: number of reps (number of y-values for each x)              ##
	##   conf: confidence level (in point-wise intervals)               ##
	##                                                                  ##
	## old comments have ###                                            ##
	## new version created by MM on 4/14/03                             ##
	######################################################################
	min.x <- min(out@x)
	max.x <- max(out@x)
	b <- out@coefficients
	cov.un <- out@cov.unscaled
	s <- out@sigma
	df.den <- out@df.residual
	tcrit <- qt((1 + conf)/2, df.den)      
	f0b <- out@coefficients[1]
	var.f0b <-  s^2 * (((f0b^2)^out@theta)/m + (cov.un[1, 1]))
	
	if(out@type == "lin"){
		lim0 <- f0b + (sign(b[2]))*tcrit*sqrt(var.f0b)
		mdcval <- (lim0 - b[1])/b[2]
	}
	
	## The only valid quadratic models are those whose critical point
	## is less than 0 or higher than max(x)
	## Whether we use the upper or lower interval is 0 is determined
	## by whether the curve is increasing or decreasing - it won't be
	## both because the critical point outside the range.  Hence
	## it is increasing if f(1) > f(0) which means b + 2c > 0
	
	if(out@type == "quad"){
		lim0 <- f0b + (sign(b[2] + 2 * b[3])) * tcrit * sqrt(var.f0b)
		mdc1 <- (-b[2] - sqrt(b[2]^2 - 4 * b[3] * (b[1] - lim0)))/(2 * b[3])
		mdc2 <- (-b[2] + sqrt(b[2]^2 - 4 * b[3] * (b[1] - lim0)))/(2 * b[3])
		mdcval <- ifelse((mdc1 > min.x) & (mdc1 < max.x), mdc1, mdc2)
	}
	
	
	## 06-30-08: This had previously been set up to output a list, however,
	## all of the other mdc function return a data frame. For consistency
	## having it return a data frame. DVS
#	mdcout <- list(mdc = as.vector(mdcval), m = m)  
	### I left this as it was - not sure why
#	return(mdcout)
#	## it is a list here while the other was made into a df
	
	## 07-01-08: I'm not really sure why mdc returns a data frame. I don't think
	## that this is necessary. Having mdc return mdcval from above. DVS
#	mdc.df <- data.frame(mdc = mdcval, m = m)
#	return(mdc.df)
	return(mdcval)
}

## 11-17-08: Have a method for mdc, do not need this anymore. DVS
#mdc.calc <- function(out, conf = out@conf.level, m = out@m){
#rdl.calc <- function(out, conf = out@conf.level, m = out@m){	
#	if((out@type == "lin") | (out@type == "quad")){
#		rdl.out <- rdl.lin.pom(out, m = m, conf = conf) 
#	}
#	if(out@type == "thpl"){
#		rdl.out <- rdl.thpl.pom(out, m = m, conf = conf)
#	}
#	if(out@type == "fpl"){
#		rdl.out <- rdl.fpl.pom(out, m = m, conf = conf)
#	}
#	if(out@type == "tpl"){
#		rdl.out <- rdl.tpl.pom(out, m = m, conf = conf)
#	}
#	rdl.out
#
#}

######################################################################
## This function computes the Reliable Detection Limit              ##
##  for the FPL  POM model.                                         ##
##                                                                  ##
## input variables:                                                 ##
##   fpl.out: the output from the fitting function                  ##
##   m: number of reps (number of y-values for each x)              ##
##   xu: upper x value                                              ##
##   conf: confidence level (in point-wise intervals)               ##
##   mit, toler: used if dos=F                                      ##
##                                                                  ##
## comments added by Matthew Mitchell 3/12/03  have ##              ##
## the original comments have ###                                   ##
######################################################################
rdl.fpl.pom <- function(out, m = out@m, conf = out@conf.level) 
{
	## this function takes output from a regression 
	## and  returns a value for Reliable Detection Limit
	cov.un <- out@cov.unscaled
	b <- out@coefficients
	
	f0b <- b[1]
	var.f0b <- (((f0b^2)^out@theta) * out@sigma * out@sigma)/m + 
			(out@sigma * out@sigma) * cov.un[1,1]
	
	tcrit <- qt((1 + conf)/2, out@df.residual)
	lim0 <- f0b + (sign(b[2] - b[1])) * tcrit * sqrt(var.f0b)
	xl <- fpl.inverse(out, lim0)
	b50 <- ifelse(out@logParm, exp(b[3]), b[3])
	xp <- seq(xl, max(out@x), len = 450)
	
	ypred <- fpl.model(xp, b, logParm = out@logParm)
	qn1 <- ((as.vector(ypred)^2)^out@theta)/m
	
	qn2 <- diag(attributes(ypred)$gradient %*% 
					cov.un %*% t(attributes(ypred)$gradient))
	yp <- as.vector(ypred) - sign(b[2] - b[1]) * tcrit * 
			out@sigma * sqrt(qn1 + qn2)
	## 01-25-07 Added ties="mean" to approx function below to avoid error messages
	## when there are ties. DVS
	rdlval <- approx(yp, xp, lim0, ties = "mean")$y
	
	## 07-01-08: I'm not really sure why rdl returns a data frame. I don't think
	## that this is necessary. Having rdl return rdlval from above. DVS
#	rdlout <- data.frame(rdl = rdlval, m = m)
#	return(rdlout)
	return(rdlval)
	
}

######################################################################
## This function computes the Reliable Detection Limit              ##
##  for the THPL  POM model.                                        ##
##                                                                  ##
## input variables:                                                 ##
##   thpl.out: the output from the fitting function                 ##
##   m: number of reps (number of y-values for each x)              ##
##   xu: upper x value                                              ##
##   conf: confidence level (in point-wise intervals)               ##
##   mit, toler: used if dos=F                                      ##
##                                                                  ##
## written by Matthew W. Mitchell on 5/20/03                        ##
## the original comments have ###                                   ##
######################################################################
rdl.thpl.pom <- function(out, m = out@m, conf = out@conf.level)
{
	## this function takes output from a regression 
	## and  returns a value for Reliable Detection Limit
	cov.un <- out@cov.unscaled
	b <- out@coefficients
	n <- length(out@x)
	
	bigs.vec <- c(cov.un[, 1], cov.un[, 2], cov.un[, 3])
	f0b <- b[1]
	var.f0b <- (((f0b^2)^out@theta) * out@sigma^2)/m + 
			(out@sigma^2) * cov.un[1,1]
	## The computation of the critical value is incorrect.  It should be
	## 1 - (1-conf)/2 = (1+conf)/2.  Correction made by MM on 4/9/03.
	##	tcrit <- qt(conf, thpl.out$df.res)
	tcrit <- qt((1 + conf)/2, out@df.residual)
	lim0 <- f0b + (sign(b[2] - b[1])) * tcrit * sqrt(var.f0b)
	xl <- thpl.inverse(out, lim0)
	
	b50 <- ifelse(out@logParm, exp(b[3]), b[3])
	xp <- c(seq(xl, b50, len = 350), 
			seq(b50, max(out@x), len = 100))
	ypred <- thpl.model(xp, b, logParm = out@logParm)
	qn1 <- ((as.vector(ypred)^2)^out@theta)/m
	## 04-03-07: Changing to above DVS
	#qn1 <- ((as.vector(ypred))^(2 * thpl.out@theta))/m
	qn2 <- diag(attributes(ypred)$gradient %*% cov.un %*% 
					t(attributes(ypred)$gradient))
	yp <- as.vector(ypred) - sign(b[2] - b[1]) * tcrit * 
			out@sigma * sqrt(qn1 + qn2)
	## 01-25-07 Added ties="mean" to approx function below to avoid error messages
	## when there are ties. DVS
	rdlval <- approx(yp, xp, lim0, ties = "mean")$y
	
	## 07-01-08: I'm not really sure why rdl returns a data frame. I don't think
	## that this is necessary. Having rdl return rdlval from above. DVS
#	rdlout <- data.frame(rdl = rdlval, m = m)
#	return(rdlout)
	return(rdlval)
	
}

rdl.tpl.pom <- function(out, m = out@m, conf = out@conf.level)  

{
	## this function takes output from a regression 
	## and  returns a value for Reliable Detection Limit
	cov.un <- out@cov.unscaled
	b <- out@coefficients
	n <- length(out@x)
	
	f0b <- b[1]
	## 08-08-08: Making a correction here. With the TPL model the variance  
	## at zero concentration will just be var(eps) where eps is underlying error,
	## so var.f0b should just be out@sigma^2/m. DVS
#	var.f0b <- (((f0b^2)^out@theta) * out@sigma^2)/m
#	var.f0b <- out@sigma^2/m
	var.f0b <- (((f0b^2)^out@theta) * out@sigma^2)/m
	
	tcrit <- qt((1 + conf)/2, out@df.residual)
	lim0 <- f0b + (sign(b[2] - b[1])) * tcrit * sqrt(var.f0b)
	xl <- tpl.inverse(out, lim0)
	if(is.na(xl)){
		return(NULL)
		break
	}
	
	b50 <- ifelse(out@logParm, exp(b[3]), b[3])
	xp <- seq(xl, max(out@x), len = 450)
	
	ypred <- tpl.model(xp, b, logParm = out@logParm)
	qn1 <- ((as.vector(ypred)^2)^out@theta)/m
	
	qn2 <- diag(attributes(ypred)$gradient %*% 
					cov.un %*% t(attributes(ypred)$gradient))
	yp <- as.vector(ypred) - sign(b[2] - b[1]) * tcrit * 
			out@sigma * sqrt(qn1 + qn2)
	
	rdlval <- approx(yp, xp, lim0, ties = "mean")$y
	
	## 07-01-08: I'm not really sure why rdl returns a data frame. I don't think
	## that this is necessary. Having rdl return rdlval from above. DVS
#	rdlout <- data.frame(rdl = rdlval, m = m)
#	return(rdlout)
	return(rdlval)
}




rdl.lin.pom <- function(out, m = out@m, conf = 0.95, tol = 0.01)
{
	############################################################################
	##	Calculates limit of quantitation (RDL)
	##	with a data set and the output from glsvfe.lin.out()
	##
	##	RDL = Reliable Detection Limit
	##		= lowest concentration that can reliabily be expected
	##		to give a reading above the minimum detectable concentration
	##		- i.e., lower conf. limit at x=RDL > upper conf. limit at x=0
	##	m is # reps of the new unknown
	##	correct df=df.den can be used for common theta/sigma cases
	##		else if df.den is missing, it is gotten from lin.out
	##	conf=.95 is the default that goes along with 90% two-sided
	##		confidence limits about the curve
	##	the default is to get mdc from pom, it can also be given explicitly
	############################################################################
	theta <- out@theta
	std.dev <- out@sigma
	cov.unscaled <- out@cov.unscaled
	df.den <- out@df.residual
	
	tcrit <- qt((1+conf)/2, df.den)
	
	if(!is.null(out@mdc))
		mdcl <- out@mdc
	else 
		mdcl <- mdc(out, m = m)
	
	xmax <- max(out@x, na.rm = T)
	if(is.finite(mdcl)){
		xmin <- mdcl
		b0 <- out@coefficients[1]
		b1 <- out@coefficients[2]
		if(out@type == "lin"){
			svec <- c(cov.unscaled[, 1], cov.unscaled[, 2])
			mdc.y <- b0 + b1 * mdcl
		}
		if(out@type == "quad"){
			b2 <- out@coefficients[3]
			svec <- c(cov.unscaled[, 1], cov.unscaled[, 2], cov.unscaled[, 3])
			mdc.y <- b0 + b1 * mdcl + b2 * mdcl^2
		}
		p <- length(out@coefficients)
		pred.limit <- function(out, xp, m, tcrit, upper = T)
		{
			if(missing(m))
				m <- out@m
			if(out@type == "lin") {
				xmat <- cbind(1, xp)
				yp <- xmat %*% out@coefficients
#				yp <- out@coefficients[1] + xp * out@coefficients[2]
			}
			if(out@type == "quad") {
				xmat <- cbind(1, xp, xp^2)
				yp <- xmat %*% out@coefficients
#				yp <- out@coefficients[1] + xp * out@coefficients[2] + 
#						xp^2 * out@coefficients[3]
			}
			qn.term1 <- ((yp^2)^theta)/m
			qn.term2 <- xmat %*% cov.unscaled %*% t(xmat)
			qn.term2 <- diag(qn.term2)
			qn.p <- out@sigma * sqrt(qn.term1 + qn.term2)
			if(upper)
				cl <- yp + tcrit * qn.p
			else 
				cl <- yp - tcrit * qn.p
			cl
		}
		
		xp <- exp(seq(log(xmin), log(max(out@x)), 
						length = 500))
		if(b1 > 0) {
			y.mdclim <- pred.limit(out, 0, tcrit = tcrit, upper = T)
			ylim <- pred.limit(out, xp, tcrit = tcrit, upper = F)
		}
		else {
			y.mdclim <- pred.limit(out, 0, tcrit = tcrit, upper = F)
			ylim <- pred.limit(out, xp, tcrit = tcrit, upper = T)
		}
		
		rdl.c <- approx(ylim, xp, y.mdclim, ties = "mean")$y
		
	}
	else 
		rdl.c <- NA
	
	## 07-01-08: I'm not really sure why rdl returns a data frame. I don't think
	## that this is necessary. Having rdl return rdl.c from above. DVS
#	rdlout <- list(rdl = as.vector(rdl.c), m = m)
#	return(rdlout)
	return(rdl.c)
}

## 11-17-08: Have a method for mdc, do not need this anymore. DVS
#mdc.calc <- function(out, conf = out@conf.level, m = out@m){
#loq.calc <- function(out, m = out@m, cv = out@cv, vlen = 700, mit = 10000, toler = .001)
###################################################################
### This computes the Limits of Quantitation                     ##
### UseMethod means if we do loq(object) that the function       ##
###  will do loq.objectclass(object).  For example, if the       ##
###  object class is fpl.pom loq(object)                         ##
###     is same as loq.pom.fpl(object).                          ##
###                                                              ##
###################################################################
#{
#	if((out@type == "lin") | (out@type == "quad")){
#		loq.out <- loq.lin.pom(out, m = m, cv = cv,
#				vlen = vlen)
#	}
#	if(out@type == "thpl"){
#		loq.out <- loq.thpl.pom(out, m = m, cv = cv,
#				vlen = vlen, mit = mit, toler = toler)
#	}
#	if(out@type == "fpl"){
#		loq.out <- loq.fpl.pom(out, m = m, cv = cv,
#				vlen = vlen, mit = mit, toler = toler)
#	}
#	if(out@type == "tpl"){
#		loq.out <- loq.tpl.pom(out, m = m, cv = cv,
#				vlen = vlen, mit = mit, toler = toler)
#	}
#	
#	return(loq.out)
#}

loq.fpl.pom <- function(out, m = out@m, cv = out@cv, vlen = 700, mit = 10000, toler = 0.001)
######################################################################
## This function computes the Limits of Quantitation                ##
##  for FPL POM model.                                              ##
##                                                                  ##
## input variables:                                                 ##
##   fpl.out: the output from the fitting function                  ##
##   m: number of reps (number of y-values for each x)              ##
##   cv: the minimum coefficient of variation desired               ##
##   mit, toler: used if dos=F                                      ##
##   vlen: used in computing a range of x below                     ##
##                                                                  ##
## comments added by Matthew Mitchell 3/12/03  have ##              ##
## the original comments have ###                                   ##
######################################################################
{
	cov.un <- out@cov.unscaled
	b <- out@coefficients
	x <- out@x
	n <- length(x)
	b50 <- ifelse(out@logParm, exp(b[3]), b[3])
	xpstart <- min(c(5e-4, min(x[x > 0])))
	xp <- c(seq(xpstart, b50, length = round(vlen/2, 0)), seq(b50, 
					max(x), length = round(vlen/2, 0)))
	yp <- as.vector(fpl.model(xp, b, logParm = out@logParm))
	## Getting the gradient of x as a function of y
	if(!out@logParm){
		dh.dy <- xp * ((b[2] - b[1])/(b[4] * (b[1] - yp) * (yp - 
							b[2])))
		dh.db1 <- xp/(b[4] * (b[1] - yp))
		dh.db2 <- xp/(b[4] * (yp - b[2]))
		dh.db3 <- xp/b[3]
		dh.db4 <- ( - xp/(b[4] * b[4])) * log((b[1] - yp)/(yp - 
							b[2]))
	}
	else {
		dh.dy <- (xp * (b[2] - b[1]))/(b[4] * (b[1] - yp) * (yp -
						b[2]))
		dh.db1 <- xp/(b[4] * (b[1] - yp))
		dh.db2 <- xp/(b[4] * (yp - b[2]))
		dh.db3 <- xp
		dh.db4 <- ( - xp/(b[4] * b[4])) * log((b[1] - yp)/(yp - 
							b[2]))
	}
	sigma2 <- out@sigma^2
	## Approximating variance of x using Wald's method.
	dh.db <- cbind(dh.db1, dh.db2, dh.db3, dh.db4)
	var.xnot.hat <- sigma2 * (apply(dh.db, 1, function(z) (t(z) %*% cov.un %*% z)) + 
				((dh.dy^2) * ((yp^2)^out@theta))/m)
	sd <- sqrt(var.xnot.hat)
	xp <- xp[is.finite(sd)]
	sd <- sd[is.finite(sd)]
	xmin <- xp[sd/xp == min(sd/xp)]
	
	if(cv > min(sd/xp[xp != 0]))
		## 01-25-07 Added ties="mean" to approx to avoid error message generated
		## from ties. DVS
		loq.out <- approx(((sd/xp) - cv)[xp < xmin], 
				xp[xp < xmin], 0,ties="mean")$y
	else 
		loq.out <- NA
	## 07-28-08: Not sure why this is returning a data frame. m and cv are contined in the calib.fit object. DVS
#	loq.out <- data.frame(loq = loq.df, m = m, cv = cv)
#	loq.out
	return(loq.out)
}

loq.thpl.pom <- function(out, m = out@m, cv = out@cv, vlen = 700, mit = 
				10000, toler = 0.001)
{
	######################################################################
	## This function computes the Limits of Quantitation                ##
	##  for THPL POM model.                                             ##
	##                                                                  ##
	## input variables:                                                 ##
	##   out: the output from the fitting function                 		##
	##   m: number of reps (number of y-values for each x)              ##
	##   cv: the minimum coefficient of variation desired               ##
	##   mit, toler: used if dos=F                                      ##
	##   vlen: used in computing a range of x below                     ##
	##                                                                  ##
	## comments added by Matthew Mitchell 3/12/03  have ##              ##
	## the original comments have ###                                   ##
	######################################################################
	
	cov.un <- out@cov.unscaled
	b <- out@coefficients
	x <- out@x
	n <- length(x)
#	b50 <-  b[3]
	b50 <- ifelse(out@logParm, exp(b[3]), b[3])
	xpstart <- min(c(5e-4, min(x[x > 0])))
	xp <- c(seq(xpstart, b50, length = round(vlen/2, 0)), seq(b50, 
					max(x), length = round(vlen/2, 0)))
	yp <- as.vector(thpl.model(xp, b))
	## Getting the gradient of x as a function of y
	if(!out@logParm) {
		dh.dy <- xp * ((b[2] - b[1])/((b[1] - yp) * (yp - 
							b[2])))
		dh.db1 <- xp/((b[1] - yp))
		dh.db2 <- xp/((yp - b[2]))
		dh.db3 <- xp/b[3]
	}
	else {
		dh.dy <- (xp * (b[2] - b[1]))/((b[1] - yp) * (yp -
						b[2]))
		dh.db1 <- xp/((b[1] - yp))
		dh.db2 <- xp/((yp - b[2]))
		dh.db3 <- xp
	}
	sigma2 <- out@sigma^2
	## Approximating variance of x using Wald's method.
	dh.db <- cbind(dh.db1, dh.db2, dh.db3)
	var.xnot.hat <- sigma2 * (apply(dh.db, 1, function(z) (t(z) %*% cov.un %*% z)) + 
				((dh.dy^2) * ((yp^2)^out@theta))/m)
#	var.xnot.hat <- (((dh.dy * dh.dy) * sigma2 * ((yp^2)^out@
#						theta))/m + sigma2 * (dh.db1 * (dh.db1 * cov.un[1, 1] + 
#						dh.db2 * cov.un[2, 1] + dh.db3 * cov.un[3, 1]) +
#					dh.db2 * (dh.db1 * cov.un[1, 2] + 
#						dh.db2 * cov.un[2, 2] + dh.db3 * cov.un[3, 2])
#					+ dh.db3 * (dh.db1 * cov.un[1, 3] + 
#						dh.db2 * cov.un[2, 3] + dh.db3 * cov.un[3, 3])))
	sd <- sqrt(var.xnot.hat)
	xp <- xp[is.finite(sd)]
	sd <- sd[is.finite(sd)]
	xmin <- xp[sd/xp == min(sd/xp)]
	#if(dos) {
	if(cv > min(sd/xp[xp != 0]))
		## 01-25-07 Added ties="mean" to approx to avoid error message generated when
		## there are ties. DVS
		loq.out <- approx(((sd/xp) - cv)[xp < xmin], 
				xp[xp < xmin], 0, ties="mean")$y
	else 
		loq.out <- NA
	
#	loq.out <- data.frame(loq = loq.df, m = m, cv = cv)
#	loq.out
	return(loq.out)
}


loq.tpl.pom <- function(out, m = out@m, cv = out@cv, vlen = 700, 
		mit = 10000, toler = 0.001)
{
	######################################################################
	## This function computes the Limits of Quantitation                ##
	##  for FPL POM model.                                              ##
	##                                                                  ##
	## input variables:                                                 ##
	##   fpl.out: the output from the fitting function                  ##
	##   m: number of reps (number of y-values for each x)              ##
	##   cv: the minimum coefficient of variation desired               ##
	##   mit, toler: used if dos=F                                      ##
	##   vlen: used in computing a range of x below                     ##
	##   dos: is T or F, used to change the data type                   ##
	##                                                                  ##
	## comments added by Matthew Mitchell 3/12/03  have ##              ##
	## the original comments have ###                                   ##
	######################################################################
	
	## Collecting variable information
	cov.un <- out@cov.unscaled
	b <- out@coefficients
	x <- out@x
	n <- length(x)
	b50 <- ifelse(out@logParm, exp(b[3]), b[3])
	sigma2 <- out@sigma^2
	
	xpstart <- min(c(5e-4, min(x[x > 0])))
	xp <- c(seq(xpstart, b50, length = round(vlen/2, 0)), 
			seq(b50, max(x), length = round(vlen/2, 0)))
	yp <- as.vector(tpl.model(xp, b, logParm = out@logParm))
	## Getting the gradient of x as a function of y
	if(out@logParm) {
		dh.dy <- xp * ((b[2] - b[1])/(b[4] * (b[1] - yp) * (yp - b[2])))
		dh.db3 <- xp/b[3]
		dh.db4 <- ( - xp/(b[4] * b[4])) * log((b[1] - yp)/(yp - b[2]))
	}
	else {
		dh.dy <- (xp * (b[2] - b[1]))/(b[4] * (b[1] - yp) * (yp - b[2]))
		dh.db3 <- xp
		dh.db4 <- ( - xp/(b[4] * b[4])) * log((b[1] - yp)/(yp - b[2]))
	}
	## Approximating variance of x using Wald's method.
	dh.db <- cbind(dh.db3, dh.db4)
	var.xnot.hat <- sigma2 * (apply(dh.db, 1, function(z) (t(z) %*% cov.un %*% z)) + 
				((dh.dy^2) * ((yp^2)^out@theta))/m)

	sd <- sqrt(var.xnot.hat)
	xp <- xp[is.finite(sd)]
	sd <- sd[is.finite(sd)]
	xmin <- xp[sd/xp == min(sd/xp)]
	
	if(cv > min(sd/xp[xp != 0]))
		## 01-25-07 Added ties="mean" to approx to avoid error message generated
		## from ties. DVS
		loq.out <- approx(((sd/xp) - cv)[xp < xmin], xp[xp < xmin], 0, ties = "mean")$y
	else 
		loq.out <- NA
	
	return(loq.out)
}




## The derivatives for the quadratic are incorrect - corrections made below by MM 4/15/03
## The mistakes mainly resulted from switching a and c so b[1] and b[3] are switched
## The model is a + bx +cx^2 and a=b[1], b=b[2], c=b[3] and dx/dy reciprocated one part
## incorrectly.  The previous versions have been commented out
## The expressions have simpler forms also when using such facts as
## y=a + bx + cx^2 so y-a = bx + cx^2 and sqrt(b2 - 4ac + 4cy) = b + 2cx

loq.lin.pom <- function(out, m = out@m, cv = out@cv, vlen = 500)
{
	theta <- out@theta
	b <- out@coefficients
	x <- out@x
	cov.un <- out@cov.unscaled
	sigma <- out@sigma
	xpstart <- min(c(5e-4, min(out@x[out@x > 0])))
	if(missing(m))
		m <- out@m
	if(missing(cv)) 
		cv <- out@cv	
	## set up a grid of x values to determine the upper boundary for the search 
	##   algorithm, this is the value where CV is at a minimum (called xmin)
	xp <- exp(seq(log(xpstart), log(max(x)), length = vlen))
	if(out@type == "lin") {
		yp <- b[1] + b[2] * xp
		var.xnot.hat <- (sigma^2) * ((yp^2)^theta)/((b[2]^2) * m) + 
				(cov.un[1, 1] + 2 * cov.un[1,2] * xp + cov.un[2, 2] * xp * 
					xp)/(b[2]^2)
	}
	else {
		yp <- b[1] + b[2] * xp + b[3] * xp * xp
		
		dh.dy <- 1/(b[2] + 2*b[3]*xp)
		
		dh.da <- -1/(b[2] + 2*b[3]*xp)
		
		dh.db <- -xp/(b[2] + 2*b[3]*xp)
		
		dh.dc <- -(xp*xp)/(b[2] + 2*b[3]*xp)
		var.xnot.hat <- ((dh.dy * dh.dy) * sigma * sigma * ((yp^2)^
						out@theta))/m + sigma * sigma * (dh.da *
					(dh.da * cov.un[1, 1] + dh.db * cov.un[1, 2] + 
						dh.dc * cov.un[1, 3]) + dh.db * (dh.da * cov.un[
								1, 2] + dh.db * cov.un[2, 2] + dh.dc * cov.un[2,
								3]) + dh.dc * (dh.da * cov.un[1, 3] + dh.db * 
						cov.un[2, 3] + dh.dc * cov.un[3, 3]))
	}
	sd <- sqrt(var.xnot.hat)
	xmin <- xp[sd/xp == min(sd/xp)]
	xl <- min(xp)
	xu <- xmin	### perform the search for the LOQ
	svec <- as.vector(cov.un)
	if(cv > min(sd/xp[xp != 0])) {
		## 01-25-07 Added ties="mean" to approx function to avoid error messages
		## generated when there are ties. DVS
#		loq.df <- approx(((sd/xp) - cv)[xp < xmin], xp[xp < xmin], 
#				0, ties = "mean")$y
#		loq.out <- list(loq = loq.df, m = m, cv = cv)
		loq.out <- approx(((sd/xp) - cv)[xp < xmin], xp[xp < xmin], 
				0, ties = "mean")$y
	}
	else 
		loq.out <- NA
#		loq.out <- list(loq = NA, m = m, cv = cv)
	
	return(loq.out)
}

lof.test <- function(out)
{
	########################################################
	## This computes the lack of fit statistics from the  ##
	##  nls POM fit.                                      ##
	## This function gets called as part of the fpl.pom   ##
	##  function.                                         ##
	##                                                    ##
	## comments added by Matthew Mitchell on 3/12/03      ##
	########################################################
	
	x <- out@x[is.finite(out@y)]
	y <- out@y[is.finite(out@y)]
	b <- out@coefficients
	nx <- tapply(x, list(x), length)
	vars <- tapply(y, list(x), var)[nx >= 2]
	wts <- tapply(((out@fitted.values^2)^out@theta), 
			list(out@x), mean)[nx >= 2]
	nx <- nx[nx >= 2]
	pure <- (((nx - 1) * vars)/wts)
	pure <- sum(pure)
	df.pure <- sum(nx - 1)
	sse <- sum((out@residuals^2)/
					((out@fitted.values^2)^out@theta))
	df.den <- out@df.residual
	lofss <- sse - pure
	df.lof <- df.den - df.pure
	lof.stat <- (lofss/df.lof)/(pure/df.pure)
	pv <- 1 - pf(lof.stat, df.lof, df.pure)
	## 07-28-08: I don't think that this should be its own class. It would be more appropriate to have this as a 
	## data frame object. DVS
#	lof.out <- new("lof.test", Fstat = lof.stat, p.value = 
#					pv, lofss = lofss, df.lof = df.lof, pure.error = pure, 
#			df.pure.error = df.pure, sse = sse, df.sse = df.den)
	lof.df <- data.frame(Fstat = lof.stat, p.value = pv, lofss = lofss,
			df.lof = df.lof, pure.error = pure, df.pure.error = df.pure, sse = sse,
			df.sse = df.den)
	return(lof.df)
}

