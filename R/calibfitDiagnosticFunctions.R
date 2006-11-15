## calibfit-diagnostic-functions.R
## 
## Functions that are called by other calibfit functions. 
## These are not accessible to the user.
##
## Author: samarov
###############################################################################

mdc.calc <- function(calib.fit.out,conf=0.95,mit=1000,
	toler=1e-04,type=calib.fit.out@type){
	if (type=="lin"){
    	mdc.out <- mdc.lin.pom(calib.fit.out,m=calib.fit.out@m,
	    	conf=conf)#, dos=calib.fit.out@dos)
	}
	else {if(type=="thpl"){
    	mdc.out <- mdc.thpl.pom(calib.fit.out,m=calib.fit.out@m,
			xu=max(calib.fit.out@x),conf=conf,mit=mit,toler=toler)#,
	        #dos=calib.fit.out@dos)
	}
	else {if(type=="fpl"){
    	mdc.out <- mdc.fpl.pom(calib.fit.out,m=calib.fit.out@m,
			conf=conf)#, dos=calib.fit.out@dos)
		}
	}
	}
}

mdc.fpl.pom <- function(fpl.out, m = fpl.out@m, 
	xu = max(fpl.out@x), conf = 0.95, mit
	= 1000, toler = 0.0001)#, dos = fpl.out@dos)
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
##   mit, toler: used if dos=F                                      ##
##   dos: is T or F, used to change the data type                   ##
##                                                                  ##
## old comments have ###                                            ##
## new version created by MM on 4/14/03                             ##
######################################################################
  
	## reads in the necessary parameters
	cov.un <- fpl.out@cov.unscaled
	b <- fpl.out@coefficients
    s <- fpl.out@sigma

	## mdc is computed by finding the x that corresponds to the upper
	## prediction limit of y when x=0 (for increasing curve - b1 < b2)
	## or the lower prediction limit of y when x=0 (decreasing curve b1 > b2)
	## The predicted value of y when x=0 is b1. So we find the prediction interval for b1.

	f0b <- b[1]
	var.f0b <-  s^2*(((f0b^2)^fpl.out@theta)/m + (cov.un[1, 1]))
	## 04-03-07: Changed from below to above
	#var.f0b <- s*s*((f0b^(2*fpl.out@theta))/m + (cov.un[1, 1]))
    tcrit <- qt( (1+conf)/2, fpl.out@df.residual)
    lim0 <- f0b + (sign(b[2]-b[1]))*tcrit*sqrt(var.f0b)
    mdcval <- fpl.inverse(fpl.out,lim0)
       
	### if not on a PC, call C function to search for MDC
    ## This part probably needs to be modified if this program is not
    ## going to be run on a PC
	#	if(!dos) {
	#		if(!is.loaded(symbol.C("fplpommdc")))
	#			dyn.load("/local/lib/calib/sun4/src/fplpom.o")
	#		mdc.out <- .C("fplpommdc",
	#			xl = as.double(0),
	#			xu = as.double(xu),
	#			as.double(b),
	#			as.integer(mit),
	#			as.double(toler),
	#			as.double(fpl.out@theta),
	#			as.double(tcrit),
	#			as.double(fpl.out@sigma),
	#                       as.double(f0b),
	#			as.double(var.f0b),
	#			as.integer(m),
	#			as.integer(fpl.out@parm))
	#		mdcval <- mdc.out$xl
	#	}
        
    mdc.df <- data.frame(mdc=mdcval, m=m)    
}

mdc.thpl.pom <- function(thpl.out, m = thpl.out@m, xu = max(thpl.out@x), conf = 0.95, mit
	 = 1000, toler = 0.0001)#, dos = thpl.out@dos)
{
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
##   dos: is T or F, used to change the data type                   ##
##                                                                  ##
## written by Matthew Mitchell on 5/20/03                           ##
######################################################################
  
	## reads in the necessary parameters
	cov.un <- thpl.out@cov.unscaled
	b <- thpl.out@coefficients
    s <- thpl.out@sigma

	## mdc is computed by finding the x that corresponds to the upper
	## prediction limit of y when x=0 (for increasing curve - b1 < b2)
	## or the lower prediction limit of y when x=0 (decreasing curve b1 > b2)
	## The predicted value of y when x=0 is b1. So we find the prediction interval for b1.

	f0b <- b[1]
	var.f0b <-  s^2*(((f0b^2)^thpl.out@theta)/m + (cov.un[1, 1]))
	## 04-03-07: Changed from below to above
	#var.f0b <- s*s*((f0b^(2*thpl.out@theta))/m + (cov.un[1, 1]))
	tcrit <- qt( (1+conf)/2, thpl.out@df.residual)
	lim0 <- f0b + (sign(b[2]-b[1]))*tcrit*sqrt(var.f0b)
	mdcval <- thpl.inverse(thpl.out,lim0)
       
	### if not on a PC, call C function to search for MDC
        ## This part probably needs to be modified if this program is not
        ## going to be run on a PC
	## 02-08-07: This call to C does not seem relavent. DVS
	#	if(!dos) {
	#		if(!is.loaded(symbol.C("thplpommdc")))
	#			dyn.load("/local/lib/calib/sun4/src/thplpom.o")
	#		mdc.out <- .C("thplpommdc",
	#			xl = as.double(0),
	#			xu = as.double(xu),
	#			as.double(b),
	#			as.integer(mit),
	#			as.double(toler),
	#			as.double(thpl.out@theta),
	#			as.double(tcrit),
	#			as.double(thpl.out@sigma),
	#                       as.double(f0b),
	#			as.double(var.f0b),
	#			as.integer(m))
	#		mdcval <- mdc.out$xl
	#	}
        
    mdc.df <- data.frame(mdc=mdcval, m=m)    
}

mdc.lin.pom <- function(lin.out, m = lin.out@m, conf = 0.95)#, 
						#dos = lin.out@dos)
{
######################################################################
## This function computes the Mimimum Detectable Concentration      ##
##  for the Linear POM model.                                       ##
##                                                                  ##
## input variables:                                                 ##
##   lin.out: the output from the fitting function                  ##
##   m: number of reps (number of y-values for each x)              ##
##   conf: confidence level (in point-wise intervals)               ##
##   dos: is T or F, used to change the data type                   ##
##                                                                  ##
## old comments have ###                                            ##
## new version created by MM on 4/14/03                             ##
######################################################################
	min.x <- min(lin.out@x)
	max.x <- max(lin.out@x)
	b <- lin.out@coefficients
	cov.un <- lin.out@cov.unscaled
	s <- lin.out@sigma
	df.den <- lin.out@df.residual
	tcrit <- qt( (1+conf)/2, df.den)      
	xu <- max(lin.out@x)
	f0b <- lin.out@coefficients[1]
	var.f0b <-  s*s*(((f0b^2)^lin.out@theta)/m + (cov.un[1, 1]))
	## 04-03-07: Changed from below to above
	#var.f0b <- s*s*((f0b^(2*lin.out@theta))/m + (cov.un[1, 1]))
     
	if(lin.out@mmod == "lin"){
		lim0 <- f0b + (sign(b[2]))*tcrit*sqrt(var.f0b)
        mdcval <- (lim0 - b[1])/b[2]
        }

        ## The only valid quadratic models are those whose critical point
        ## is less than 0 or higher than max(x)
        ## Whether we use the upper or lower interval is 0 is determined
        ## by whether the curve is increasing or decreasing - it won't be
        ## both because the critical point outside the range.  Hence
        ## it is increasing if f(1) > f(0) which means b + 2c > 0
         
	if(lin.out@mmod == "quad"){
		lim0 <- f0b + (sign(b[2]+ 2*b[3]))*tcrit*sqrt(var.f0b)
		mdc1 <- (-b[2] - sqrt(b[2]^2 - 4*b[3]*(b[1]-lim0)))/(2*b[3])
		mdc2 <- (-b[2] + sqrt(b[2]^2 - 4*b[3]*(b[1]-lim0)))/(2*b[3])
		mdcval <- ifelse(mdc1 > min.x & mdc1 < max.x, mdc1, mdc2)
	    }
	list(mdc = as.vector(mdcval), m = m)  ## I left this as it was - not sure why
        ## it is a list here while the other was made into a df
}

rdl.calc <- function(calib.fit.out,conf=.95,tol=1e-04,xu,mit=1000,type="fpl"){
	if(type=="lin"){
		rdl.out <- rdl.lin.pom(calib.fit.out, m = calib.fit.out@m,
			conf = conf, tol = tol)#, dos = calib.fit.out@dos)
	}
	else {if(type=="thpl"){
    	rdl.out <- rdl.thpl.pom(calib.fit.out, m=calib.fit.out@m, xu, conf = conf,
		    mit = mit, toler = tol)#, dos = calib.fit.out@dos)
	}
	else {if(type=="fpl"){
    	rdl.out <- rdl.fpl.pom(calib.fit.out, m=calib.fit.out@m, xu, conf = conf,
		    mit = mit, toler = tol)#, dos = calib.fit.out@dos)
	}
    }
    }
}

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
##   dos: is T or F, used to change the data type                   ##
##                                                                  ##
## comments added by Matthew Mitchell 3/12/03  have ##              ##
## the original comments have ###                                   ##
######################################################################

rdl.fpl.pom <- function(fpl.out, m, xu, conf = 0.95, mit = 1000, toler = 0.0001)
#, dos = fpl.out@dos) # No longer using dos
{
### this function takes output from a regression 
### and  returns a value for Reliable Detection Limit
	cov.un <- fpl.out@cov.unscaled
	b <- fpl.out@coefficients
	n <- length(fpl.out@x)
	if(missing(m))
		m <- fpl.out@m
	## 03-12-07: No longer needs to be calculated. DVS
	#bigs.vec <- c(cov.un[, 1], cov.un[, 2], cov.un[, 3], cov.un[, 4
	#	])
	f0b <- b[1]
	var.f0b <- (((f0b^2)^fpl.out@theta) * fpl.out@sigma * fpl.out@sigma)/m + 
		(fpl.out@sigma * fpl.out@sigma) * cov.un[1,1]
## The computation of the critical value is incorrect.  It should be
## 1 - (1-conf)/2 = (1+conf)/2.  Correction made by MM on 4/9/03.

    tcrit <- qt((1+conf)/2, fpl.out@df.residual)
	lim0 <- f0b + (sign(b[2] - b[1])) * tcrit * sqrt(var.f0b)
	xl <- fpl.inverse(fpl.out, lim0)
#	if(!dos) {
#		if(!is.loaded(symbol.C("fplpomcal")))
#			dyn.load("/local/lib/calib/sun4/src/fplpom.o")
#		if(missing(xu))
#			xu <- max(fpl.out$x)
#		rdl.out <- .C("fplpomcal",
#			xl = as.double(xl),
#			xu = as.double(xu),
#			as.double(b),
#			as.double(bigs.vec),
#			as.integer(mit),
#			as.double(toler),
#			as.double(fpl.out@theta),
#			as.double(tcrit),
#			as.double(fpl.out@sigma),
#			as.double(lim0),
#			as.integer(sign(b[2] - b[1])),
#			as.integer(1),
#			as.integer(m),
#			as.integer(fpl.out@parm),
#			as.double(as.integer(sign(b[1] - b[2]))))
#		rdlval <- rdl.out@xl
#	}
#	else {
	b50 <- ifelse(fpl.out@parm == 1, b[3], exp(b[3]))
	xp <- seq(xl, max(fpl.out@x), len = 450)
		## 04-17-07: Seeing if this correction will do a better job
		## estimating RDL. Reason for change is that b50 is at times
		## out of the range of concentrations. This seems to generally be
		## happening when attempting to fit fpl to data where a more
		## appropriate fit might be a linear or quadrattic fit. DVS
#		xp <- c(seq(xl, b50, len = 350), 
#				seq(b50, max(fpl.out@x), len = 100))
	ypred <- fpl.model(xp, b, parm = fpl.out@parm)
	qn1 <- ((as.vector(ypred)^2)^fpl.out@theta)/m
		## 04-03-07: Changing to above DVS
		#qn1 <- ((as.vector(ypred))^(2 * thpl.out@theta))/m
	qn2 <- diag(attributes(ypred)$gradient %*% 
		cov.un %*% t(attributes(ypred)$gradient))
	yp <- as.vector(ypred) - sign(b[2] - b[1]) * tcrit * 
		fpl.out@sigma * sqrt(qn1 + qn2)
		## 01-25-07 Added ties="mean" to approx function below to avoid error messages
		## when there are ties. DVS
	rdlval <- approx(yp, xp, lim0, ties="mean")$y
#	}
	rdlout <- data.frame(rdl = rdlval, m = m)
	rdlout
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
##   dos: is T or F, used to change the data type                   ##
##                                                                  ##
## written by Matthew W. Mitchell on 5/20/03                        ##
## the original comments have ###                                   ##
######################################################################

rdl.thpl.pom <- function(thpl.out, m, xu, conf = 0.95, mit = 1000, toler = 0.0001)
#, dos = thpl.out@dos)
{
### this function takes output from a regression 
### and  returns a value for Reliable Detection Limit
	cov.un <- thpl.out@cov.unscaled
	b <- thpl.out@coefficients
	n <- length(thpl.out@x)
	if(missing(m))
		m <- thpl.out@m
	bigs.vec <- c(cov.un[, 1], cov.un[, 2], cov.un[, 3])
	f0b <- b[1]
	var.f0b <- (((f0b^2)^thpl.out@theta) * thpl.out@sigma^2)/m + 
		(thpl.out@sigma^2) * cov.un[1,1]
## The computation of the critical value is incorrect.  It should be
## 1 - (1-conf)/2 = (1+conf)/2.  Correction made by MM on 4/9/03.
##	tcrit <- qt(conf, thpl.out$df.res)
    tcrit <- qt((1+conf)/2, thpl.out@df.residual)
	lim0 <- f0b + (sign(b[2] - b[1])) * tcrit * sqrt(var.f0b)
	xl <- thpl.inverse(thpl.out, lim0)
#	if(!dos) {
#		if(!is.loaded(symbol.C("thplpomcal")))
#			dyn.load("/local/lib/calib/sun4/src/thplpom.o")
#		if(missing(xu))
#			xu <- max(thpl.out@x)
#		rdl.out <- .C("thplpomcal",
#			xl = as.double(xl),
#			xu = as.double(xu),
#			as.double(b),
#			as.double(bigs.vec),
#			as.integer(mit),
#			as.double(toler),
#			as.double(thpl.out@theta),
#			as.double(tcrit),
#			as.double(thpl.out@sigma),
#			as.double(lim0),
#			as.integer(sign(b[2] - b[1])),
#			as.integer(1),
#			as.integer(m),
#		        as.double(as.integer(sign(b[1] - b[2]))))
#		rdlval <- rdl.out$xl
#	}
#	else {
	b50 <- ifelse(thpl.out@parm, b[3], exp(b[3]))
	xp <- c(seq(xl, b50, len = 350), 
		seq(b50, max(thpl.out@x), len = 100))
	ypred <- thpl.model(xp, b, parm = thpl.out@parm)
	qn1 <- ((as.vector(ypred)^2)^thpl.out@theta)/m
	## 04-03-07: Changing to above DVS
	#qn1 <- ((as.vector(ypred))^(2 * thpl.out@theta))/m
	qn2 <- diag(attributes(ypred)$gradient %*% cov.un %*% t(
		attributes(ypred)$gradient))
	yp <- as.vector(ypred) - sign(b[2] - b[1]) * tcrit * 
		thpl.out@sigma * sqrt(qn1 + qn2)
	## 01-25-07 Added ties="mean" to approx function below to avoid error messages
	## when there are ties. DVS
	rdlval <- approx(yp, xp, lim0, ties="mean")$y
#	}
	rdlout <- data.frame(rdl = rdlval, m = m)
	rdlout
}

rdl.lin.pom <- function(lin.out, m = lin.out@m, conf = 0.95, tol = 0.01)
#, dos = lin.out@dos)
{
####	Calculates limit of quantitation (RDL)
###	with a data set and the output from glsvfe.lin.out()
###
###	RDL = Reliable Detection Limit
###		= lowest concentration that can reliabily be expected
###		to give a reading above the minimum detectable concentration
###		- i.e., lower conf. limit at x=RDL > upper conf. limit at x=0
###	m is # reps of the new unknown
###	correct df=df.den can be used for common theta/sigma cases
###		else if df.den is missing, it is gotten from lin.out
###	conf=.95 is the default that goes along with 90% two-sided
###		confidence limits about the curve
###	the default is to get mdc from pom, it can also be given explicitly

## 02-08-07: This call to C does not seem relavent. DVS
#	if(!is.loaded(symbol.C("linlow")) & dos == F) {
#		cat("Loading linpom.o\n")
#		dyn.load("/local/lib/calib/sun4/src/linpom.o")
#	}
	theta <- lin.out@theta
	std.dev <- lin.out@sigma
	sg.inv <- lin.out@cov.unscaled
	df.den <- lin.out@df.residual
## line below is incorrect. Corrected by MM on 4/14/03
##      tcrit <- qt(conf, df.den)
	tcrit <- qt( (1+conf)/2, df.den)
    if(m == lin.out@m)
	if(!is.null(lin.out@mdc))
		mdcl <- lin.out@mdc
	else mdcl <- mdc(lin.out)
	else mdcl <- mdc(lin.out, m = m)
	xmax <- max(lin.out@x, na.rm = T)
	if(is.finite(mdcl)){
		xmin <- mdcl
		b0 <- lin.out@coefficients[1]
		b1 <- lin.out@coefficients[2]
		if(lin.out@mmod == "lin"){
			svec <- c(sg.inv[, 1], sg.inv[, 2])
			mdc.y <- b0 + b1 * mdcl
		}
		if(lin.out@mmod == "quad"){
			b2 <- lin.out@coefficients[3]
			svec <- c(sg.inv[, 1], sg.inv[, 2], sg.inv[, 3])
			mdc.y <- b0 + b1 * mdcl + b2 * mdcl^2
		}
		p <- length(lin.out@coefficients)
		pred.limit <- function(pom, xp, m, tcrit, upper = T)
		{
			if(missing(m))
				m <- pom@m
			if(pom@mmod == "lin") {
				yp <- pom@coefficients[1] + xp * pom@coefficients[2]
				xmat <- cbind(1, xp)
			}
			if(attr(pom, "mmod") == "quad") {
				yp <- pom@coefficients[1] + xp * pom@coefficients[2] + 
				  xp * xp * pom@coefficients[3]
				xmat <- cbind(1, xp, xp * xp)
			}
			qn.term1 <- ((yp^2)^pom@theta)/m
			qn.term2 <- xmat %*% pom@cov.unscaled %*% t(xmat)
			qn.term2 <- diag(qn.term2)
			qn.p <- pom@sigma * sqrt(qn.term1 + qn.term2)
			if(upper)
				cl <- yp + tcrit * qn.p
			else cl <- yp - tcrit * qn.p
			cl
		}
#		if(!dos) {
#			if(b1 > 0) {
#				lout <- .C("linlow",
#				  xl = as.double(xmin),
#				  xu = as.double(xmax),
#				  as.double(lin.out@coefficients),
#				  as.double(svec),
#				  as.integer(10000),
#				  as.double(0.001),
#				  as.double(theta),
#				  as.double(tcrit),
#				  as.double(std.dev),
#				  as.double(mdc.y),
#				  as.integer(1),
#				  as.integer(m),
#				  as.integer(p))
#				rdl.c <- lout$xl
#			}
#			if(b1 < 0) {
#				lout <- .C("linup",
#				  xl = as.double(xmin),
#				  xu = as.double(xmax),
#				  as.double(lin.out@coef),
#				  as.double(svec),
#				  as.integer(10000),
#				  as.double(0.001),
#				  as.double(theta),
#				  as.double(tcrit),
#				  as.double(std.dev),
#				  as.double(mdc.y),
#				  as.integer(1),
#				  as.integer(m),
#				  as.integer(p))
#				rdl.c <- lout$xu
#			}
#		}
#		else {
		xp <- exp(seq(log(xmin), log(max(lin.out@x)), 
			length = 500))
		if(b1 > 0) {
			y.mdclim <- pred.limit(lin.out, 0, 
				tcrit = tcrit, upper = T)
			ylim <- pred.limit(lin.out, xp, tcrit
				= tcrit, upper = F)
			}
		else {
			y.mdclim <- pred.limit(lin.out, 0, 
				tcrit = tcrit, upper = F)
			ylim <- pred.limit(lin.out, xp, tcrit
				= tcrit, upper = T)
			}
			## 01-25-07: Added ties="mean" for the interpolating function
			## approx, if this is not added than an error message is generated. DVS
			rdl.c <- approx(ylim, xp, y.mdclim, ties="mean")$y
#		}
	}
	else rdl.c <- NA
	list(rdl = as.vector(rdl.c), m = m)
}

loq.calc <- function(calib.fit.out,vlen=700,mit=10000,toler=.001,type="fpl")
##################################################################
## This computes the Limits of Quantitation                     ##
## UseMethod means if we do loq(object) that the function       ##
##  will do loq.objectclass(object).  For example, if the       ##
##  object class is fpl.pom loq(object)                         ##
##     is same as loq.pom.fpl(object).                          ##
##                                                              ##
##################################################################
{
	if(type=="lin"){
    	loq.out <- loq.lin.pom(calib.fit.out, m = calib.fit.out@m, cv = calib.fit.out@cv,
		    vlen = vlen)
	}
	else {if(type=="thpl"){
    	loq.out <- loq.thpl.pom(calib.fit.out, m = calib.fit.out@m, cv = calib.fit.out@cv,
			vlen = vlen, mit = mit, toler = toler)#, dos = calib.fit.out@dos)
	}
	else {if(type=="fpl"){
		loq.out <- loq.fpl.pom(calib.fit.out, m = calib.fit.out@m, cv = calib.fit.out@cv,
			vlen = vlen, mit = mit, toler = toler)#, dos = calib.fit.out@dos)
	}
    }
    }
}

loq.fpl.pom <- function(fpl.out, m = fpl.out@m, cv = fpl.out@cv, vlen = 700, mit = 
	10000, toler = 0.001)#, dos = fpl.out@dos)
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
  
	cov.un <- fpl.out@cov.unscaled
	b <- fpl.out@coefficients
	n <- length(fpl.out@x)
	b50 <- ifelse(fpl.out@parm == 1, b[3], exp(b[3]))
	xpstart <- min(c(0.0005, min(fpl.out@x[fpl.out@x > 0])))
	xp <- c(seq(xpstart, b50, length = round(vlen/2, 0)), seq(b50, 
		max(fpl.out@x), length = round(vlen/2, 0)))
	yp <- as.vector(fpl.model(xp, b, parm = fpl.out@parm))
	if(fpl.out@parm == 1) {
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
	sigma2 <- fpl.out@sigma^2 # * fpl.out@sigma
	var.xnot.hat <- (((dh.dy * dh.dy) * sigma2 * ((yp^2)^fpl.out@
		theta))/m + sigma2 * (dh.db1 * (dh.db1 * cov.un[1, 1] + 
		dh.db2 * cov.un[2, 1] + dh.db3 * cov.un[3, 1] + dh.db4 * 
		cov.un[4, 1]) + dh.db2 * (dh.db1 * cov.un[1, 2] + 
		dh.db2 * cov.un[2, 2] + dh.db3 * cov.un[3, 2] + dh.db4 * 
		cov.un[4, 2]) + dh.db3 * (dh.db1 * cov.un[1, 3] + 
		dh.db2 * cov.un[2, 3] + dh.db3 * cov.un[3, 3] + dh.db4 * 
		cov.un[4, 3]) + dh.db4 * (dh.db1 * cov.un[1, 4] + 
		dh.db2 * cov.un[2, 4] + dh.db3 * cov.un[3, 4] + dh.db4 * 
		cov.un[4, 4])))
	sd <- sqrt(var.xnot.hat)
	xp <- xp[is.finite(sd)]
	sd <- sd[is.finite(sd)]
	xmin <- xp[sd/xp == min(sd/xp)]
	#if(dos) {
	if(cv > min(sd/xp[xp != 0]))
	## 01-25-07 Added ties="mean" to approx to avoid error message generated
		## from ties. DVS
		loq.df <- approx(((sd/xp) - cv)[xp < xmin], 
			xp[xp < xmin], 0,ties="mean")$y
	else loq.df <- NA
	#}
#	else {
#		xl <- min(xp)
#		xu <- xmin
#		svec <- c(cov.un[, 1], cov.un[, 2], cov.un[, 3], cov.un[
#			, 4])
#		if(!is.loaded(symbol.C("fplpomprec")))
#			dyn.load("/local/lib/calib/sun4/src/fplpom.o")
#		if(cv > min(sd/xp[xp != 0])) {
#			bp.out <- .C("fplpomprec",
#				xl = as.double(xl),
#				xu = as.double(xu),
#				as.double(cv),
#				as.double(b),
#				as.double(svec),
#				as.integer(rep(mit, length(xl))),
#				as.double(toler),
#				as.double(fpl.out@theta),
#				as.double(fpl.out@sigma),
#				as.integer(1),
#				as.integer(m),
#				as.integer(fpl.out@parm))
#			loq.df <- min(bp.out$xl)
#		}
#		else loq.df <- NA
#	}
	loq.out <- data.frame(loq = loq.df, m = m, cv = cv)
	loq.out
}

loq.thpl.pom <- function(thpl.out, m = thpl.out@m, cv = thpl.out@cv, vlen = 700, mit = 
	10000, toler = 0.001)#, dos = thpl.out@dos)
{
######################################################################
## This function computes the Limits of Quantitation                ##
##  for THPL POM model.                                              ##
##                                                                  ##
## input variables:                                                 ##
##   thpl.out: the output from the fitting function                  ##
##   m: number of reps (number of y-values for each x)              ##
##   cv: the minimum coefficient of variation desired               ##
##   mit, toler: used if dos=F                                      ##
##   vlen: used in computing a range of x below                     ##
##   dos: is T or F, used to change the data type                   ##
##                                                                  ##
## comments added by Matthew Mitchell 3/12/03  have ##              ##
## the original comments have ###                                   ##
######################################################################
  
	cov.un <- thpl.out@cov.unscaled
	b <- thpl.out@coefficients
	n <- length(thpl.out@x)
	b50 <-  b[3]
	xpstart <- min(c(0.0005, min(thpl.out@x[thpl.out@x > 0])))
	xp <- c(seq(xpstart, b50, length = round(vlen/2, 0)), seq(b50, 
		max(thpl.out@x), length = round(vlen/2, 0)))
	yp <- as.vector(thpl.model(xp, b))
		dh.dy <- xp * ((b[2] - b[1])/((b[1] - yp) * (yp - b[2])))
		dh.db1 <- xp/(b[1] - yp)
		dh.db2 <- xp/(yp - b[2])
		dh.db3 <- xp/b[3]
	
	sigma2 <- thpl.out@sigma * thpl.out@sigma
	var.xnot.hat <- (((dh.dy * dh.dy) * sigma2 * ((yp^2)^thpl.out@
		theta))/m + sigma2 * (dh.db1 * (dh.db1 * cov.un[1, 1] + 
		dh.db2 * cov.un[2, 1] + dh.db3 * cov.un[3, 1]) +
                dh.db2 * (dh.db1 * cov.un[1, 2] + 
		dh.db2 * cov.un[2, 2] + dh.db3 * cov.un[3, 2])
                + dh.db3 * (dh.db1 * cov.un[1, 3] + 
		dh.db2 * cov.un[2, 3] + dh.db3 * cov.un[3, 3])))
	sd <- sqrt(var.xnot.hat)
	xp <- xp[is.finite(sd)]
	sd <- sd[is.finite(sd)]
	xmin <- xp[sd/xp == min(sd/xp)]
	#if(dos) {
	if(cv > min(sd/xp[xp != 0]))
		## 01-25-07 Added ties="mean" to approx to avoid error message generated when
		## there are ties. DVS
		loq.df <- approx(((sd/xp) - cv)[xp < xmin], 
			xp[xp < xmin], 0, ties="mean")$y
	else loq.df <- NA
	#}
## 02-08-07: This call to C does not seem relavent. DVS
#	else {
#		xl <- min(xp)
#		xu <- xmin
#		svec <- c(cov.un[, 1], cov.un[, 2], cov.un[, 3])
#		if(!is.loaded(symbol.C("thplpomprec")))
#			dyn.load("/local/lib/calib/sun4/src/thplpom.o")
#		if(cv > min(sd/xp[xp != 0])) {
#			bp.out <- .C("thplpomprec",
#				xl = as.double(xl),
#				xu = as.double(xu),
#				as.double(cv),
#				as.double(b),
#				as.double(svec),
#				as.integer(rep(mit, length(xl))),
#				as.double(toler),
#				as.double(thpl.out@theta),
#				as.double(thpl.out@sigma),
#				as.integer(1),
#				as.integer(m))
#			        loq.df <- min(bp.out$xl)
#		}
#		else loq.df <- NA
#              }
	loq.out <- data.frame(loq = loq.df, m = m, cv = cv)
	loq.out
}

## The derivatives for the quadratic are incorrect - corrections made below by MM 4/15/03
## The mistakes mainly resulted from switching a and c so b[1] and b[3] are switched
## The model is a + bx +cx^2 and a=b[1], b=b[2], c=b[3] and dx/dy reciprocated one part
## incorrectly.  The previous versions have been commented out
## The expressions have simpler forms also when using such facts as
## y=a + bx + cx^2 so y-a = bx + cx^2 and sqrt(b2 - 4ac + 4cy) = b + 2cx

loq.lin.pom <- function(lin.out, m = lin.out@m, cv = lin.out@cv, vlen = 500)
{
	theta <- lin.out@theta
	b <- lin.out@coefficients
	cov.un <- lin.out@cov.unscaled
	sigma <- lin.out@sigma
	xlow <- min(c(0.0005, min(lin.out@x[lin.out@x > 0])))
	if(missing(m))
		m <- lin.out@m
	if(missing(cv)) cv <- lin.out@cv	
	### set up a grid of x values to determine the upper boundary for the search 
###   algorithm, this is the value where CV is at a minimum (called xmin)
	xp <- exp(seq(log(xlow), log(max(lin.out@x)), length = vlen))
	if(lin.out@mmod == "lin") {
		yp <- b[1] + b[2] * xp
		var.xnot.hat <- (sigma * sigma) * ((yp^2)^theta)/((b[2] * b[2]) * m) + 
			(cov.un[1, 1] + 2 * cov.un[1,2] * xp + cov.un[2, 2] * xp * 
			xp)/(b[2]^2)
	}
	else {
		yp <- b[1] + b[2] * xp + b[3] * xp * xp
##		dh.dy <- 1/(xp + (b[2]/(2 * b[1])))
		dh.dy <- 1/(b[2] + 2*b[3]*xp)
##		dh.da <- b[2]/(2 * b[1] * b[1]) + ((yp - b[3])/b[1]) * (
##			1/(xp + (b[2]/(2 * b[1])))) - (xp + (b[2]/(2 * 
##			b[1])))/(2 * b[1] * b[1])
        dh.da <- -1/(b[2] + 2*b[3]*xp)
##		dh.db <- -1/(2 * b[1]) + (b[2]/(2 * b[1])) * (1/(xp + (
##			b[2]/(2 * b[1]))))
        dh.db <- -xp/(b[2] + 2*b[3]*xp)
##		dh.dc <- -1/(xp + (b[2]/(2 * b[1])))
        dh.dc <- -(xp*xp)/(b[2] + 2*b[3]*xp)
		var.xnot.hat <- ((dh.dy * dh.dy) * sigma * sigma * ((yp^2)^
			lin.out@theta))/m + sigma * sigma * (dh.da *
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
		loq.df <- approx(((sd/xp) - cv)[xp < xmin], xp[xp < xmin], 
			0, ties="mean")$y
		loqout <- list(loq = loq.df, m = m, cv = cv)
	}
	else loqout <- list(loq = NA, m = m, cv = cv)
	loqout
}

lof.test <- function(calib.fit.out)
{
########################################################
## This computes the lack of fit statistics from the  ##
##  nls POM fit.                                      ##
## This function gets called as part of the fpl.pom   ##
##  function.                                         ##
##                                                    ##
## comments added by Matthew Mitchell on 3/12/03      ##
########################################################
  
	x <- calib.fit.out@x[is.finite(calib.fit.out@y)]
	y <- calib.fit.out@y[is.finite(calib.fit.out@y)]
	b <- calib.fit.out@coefficients
	nx <- tapply(x, list(x), length)
	vars <- tapply(y, list(x), var)[nx >= 2]
	wts <- tapply(((calib.fit.out@fitted.values^2)^calib.fit.out@theta), 
					list(calib.fit.out@x), mean)[nx >= 2]
	nx <- nx[nx >= 2]
	pure <- (((nx - 1) * vars)/wts)
	pure <- sum(pure)
	df.pure <- sum(nx - 1)
	sse <- sum((calib.fit.out@residuals * calib.fit.out@residuals)/
                   ((calib.fit.out@fitted.values^2)^calib.fit.out@theta))
	df.den <- calib.fit.out@df.residual
	lofss <- sse - pure
	df.lof <- df.den - df.pure
	lof.stat <- (lofss/df.lof)/(pure/df.pure)
	pv <- 1 - pf(lof.stat, df.lof, df.pure)
	lof.out <- new("lof.test", Fstat = lof.stat, p.value = 
		pv, lofss = lofss, df.lof = df.lof, pure.error = pure, 
		df.pure.error = df.pure, sse = sse, df.sse = df.den)
}

