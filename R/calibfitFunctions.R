## calibfit-functions.R
## 
## All class definitions for the package calib.
##
## Author: samarov
###############################################################################
calib.fit <- function(x,y,m=2,cv=0.2,conf=0.95,b3start,b4start,mx=50,ap.calc=T,lof.calc=T,
	#dos=T, 04-17-07: Dropping this from the arguments. Does not seem to serve any purpose. DVS
	theta, method="squared",kmax=10,theta.init=1, digits=3, 
	## 02-08-07: Removed print=T from options. DVS
	#mmod="lin", 05-03-07 Does not appear as though we need this anymore. DVS
	type=c("log.fpl.pom","fpl.pom","log.fpl","fpl","lin.quad.pom","log.thpl.pom","thpl.pom","log.thpl","thpl","lin.pom"),
	intercept=T){
  
	## 01-30-07 This function combines all the various fitting methods for calibration, listed below.
	## Variable definitions can be found in the models below.

	n.type <- length(type)

	## 01-25-07 Making a few changes to function here. Will have an trycatch structure set up where
	## the default in the model fitting process will be
	## (1) fpl.pom 	alt=T
	## (2) fpl.pom 	alt=F
	## (3) fpl	 	alt=T
	## (4) fpl	 	alt=F
	## (5) lin.quad.pom  Added 04-10-07
	## (6) thpl.pom	
	## (7) thpl		
	## (8) lin.pom	
	##
	## 01-30-07 The function will default to the above but the user can specify precisely what they want
	## and in what order. Note that not all options need to be used. DVS
 
  	for(i in 1:n.type){
  	## 04-17-07: Note removed "dos" argument from all functions 
  	
  			if(type[i]=="log.thpl.pom"){
  			
  				calib.fit.out <- try(thpl.pom(x=x, y=y, theta=theta, method=method, m=m, cv=cv, conf=conf,
		        	kmax=kmax, theta.init=theta.init, b3start=b3start, mx=mx, digits=digits,
			        ap.calc=ap.calc, lof.calc=lof.calc,alt=T),TRUE)
  					
  				if(class(calib.fit.out)!="try-error"){
  				
  					if(i!=1)
  						cat(paste("Warning",type[1],"produced error, used",type[i],"instead"),"\n")
  					break
  				}
  			}
			
	
  			if(type[i]=="log.thpl"){
  			
  				calib.fit.out <- try(thpl(x=x, y=y, m=m, cv=cv, conf=conf,
	        		b3start=b3start, mx=mx, ap.calc=ap.calc, lof.calc=lof.calc,alt=T),TRUE)
  					
  				if(class(calib.fit.out)!="try-error"){
  				
  					if(i!=1)
  						cat(paste("Warning",type[1],"produced error, used",type[i],"instead"),"\n")
  					break
  				}
  			}
			
	
  			if(type[i]=="lin.quad.pom"){
  			
  				calib.fit.out <- try(lin.pom(x=x, y=y, theta=theta, method=method, m=m, cv=cv, conf=conf,
  					theta.init=theta.init, mmod="quad", intercept=intercept, kmax=kmax, mx=mx,
  					ap.calc=ap.calc, lof.calc=lof.calc),TRUE)
  					
  				if(class(calib.fit.out)!="try-error"){
  				
  					if(i!=1)
  						cat(paste("Warning",type[1],"produced error, used",type[i],"instead"),"\n")
  					break
  				}
  			}

    	    if(type[i]=="lin.pom"){

    	    	calib.fit.out <- try(lin.pom(x=x, y=y, theta=theta, method=method, m=m, cv=cv, conf=conf,
	        		theta.init=theta.init, mmod="lin", intercept=intercept, kmax=kmax, mx=mx,
			        ap.calc=ap.calc, lof.calc=lof.calc),TRUE)
			    
			    if(class(calib.fit.out)!="try-error"){
        		
        		if(i!=1)
	        		cat(paste("Warning",type[1],"produced error, used",type[i],"instead"),"\n")
		       	break  
		       	}      	
        	}

      		if(type[i]=="thpl.pom"){
      		
      			calib.fit.out <- try(thpl.pom(x=x, y=y, theta=theta, method=method, m=m, cv=cv, conf=conf,
		        	kmax=kmax, theta.init=theta.init, b3start=b3start, mx=mx, digits=digits,
			        ap.calc=ap.calc, lof.calc=lof.calc),TRUE)
			        
			    if(class(calib.fit.out)!="try-error"){
			        
        		if(i!=1)
	        		cat(paste("Warning",type[1],"produced error, used",type[i],"instead"),"\n")
		       	break        	
		       	}
        	}

    	   if(type[i]=="fpl.pom"){
    	   
	    	   calib.fit.out <- try(fpl.pom(x=x, y=y, theta=theta, alt=F, method=method, m=m, cv=cv,
	        		conf=conf, kmax=kmax, theta.init=theta.init, b3start=b3start, b4start=b4start,
			        mx=mx, digits=digits, ap.calc=ap.calc, lof.calc=lof.calc),TRUE)
        			
        		if(class(calib.fit.out)!="try-error"){
        			
        		if(i!=1)
	        		cat(paste("Warning",type[1],"produced error, used",type[i],"instead"),"\n")
		       	break        
		       	}	
        	}

    	   if(type[i]=="log.fpl.pom"){
    	   
	    	   calib.fit.out <- try(fpl.pom(x=x, y=y, theta=theta, alt=T, method=method, m=m, cv=cv,
    	    		conf=conf, kmax=kmax, theta.init=theta.init, b3start=b3start, b4start=b4start,
			        mx=mx, digits=digits, ap.calc=ap.calc, lof.calc=lof.calc),TRUE)
	        	
	        	if(class(calib.fit.out)!="try-error"){
	        	
        		if(i!=1)
	        		cat(paste("Warning",type[1],"produced error, used",type[i],"instead"),"\n")
		       	break        	
		       	}
        	}

    	   if(type[i]== "thpl"){
    	   
    		   calib.fit.out <- try(thpl(x=x, y=y, m=m, cv=cv, conf=conf,
	        		b3start=b3start, mx=mx, ap.calc=ap.calc, lof.calc=lof.calc),TRUE)
        		
        		if(class(calib.fit.out)!="try-error"){

        		if(i!=1)
	        		cat(paste("Warning",type[1],"produced error, used",type[i],"instead"),"\n")
		       	break        	
		       	}
        	}

    	   if(type[i]=="log.fpl"){
    	   	  calib.fit.out <- try(fpl(x=x, y=y, alt=T, m=m, cv=cv,
		       			conf=conf, b3start=b3start, b4start=b4start,
		    		    mx=mx, ap.calc=ap.calc, lof.calc=lof.calc),TRUE)
        				
        	  if(class(calib.fit.out)!="try-error"){
       		 	
        		if(i!=1)
	        		cat(paste("Warning",type[1],"produced error, used",type[i],"instead"),"\n")
	        	break        	
    	    	}
        	}

    	   if(type[i]=="fpl"){
    	   
	    	   calib.fit.out <- try(fpl(x=x, y=y, alt=F, m=m, cv=cv,
			        conf=conf, b3start=b3start, b4start=b4start,
    			    mx=mx, ap.calc=ap.calc, lof.calc=lof.calc),TRUE)
        		
        		if(class(calib.fit.out)!="try-error"){
        		
        		if(i!=1)
	        		cat(paste("Warning",type[1:(i-1)],"produced error, used",type[i],"instead"),"\n")
		       	break        	
		       	}
        	}

			if((i==n.type)&(class(calib.fit.out)=="try-error"))
	    			stop("calib is unable to fit a model to the data")
    
  	}
	calib.fit.out
}
  


fpl <- function(x, y, alt = F, m = 2, cv = 0.2, conf = 0.95, 
	b3start, b4start, mx = 50, ap.calc = T, lof.calc = T, lowLim = 10e-6)#, dos = T) 04-17-07
	## Removing dos argument
{
####################################################################
## The 4 Parameter Logistic Function used in Calibration Analysis ##
## modified by Matthew Mitchell on 3/11/03                        ##
##                                                                ##
## input variables:                                               ##
##  x: the independent variable (usually dosage here)             ##
##  y: the response variable (example: OD, optical density)       ##
##  alt:  If alt=F the nonlinear function fitted is               ##
##        y = b2 + (b1-b2)/(1 + (x/b3)^b4)                        ##
##        If alt=T the nonlinear function fitted is               ##
##        y = b2 + (b1-b2)/(1 + exp(b4(logx - b3)),               ##       
##        which is reparametrization of the other form            ##
##        since b3 from alt=T is log(b3 from alt=F).              ##
##  m: This is the number of replicates (number of y's at         ##
##     each x).                                                   ##
##  cv: This is the acceptable coefficient of variation.          ##
##      The limits of quantitation are calculated with            ##
##      this constraint.                                          ##
##  conf: The confidence level used for any prediction interval   ##
##  b3start: the initial value used for b3 which corresponds      ##
##          to the model when alt=F.  This gets converted         ##
##          to log(b3start) if alt=T.                             ##
##  b4start: the initial value for b4                             ##
##  mx: the maximum number of iterations used in the non-linear   ##
##     least-squares fit.                                         ##
##  ap.calc: if T, the assay performance measures will be         ##
##          computed (mdc, rdl, loq)                              ##
##  lof.calc: if T, the lack of fit statistics are computed       ##
##  dos: used to make sure the data types are okay                ##
##                                                                ##
## The object class of the output is fpl.pom,pom                  ##
## output variables:                                              ##
##  coefficients: the estimated values for b1,b2,b3,b4            ##
##  se.coefficients: the standard errors of these                 ##
##  sigma: the residual standard error                            ##
##  cov.unscaled: this is the variance-covariance matrix          ##
##   of the estimates of b1,b2,b3,b4 without the sigma^2          ##
##  df.residual: the residual degrees of freedom                  ##
##  residuals: all the residuals from the fit                     ##
##  fitted.values: the fitted values from the fit                 ##
##  theta: This is the power of the mean (POM) parameter          ##
##  and is equal to 0 since OLS is used here.                     ##
##  vfe.method: variance function estimation method.  Here it     ##
##   is OLS (ordinary least squares)                              ##
##  status: indicates whether the non-linear least squares        ##
##          converged                                             ##
##  kused: number of iterations used to estimate theta.           ##
##         Here it is 0 because OLS was used.                     ##
##  x,y,m,cv: the same values that were input above               ##
##  parm: if alt=F this is 1, if alt=T this is 2                  ##
##  cf: some measure of the assay performace?                     ##
##  mdc: the minimum detectable concentration                     ##
##  rdl: the reliable detection limits                            ##
##  loq: the limits of quantitation                               ##
##                                                                ##
## Functions called within this function:                         ##
##  fpl.model, fpl.inverse, mdc, rdl, loq                         ##
## Original comments are ###, added are ##                        ##
####################################################################
    require(stats)

	## added by MM on 4/15/03 to prevent taking log of 0
    x[x == 0] <- lowLim
   
    fpl.frame <- data.frame(x = x, y = y)
    
	##  set up initial values for parameters
    xmax <- max(x)
    b1start <- mean(y[x == min(x)])
    b2start <- mean(y[x == xmax])
    if(missing(b4start)) b4start <- 1
    if(missing(b3start)) {
        ymid <- (b1start + b2start)/2
        y.diff <- abs(y - ymid)
        min.diff <- min(y.diff)
        ymid <- (y[y.diff == min.diff])[1]
        xmid <- (x[y.diff == min.diff])[1]
        b3start <- xmid/(((b1start - b2start)/(ymid - b2start) - 1)^(1/b4start))
        if(alt) b3start <- log(b3start)
    }
    attr(fpl.frame, "param") <- list(b1=b1start, b2=b2start, b3=b3start, b4=b4start)

	## The non-linear fit is performed
    if(!alt) {
        nlsfpl <- nls(y ~ fpl.model(x, b1, b2, b3, b4, parm = 1),
            fpl.frame, control = nls.control(maxiter = mx),
            start=attr(fpl.frame, "param"))
        parm <- 1
    }
    else {
        nlsfpl <- nls(y ~ fpl.model(x, b1, b2, b3, b4, parm = 2),
            fpl.frame, control = nls.control(maxiter = mx),
            start=attr(fpl.frame, "param"))
        parm <- 2
    }
	## summary is  formula, parameters(est, se, t, pvalue), residual se,
	##  and correlation of parameter estimates    
    summ.fpl <- summary(nlsfpl)
    b <- coef(nlsfpl)
    fitted <- fpl.model(x, b, parm = parm)
    grad <- attributes(fitted)$gradient ## this is the gradient matrix
    fitted <- as.vector(fitted)
    if(missing(m))
        m <- length(x[x == min(x)])
        
    ## 01-30-07 Making a few changes to the slots of objects of class of calib.fit.
    ## Does not seem necessary to have the gradient stored as this is not really used in any other
    ## part of the package.
    fpl.out <- new("calib.fit",coefficients = b, se.coefficients = 
    	as.vector(summ.fpl$parameters[, 2]), sigma = as.vector(summ.fpl$
        sigma), cov.unscaled = as.matrix(summ.fpl$cov.unscaled),
        df.residual = as.vector(summ.fpl$df)[2], residuals = y - 
        fitted, fitted.values = fitted, theta = 0, vfe.method
         = "OLS", status = "converged", kused = 0, x = x, y = y,
        m = m, cv = cv, parm = parm,type="fpl")
    #fpl.out@dos = dos 04-17-07: Removing dos arg DVS
    fpl.out@gradient <- grad
    ### if ap.calc=T, compute assay performance statistics
    #browser()
    if(ap.calc) {
	### Here's a good one: the automation interface has a problem
	### with the "single" datatype. the next three values come
	### back through as singles but I'm coercing them to doubles
	### otherwise we get really messed up values in Excel.
	### (TDS, 990104)
        fpl.out@mdc <- as.double(mdc.calc(fpl.out, conf = conf,type="fpl")$mdc)
        fpl.out@rdl <- as.double(rdl.calc(fpl.out, conf = conf,type="fpl")$rdl)

		if(fpl.out@rdl > max(x))
			fpl.out@rdlwarn <- "RDL outside range of values"
		else
			fpl.out@rdlwarn <- "RDL within range of values"

        fpl.out@loq <- as.double(loq.calc(fpl.out,type="fpl")$loq)
        fpl.out@conf.level = conf
        y.cf <- y[x != 0]
        x.cf <- x[x != 0]
        cal.x <- fpl.inverse(fpl.out, y.cf) 
        diff.x <- abs(cal.x - x.cf)
        fpl.out@cf <- 1 - (1/length(x.cf)) * sum(diff.x/x.cf)
    }
    ## 04-17-07: Does not seem necessary to assign class NULL to mdc, rdl and loq. Looking to
    ## take this out. DVS
    else {
        fpl.out@mdc <- NULL
        fpl.out@rdl <- NULL
        fpl.out@loq <- NULL
    }
     if(lof.calc) {
        if(length(x) > length(unique(x)))
            fpl.out@lof.test <- lof.test(fpl.out)
        else fpl.out@lof.test <- NULL
    }
    #else attr(fpl.out, "lof.test") <- NULL
    #xlout <<- fpl.out   
    ### added by TDS, 990104 to make this data available to Excel
    fpl.out
}

thpl <- function(x, y, m = 2, cv = 0.2, conf = 0.95, b3start, 
    mx = 25, ap.calc = T, lof.calc = T, alt = F, lowLim = 10e-6)#, dos = T)
{
####################################################################
## The 3 Parameter Logistic Function used in Calibration Analysis ##
## written by Matthew Mitchell on 5/20/03                         ##
## developed for the fluorescence curves fitted by Biosense       ##
##                                                                ##
## b1 = F_0, b2 = F_inf, b3 = K_d                                 ##
##                                                                ##
## input variables:                                               ##
##  x: the independent variable (usually glucose conc)            ##
##  y: the response variable (example: fluorescence intensity)    ##
##  m: This is the number of replicates (number of y's at         ##
##     each x).                                                   ##
##  cv: This is the acceptable coefficient of variation.          ##
##      The limits of quantitation are calculated with            ##
##      this constraint.                                          ##
##  conf: The confidence level used for any prediction interval   ##
##  b3start: the initial value used for b3                        ##
##  mx: the maximum number of iterations used in the non-linear   ##
##     least-squares fit.                                         ##
##  ap.calc: if T, the assay performance measures will be         ##
##          computed (mdc, rdl, loq)                              ##
##  lof.calc: if T, the lack of fit statistics are computed       ##
##  dos: used to make sure the data types are okay                ##
##                                                                ##
## The object class of the output is thpl.pom,pom                 ##
## output variables:                                              ##
##  coefficients: the estimated values for b1,b2,b3               ##
##  se.coefficients: the standard errors of these                 ##
##  sigma: the residual standard error                            ##
##  cov.unscaled: this is the variance-covariance matrix          ##
##   of the estimates of b1,b2,b3    without the sigma^2          ##
##  df.residual: the residual degrees of freedom                  ##
##  residuals: all the residuals from the fit                     ##
##  fitted.values: the fitted values from the fit                 ##
##  theta: This is the power of the mean (POM) parameter          ##
##  and is equal to 0 since OLS is used here.                     ##
##  vfe.method: variance function estimation method.  Here it     ##
##   is OLS (ordinary least squares)                              ##
##  status: indicates whether the non-linear least squares        ##
##          converged                                             ##
##  kused: number of iterations used to estimate theta.           ##
##         Here it is 0 because OLS was used.                     ##
##  x,y,m,cv: the same values that were input above               ##
##  cf: some measure of the assay performace?                     ##
##  mdc: the minimum detectable concentration                     ##
##  rdl: the reliable detection limits                            ##
##  loq: the limits of quantitation                               ##
##                                                                ##
## Functions called within this function:                         ##
##  thpl.model, thpl.inverse, mdc, rdl, loq                       ##
## Original comments are ###, added are ##                        ##
####################################################################
    require(stats)

	x[x == 0] <- lowLim

    thpl.frame <- data.frame(x = x, y = y)
    
	##  set up initial values for parameters
    xmax <- max(x)
    b1start <- mean(y[x == min(x)])
    b2start <- mean(y[x == xmax])
    if(missing(b3start)) {
        ymid <- (b1start + b2start)/2
        y.diff <- abs(y - ymid)
        min.diff <- min(y.diff)
        ymid <- (y[y.diff == min.diff])[1]
        xmid <- (x[y.diff == min.diff])[1]
        b3start <- xmid/(((b1start - b2start)/(ymid - b2start) - 1))
		if(alt) b3start <- log(b3start)
    }
    attr(thpl.frame, "param") <- list(b1=b1start, b2=b2start, b3=b3start)

	## The non-linear fit is performed
	if(!alt) {
    	nlsthpl <- nls(y ~ thpl.model(x, b1, b2, b3, parm = 1),
    		thpl.frame, control = nls.control(maxiter = mx),
        	start=attr(thpl.frame, "param"))
		parm <- 1
	}
	else {
		nlsthpl <- nls(y ~ thpl.model(x, b1, b2, b3, parm = 2),
			thpl.frame, control = nls.control(maxiter = mx),
			start = attr(thpl.frame, "param"))
		parm <- 2
	}
     
	## summary is  formula, parameters(est, se, t, pvalue), residual se,
	##  and correlation of parameter estimates    
    summ.thpl <- summary(nlsthpl)
    b <- coef(nlsthpl)
    fitted <- thpl.model(x, b)
    grad <- attributes(fitted)$gradient ## this is the gradient matrix
    fitted <- as.vector(fitted)
    if(missing(m))
        m <- length(x[x == min(x)])

    thpl.out <- new("calib.fit",coefficients = b, se.coefficients = as.vector(
        summ.thpl$parameters[, 2]), sigma = as.vector(summ.thpl$
        sigma), cov.unscaled = as.matrix(summ.thpl$cov.unscaled),
        df.residual = as.vector(summ.thpl$df)[2], residuals = y - 
        fitted, fitted.values = fitted, theta = 0, vfe.method
         = "OLS", status = "converged", kused = 0, x = x, y = y,
        m = m, cv = cv, type="thpl")
	if(!alt)
        thpl.out@parm <- 1
    else thpl.out@parm <- 2

    #thpl.out@dos <- dos
    thpl.out@gradient <- grad   
    ### if ap.calc=T, compute assay performance statistics
    if(ap.calc) {
	### Here's a good one: the automation interface has a problem
	### with the "single" datatype. the next three values come
	### back through as singles but I'm coercing them to doubles
	### otherwise we get really messed up values in Excel.
	### (TDS, 990104)
        thpl.out@mdc <- as.double(mdc.calc(thpl.out, conf = conf,type="thpl")$mdc)
        thpl.out@rdl <- as.double(rdl.calc(thpl.out, conf = conf,type="thpl")$rdl)
		
		if(thpl.out@rdl > max(x))
			thpl.out@rdlwarn <- "RDL outside range of values"
		else
			thpl.out@rdlwarn <- "RDL within range of values"
		
        thpl.out@loq <- as.double(loq.calc(thpl.out, type="thpl")$loq)
        thpl.out@conf.level <- conf
        y.cf <- y[x != 0]
        x.cf <- x[x != 0]
        cal.x <- thpl.inverse(thpl.out, y.cf) 
        diff.x <- abs(cal.x - x.cf)
        thpl.out@cf <- 1 - (1/length(x.cf)) * sum(diff.x/x.cf)
    }
    else {
        thpl.out@mdc <- NULL
        thpl.out@rdl <- NULL
        thpl.out@loq <- NULL
    }
    if(lof.calc) {
      if(length(x) > length(unique(x)))
          thpl.out@lof.test <- lof.test(thpl.out)
	}
    ### added by TDS, 990104 to make this data available to Excel
    thpl.out
}


fpl.pom <- function(x, y, theta, alt = F, method = "squared", m = 2, cv = 0.2, 
    conf = 0.95, kmax = 10, theta.init = 1, b3start, b4start, mx = 50,
    digits = 3, ap.calc = T, lof.calc = T, lowLim = 10e-6)#, dos = T)
    ## 02-08-07: Removed print=F from options.
{
####################################################################
## The 4 Parameter Logistic Function used in Calibration Analysis ##
## modified by Matthew Mitchell on 3/11/03                        ##
##                                                                ##
## input variables:                                               ##
##  x: the independent variable (usually dosage here)             ##
##  y: the response variable (example: OD, optical density)       ##
##  theta: if the value for the power of the mean (POM)           ##
##   is known give the value here.  Otherwise, leave              ##
##   it missing.                                                  ##
##  alt:  If alt=F the nonlinear function fitted is               ##
##        y = b2 + (b1-b2)/(1 + (x/b3)^b4)                        ##
##        If alt=T the nonlinear function fitted is               ##
##        y = b2 + (b1-b2)/(1 + exp(b4(logx - b3)),               ##       
##        which is reparametrization of the other form            ##
##        since b3 from alt=T is log(b3 from alt=F).              ##
##  method: this is the method used in the optimization           ##
##          process when estimating theta.  The values            ##
##          for method can be "absolute", "squared",              ##
##          "ols" or "reml."  The default is "squared."           ##
##          The first two methods are pseudo-likelihood           ##
##          based.                                                ##
##  m: This is the number of replicates (number of y's at         ##
##     each x).                                                   ##
##  cv: This is the acceptable coefficient of variation.          ##
##      The limits of quantitation are calculated with            ##
##      this constraint.                                          ##
##  conf: The confidence level used for any prediction interval   ##
##  kmax: the maximum number of iterations used to estimate       ##
##        theta.  The default was originally 5 but I have         ##
##        found that is usually too small so I changed it to 10.  ##
##  theta.init: the initial value set for theta, the              ##
##              POM parameter.  The default is 1.                 ##
##  b3start: the initial value used for b3 which corresponds      ##
##          to the model when alt=F.  This gets converted         ##
##          to log(b3start) if alt=T.                             ##
##  b4start: the initial value for b4                             ##
##  mx: the maximum number of iterations used in the non-linear   ##
##     least-squares fit.                                         ##
##  print: if T prints the results from the iterations            ##
##  ap.calc: if T, the assay performance measures will be         ##
##          computed (mdc, rdl, loq)                              ##
##  lof.calc: if T, the lack of fit statistics are computed       ##
##  dos: used to make sure the data types are okay                ##
##                                                                ##
## The object class of the output is fpl.pom,pom                  ##
## output variables:                                              ##
##  coefficients: the estimated values for b1,b2,b3,b4            ##
##  se.coefficients: the standard errors of these                 ##
##  sigma: the residual standard error                            ##
##  cov.unscaled: this is the variance-covariance matrix          ##
##   of the estimates of b1,b2,b3,b4 without the sigma^2          ##
##  theta: estimate of theta the POM parameter                    ##
##  df.residual: the residual degrees of freedom                  ##
##  residuals: all the residuals from the fit                     ##
##  fitted.values: the fitted values from the fit                 ##
##  theta: This is the power of the mean (POM) parameter          ##
##  and is equal to 0 since OLS is used here.                     ##
##  vfe.method: variance function estimation method.              ##
##  status: indicates whether the non-linear least squares        ##
##          converged                                             ##
##  kused: number of iterations used to estimate theta.           ##
##  x,y,m,cv: the same values that were input above               ##
##  parm: if alt=F this is 1, if alt=T this is 2                  ##
##  mdc: the minimum detectable concentration                     ##
##  rdl: the reliable detection limits                            ##
##  loq: the limits of quantitation                               ##
##  The lof (lack of fit) information is an attribute             ##
##   of type "lof.test."                                          ##
##                                                                ##
## Functions called within this function:                         ##
##  fpl, fplpom, vfemin.abs.pom, vfemin.pom, reml.pom,            ##
##  glsvfe.print, fpl.model, fpl.inverse, mdc, rdl, loq, lof      ##
## Original comments are ###, added are ##   					  ##
##																  ##
## Comments made on 02-08-07 see code below. DVS                  ##
####################################################################

	## MM added this line on 4/15/03 to prevent taking log of 0  
    x[x == 0] <- lowLim

	## 02-08-07 The following bit of code no longer seems relavent as it
	## making a call to some C code which is no longer in use. DVS
    #if(!dos) {
    #    if(!is.loaded(symbol.C("fplpomcal"))) {
    #        cat("Loading fplpom.o\n")
    #        dyn.load("/local/lib/calib/sun4/src/fplpom.o")
    #    }
    #}
    if(!all(is.finite(y)))
        stop("Missing values in y not allowed")
    if(!all(is.finite(x)))
        stop("Missing values in x not allowed")
    if(missing(m))
        m <- length(x[x == min(x)])
    if(!missing(theta))
        tfixed <- T
    else tfixed <- F
        k <- 0
    if(missing(b4start))
        b4start <- 1
    if(missing(b3start)) {
        ymid <- (mean(y[x == min(x)]) + mean(y[x == max(x)]))/2
        y.diff <- abs(y - ymid)
        min.diff <- min(y.diff)
        ymid <- (y[y.diff == min.diff])[1]
        xmid <- (x[y.diff == min.diff])[1]
        b3start <- xmid/(((mean(y[x == min(x)]) - mean(y[x == 
            max(x)]))/(ymid - mean(y[x == max(x)])) - 1)^(1/
            b4start))
    }
    if(alt)
        b3start <- log(b3start)
    
    fpl.out <- fpl(x, y, b3start = b3start, b4start = b4start, alt
         = alt, lof.calc = F, ap.calc = F)
    coef <- fpl.out@coefficients
    sigma <- fpl.out@sigma
    se.bhat <- fpl.out@se.coefficients
    ## 02-08-07: Removed print from options. DVS
    #if(print == T) {
    #    glsvfe.print(k, coef, se.bhat, sigma, digits = digits)
    #}
    wts <- 1
    thold <- 0
    sigold <- sigma
    m.target <- c("absolute", "squared", "reml", "ols", "")
    m.ind <- charmatch(method, m.target)
    method <- m.target[m.ind]
    if(is.na(m.ind))
        stop("\nmethod must be one of:  \"absolute\", \"squared\", \"reml\", or \"ols\""
            )
    ## Theta and b1,b2,b3,b4 are estimated through iteratively re-weighted least-squares
   
    if(method != "ols") {
        for(k in 1:kmax) {

            pred <- fpl.out@fitted.values

            res <- fpl.out@residuals

            if(!tfixed) {
                if(method == "absolute")
                    nl.out <- optim(theta.init, vfemin.abs.pom, method="BFGS",pred=pred,res=res)

                if(method == "squared")
                    nl.out <- optim(theta.init, vfemin.pom, method="BFGS", pred=pred, res=res)

                if(method == "reml") {
                  assign("grad", fpl.out@gradient)
                  #if(dos==F)
                  #  nl.out <- optim(theta.init, reml.pom, method="BFGS",pred=pred)
                  #else 
                  nl.out <- optim(theta.init, reml.pom, method="BFGS",pred=pred)#,max.facl=60, maxiter=30)
                  }
                  
                   ## 02-08-07 The following was not working in estimating theta using REML. Replaced code
                   ## below with singe "#" with code above. DVS
                   #if(dos == F)
                   #  nl.out <- optim(theta.init, reml.pom, method="BFGS")
                   #else nl.out <- nlmin(remld.pom, theta.init, max.fcal = 60, max.iter = 30)
	               #

                if(nl.out$convergence==0)
                  theta <- nl.out$par
                else {
                  cat(paste("theta failed to converge:", nl.out$message, "\n\n"))
                  return()
                }
                thnew <- theta
                if(abs(thnew - thold) < 0.001) {
                  status <- "converged"
                  #cat("theta converged\n\n") 
                  # 05-14-07: No need for this display. DVS
                  break
                }
                thold <- thnew
            }
            its <- k
            ## 04-03-07: So making a few changes here. It seems as though that throughout
            ## the calib program we have situations like that below. However it seems to
            ## me, at least based off of the theory that this is developed off of that the
            ## weights should be strictly greater than zero. However by having 2 * theta 
            ## like below allows for the possibility of non-positive weights. Changing
            ## from 1/pred^(2*theta) to what is now below. DVS
            wts <- 1/((pred^2)^theta)
            fpl.out <- fplpom(x, y, wt = wts, theta = theta, 
                            b3start = b3start, b4start = b4start, alt = alt, mx = mx)
            ## 01-23-07: Made a change to the line abobve allowing the user to specify the
            ## maximum number of iterations in the fitting. DVS
            
            ## 02-09-07: There does not seem to be any reason to have a print option. DVS
            #if(print == T)
            #    glsvfe.print(k, fpl.out@coefficients, fpl.out@se.coefficients,
            #                 fpl.out@sigma, theta, digits = digits)
            if(tfixed) {
                signew <- fpl.out@sigma
                if(abs(sigold - signew) < 0.0001) {
                  status <- "converged"
                  cat("sigma converged\n\n")
                  break
                }
                sigold <- signew
            }
            if(k == kmax)
                cat(paste(
                  "GLS algorithm failed to converge for",
                  kmax, "iterations\n\n"))
            status <- "failed to converge"
            ## 04-05-07: Adding the following line. Hoping that it adds some stability
            ## to the fitting algorithm.
            ## Taking out for the time being
            # theta.init <- thold
        }
    }
    else {
        theta <- 0
        status <- "OLS fit"
    }
    
    meth <- "PL: squared residuals"
    if(tfixed)
        meth <- "theta fixed"
    if(method == "absolute")
        meth <- "PL: absolute residuals"
    if(method == "reml")
        meth <- "REML"
    if(method == "ols")
        meth <- "OLS"

    fplvfe.out <- new("calib.fit",coefficients = fpl.out@coefficients, se.coefficients
         = fpl.out@se.coefficients, sigma = fpl.out@sigma, cov.unscaled
         = fpl.out@cov.unscaled, theta = theta, df.residual = fpl.out@
        df.residual, fitted.values = fpl.out@fitted.values, residuals
         = fpl.out@residuals, vfe.method = meth, kused = its, 
        status = status, x = x, y = y,type="fpl")
    #fplvfe.out@dos <- dos
    fplvfe.out@gradient <- fpl.out@gradient
    if(!alt)
        fplvfe.out@parm <- 1
    else fplvfe.out@parm <- 2
    
    fplvfe.out@m <- m
    fplvfe.out@cv <- cv
	## computation of assay performance measures here
    if(ap.calc) {
        fplvfe.out@mdc <- mdc.calc(fplvfe.out, conf = conf,type="fpl")$mdc
        fplvfe.out@rdl <- rdl.calc(fplvfe.out, conf = conf,type="fpl")$rdl
		
		if(fplvfe.out@rdl > max(x))
			fplvfe.out@rdlwarn <- "RDL outside range of values"
		else
			fplvfe.out@rdlwarn <- "RDL within range of values"
		
        fplvfe.out@loq <- loq.calc(fplvfe.out,type="fpl")$loq
    }
    else {
        fplvfe.out@mdc <- NULL
        fplvfe.out@rdl <- NULL
        fplvfe.out@loq <- NULL
    }
	## Lack of fit computations here

    if(lof.calc) {
        if(length(x) > length(unique(x)))
            fplvfe.out@lof.test <- lof.test(fplvfe.out)
       # Don't think the following line is necessary DVS 10/18/06
       #else fplvfe.out@lof.test <- NULL
    }
    fplvfe.out
}

thpl.pom <- function(x, y, theta, alt = F, method = "squared", m = 2, cv = 0.2, 
    conf = 0.95, kmax = 10, theta.init = 1, b3start, mx = 
    25, digits = 3, ap.calc = T, lof.calc = T, lowLim = 10e-6)#, dos = T)
{
## 02-08-07: Removed print from options. DVS
####################################################################
## The 3 Parameter Logistic Function used in Calibration Analysis ##
## written  by Matthew Mitchell on 5/20/03                        ##
##                                                                ##
## input variables:                                               ##
##  x: the independent variable (usually dosage here)             ##
##  y: the response variable (example: OD, optical density)       ##
##  theta: if the value for the power of the mean (POM)           ##
##   is known give the value here.  Otherwise, leave              ##
##   it missing.                                                  ##
##  method: this is the method used in the optimization           ##
##          process when estimating theta.  The values            ##
##          for method can be "absolute", "squared",              ##
##          "ols" or "reml."  The default is "squared."           ##
##          The first two methods are pseudo-likelihood           ##
##          based.                                                ##
##  m: This is the number of replicates (number of y's at         ##
##     each x).                                                   ##
##  cv: This is the acceptable coefficient of variation.          ##
##      The limits of quantitation are calculated with            ##
##      this constraint.                                          ##
##  conf: The confidence level used for any prediction interval   ##
##  kmax: the maximum number of iterations used to estimate       ##
##        theta.  The default was originally 5 but I have         ##
##        found that is usually too small so I changed it to 10.  ##
##  theta.init: the initial value set for theta, the              ##
##              POM parameter.  The default is 1.                 ##
##  b3start: the initial value used for b3                        ##
##  mx: the maximum number of iterations used in the non-linear   ##
##     least-squares fit.                                         ##
##  print: if T prints the results from the iterations            ##
##  ap.calc: if T, the assay performance measures will be         ##
##          computed (mdc, rdl, loq)                              ##
##  lof.calc: if T, the lack of fit statistics are computed       ##
##  dos: used to make sure the data types are okay                ##
##                                                                ##
## The object class of the output is thpl.pom,pom                 ##
## output variables:                                              ##
##  coefficients: the estimated values for b1,b2,b3               ##
##  se.coefficients: the standard errors of these                 ##
##  sigma: the residual standard error                            ##
##  cov.unscaled: this is the variance-covariance matrix          ##
##   of the estimates of b1,b2,b3    without the sigma^2          ##
##  theta: estimate of theta the POM parameter                    ##
##  df.residual: the residual degrees of freedom                  ##
##  residuals: all the residuals from the fit                     ##
##  fitted.values: the fitted values from the fit                 ##
##  theta: This is the power of the mean (POM) parameter          ##
##  and is equal to 0 since OLS is used here.                     ##
##  vfe.method: variance function estimation method.              ##
##  status: indicates whether the non-linear least squares        ##
##          converged                                             ##
##  kused: number of iterations used to estimate theta.           ##
##  x,y,m,cv: the same values that were input above               ##
##  mdc: the minimum detectable concentration                     ##
##  rdl: the reliable detection limits                            ##
##  loq: the limits of quantitation                               ##
##  The lof (lack of fit) information is an attribute             ##
##   of type "lof.test."                                          ##
##                                                                ##
## Functions called within this function:                         ##
##  thpl, thplpom, vfemin.abs.pom, vfemin.pom, reml.pom,          ##
##  glsvfe.print, thpl.model, thpl.inverse, mdc, rdl, loq, lof    ##
## Original comments are ###, added are ##                        ##
####################################################################
	
	x[x==0] <- lowLim
	## 02-08-07 The following bit of code no longer seems relavent as it
	## making a call to some C code which is no longer in use. DVS
    #if(!dos) {
    #    if(!is.loaded(symbol.C("thplpomcal"))) {
    #        cat("Loading thplpom.o\n")
    #        dyn.load("/local/lib/calib/sun4/src/thplpom.o")
    #    }
    #}
    if(!all(is.finite(y)))
        stop("Missing values in y not allowed")
    if(!all(is.finite(x)))
        stop("Missing values in x not allowed")
    if(missing(m))
        m <- length(x[x == min(x)])
    if(!missing(theta))
        tfixed <- T
    else tfixed <- F
    k <- 0
	
    if(missing(b3start)) {
        ymid <- (mean(y[x == min(x)]) + mean(y[x == max(x)]))/2
        y.diff <- abs(y - ymid)
        min.diff <- min(y.diff)
        ymid <- (y[y.diff == min.diff])[1]
        xmid <- (x[y.diff == min.diff])[1]
        b3start <- xmid/(((mean(y[x == min(x)]) - mean(y[x == 
            max(x)]))/(ymid - mean(y[x == max(x)])) - 1))
		if(alt) b3start <- log(b3start)
    }
    
	thpl.out <- thpl(x, y, b3start = b3start, lof.calc = F, ap.calc = F, alt = alt)
    coef <- thpl.out@coefficients
    sigma <- thpl.out@sigma
    se.bhat <- thpl.out@se.coefficients
    ## 02-08-07: Removed print from options. DVS
    #if(print == T) {
    #    glsvfe.print(k, coef, se.bhat, sigma, digits = digits, mod=3)
    #}
    wts <- 1
    thold <- 0
    sigold <- sigma
    m.target <- c("absolute", "squared", "reml", "ols", "")
    m.ind <- charmatch(method, m.target)
    method <- m.target[m.ind]
    if(is.na(m.ind))
        stop("\nmethod must be one of:  \"absolute\", \"squared\", \"reml\", or \"ols\""
            )
	## Theta and b1,b2,b3 are estimated through iteratively re-weighted least-squares
    
    if(method != "ols") {
        for(k in 1:kmax) {

            pred <- thpl.out@fitted.values

            res <- thpl.out@residuals

            if(!tfixed) {
                if(method == "absolute")

                    nl.out <- optim(theta.init, vfemin.abs.pom, method="BFGS", pred=pred, res=res)
                if(method == "squared")

                    nl.out <- optim(theta.init, vfemin.pom, method="BFGS", pred=pred, res=res)

                if(method == "reml") {
                  assign("grad", attributes(thpl.out)$gradient)
                  #if(dos==F) 04-17-07 Removing dos from options. DVS
                  #  nl.out <- optim(theta.init, reml.pom, method="BFGS",pred=pred)
                  #else 
                  nl.out <- optim(theta.init, reml.pom, method="BFGS",pred=pred)#max.facl=60, maxiter=30)
                  }
                  ## 02-08-07 The following was not working in estimating theta using REML. Replaced code
                  ## below with singe "#" with code above. DVS
                  #if(dos == F)
				  #
                  #  nl.out <- optim(theta.init, reml.pom, method="BFGS")
                  #else nl.out <- nlmin(remld.pom, theta.init, max.fcal = 60, max.iter = 30)
                  #

                if(nl.out$convergence==0)
                  theta <- nl.out$par
                else {
                  cat(paste("theta failed to converge:", nl.out$message, "\n\n"))
                  return()
                }
                thnew <- theta
                if(abs(thnew - thold) < 0.001) {
                  status <- "converged"
                  #cat("theta converged\n\n")
                  break
                }
                thold <- thnew
                # 04-03-07: Removed brack from next line DVS
            } 
            its <- k
            ## 04-03-07: So making a few changes here. It seems as though that throughout
            ## the calib program we have situations like that below. However it seems to
            ## me, at least based off of the theory that this is developed off of that the
            ## weights should be strictly greater than zero. However by having 2 * theta 
            ## like below allows for the possibility of non-positive weights. Changing
            ## from 1/pred^(2*theta) to what is now below. DVS
            wts <- 1/((pred^2)^theta)
            thpl.out <- thplpom(x, y, wt = wts, theta = theta, 
                            b3start = b3start, alt = alt, mx = mx)
            ## 02-08-07: Removed print from options. DVS
            #if(print == T)
            #    glsvfe.print(k, thpl.out@coefficients, thpl.out@
            #      se.coefficients, thpl.out@sigma, theta, digits
            #       = digits, mod=3)
            if(tfixed) {
                signew <- thpl.out@sigma
                if(abs(sigold - signew) < 0.0001) {
                  status <- "converged"
                  cat("sigma converged\n\n")
                  break
                }
                sigold <- signew
            }
            if(k == kmax)
                cat(paste(
                  "GLS algorithm failed to converge for",
                  kmax, "iterations\n\n"))
            status <- "failed to converge"
        }
    }
    else{
        theta <- 0
        status <- "OLS fit"
    }
	meth <- "PL: squared residuals"
    if(tfixed)
        meth <- "theta fixed"
    if(method == "absolute")
        meth <- "PL: absolute residuals"
    if(method == "reml")
        meth <- "REML"
    if(method == "ols")
        meth <- "OLS"
    
    thplvfe.out <- new("calib.fit",coefficients = thpl.out@coefficients, se.coefficients
         = thpl.out@se.coefficients, sigma = thpl.out@sigma, cov.unscaled
         = thpl.out@cov.unscaled, theta = theta, df.residual = thpl.out@
        df.residual, fitted.values = thpl.out@fitted.values, residuals
         = thpl.out@residuals, vfe.method = meth, kused = its, 
        status = status, x = x, y = y,type="thpl")
    #thplvfe.out@dos <- dos
    thplvfe.outgradient <- thpl.out@gradient
	if(!alt)
		thplvfe.out@parm <- 1
	else thplvfe.out@parm <- 2

    thplvfe.out@m <- m
    thplvfe.out@cv <- cv

    ##    computation of assay performance measures here
    if(ap.calc) {
        thplvfe.out@mdc <- mdc.calc(thplvfe.out, conf = conf,type="thpl")$mdc
        thplvfe.out@rdl <- rdl.calc(thplvfe.out, conf = conf,type="thpl")$rdl
		
		if(thplvfe.out@rdl > max(x))
			thplvfe.out@rdlwarn <- "RDL outside range of values"
		else
			thplvfe.out@rdlwarn <- "RDL within range of values"
		
        thplvfe.out@loq <- loq.calc(thplvfe.out,type="thpl")$loq
    }
    else {
        thplvfe.out@mdc <- NULL
        thplvfe.out@rdl <- NULL
        thplvfe.out@loq <- NULL
    }
	## Lack of fit computations here
    
    if(lof.calc) {
         if(length(x) > length(unique(x)))
            thplvfe.out@lof.test <- lof.test(thplvfe.out)
    #   else thplvfe.out@lof.test <- NULL
       }
    thplvfe.out
}

lin.pom <- function(x, y, theta, method = "squared", m = 2, cv = 0.2, conf = 0.95, 
	theta.init = 1, mmod = "lin", intercept = T, kmax = 75, mx = 25,
	ap.calc = T, lof.calc = T)#, dos = T)
## 02-08-07: Removed print from options. DVS
{
#################################################################
## Fits a linear model rather than FPL - very similar code to  ##
##  to fpl.pom.  Most inputs are the same                      ##
## Inputs in addition to those of fpl.pom                      ##
##  mmod: indicate whether you want to fit a line = "lin"       ##
##        or a quadratic = "quad"                              ##
##  intercept: T or F depending on if an intercept is to be    ##
##             included                                        ##
## default kmax is higher here than in fpl                     ##
##                                                             ##
## comments added by Matthew Mitchell on 3/18 or 3/26/03 are   ##
##  indicated with ##, original with ###                       ##
## most of the changes made by MM involved changing SPLUS      ##
## commands to R commands                                      ##
##                                                             ##
## changed the default kmax to 75 on 3/28/03                   ##
#################################################################

    ## 02-08-07: Removed print from options. DVS
	#glsvfe.print <- function(k, coef, reglsd, theta, digits = c(1, 
	#	2, 3, 4, 4, 3))
	#{
	#	cat(paste("\n k=", k, "\n"))
	#	cat(paste("\t\t\tEstimate\tStd.Error\n"))
	#	for(i in seq(along = coef))
	#		cat(paste("\tb", i, ":", "\t\t", format(round(
	#			coef[i], digits[1])), "\t\t", format(
	#			round(reglsd$std.err[i], digits[2])), 
	#			"\n", sep = ""))
	#	cat(paste("\n\tsigma=", format(round(reglsd$std.dev, 
	#		digits[5])), "\n"))
	#	if(!missing(theta))
	#		cat(paste("\ttheta=", format(round(theta, 
	#			digits[6])), "\n"))
	#}
#################################################################################
	m.target <- c("absolute", "squared", "reml", "")
	m.ind <- charmatch(method, m.target)
	method <- m.target[m.ind]
	x <- x[is.finite(y)]
	y <- y[is.finite(y)]
	if(missing(theta))
		theta.calc <- T
	else theta.calc <- F
	k <- 0
	if(mmod == "lin")
		regls <- lsfit(x, y, int = intercept)
    
    if(mmod == "quad")
		regls <- lsfit(cbind(x, x * x), y, int = intercept)
	if(mmod != "lin" & mmod != "quad")
		stop("\n invalid mean model\n")
	reglsd <- ls.diag(regls)
	coef <- regls$coef
	## 03-01-07 Took out the following line. No longer printing out results. DVS
	##if(print) glsvfe.print(k, coef, reglsd)	## now iterate
	k <- 1
	theta.old <- theta.init
	std.dev.old <- reglsd$std.dev
        
	##	repeat {   loop added by MM on 3/26/03

	### Step 2.
	### Given B(0) use (17) to obtain theta(0) and form estimated weights
	### use <<- and synchronize to make sure that pred and res have their
	### current values in frame0
	###

	## MM 3/26/03 frame is an argument for SPLUS but not R
	## made changes similar to that of fpl.pom
               
    for(k in 1:kmax){
		pred <- abs(y - regls$residuals)
		res <- regls$residuals
        if(mmod == "lin")
        	grad <- cbind(1,x)
        if(mmod == "quad")
        	grad <- cbind(1,x,x*x)
  	
	###  if theta is supplied 
	###  vfemin.pom is power of the mean variance function
	###	 call vfemin.abs.pom() if method='abs'
	###       reml.pom if method='reml'
	###	     otherwise vfemin.pom() for pl squared residuals 
	## 02-08-07: Since none of the C functions are in use removing code below
	#		if(!dos) {
	#			if(!is.loaded(symbol.C("linlow"))) {
	#				cat("Loading linpom.o\n")
	#				dyn.load(
	#				  "/local/lib/calib/sun4/src/linpom.o")
	#			}
	#		}
	## nlmin is not found in R so these lines were modified
	## adjustments followed the code for fpl.pom
	## done by MM on 3/26/03
		if(theta.calc) {
			if(method == "absolute")
            	nl.out <- optim(theta.init, vfemin.abs.pom, method="BFGS") 
			if(method == "reml")
				#{
	            #if(dos==F) 04-17-07 Removing dos argument
                # 	nl.out <- optim(theta.init, reml.pom, method="BFGS")
                #else 
                nl.out <- optim(theta.init, reml.pom, method="BFGS",max.facl=60, maxiter=30)
                #}
			if(method == "squared")
            	nl.out <- optim(theta.init, vfemin.pom, method="BFGS", pred=pred, res=res)
			## MM 3/27/03 changed line  if(nl.out$converged)
            if(nl.out$convergence==0)
				theta <- nl.out$par  ## $x changed to $par MM 3/27
			else {
				cat(paste("theta failed to converge:", 
					nl.out$conv.type, "\n"))
				return()
			}
		}
	### Step 3.
	### now do a weighted least squares regression with POM weights
	### and estimated theta

		    ## 04-03-07: So making a few changes here. It seems as though that throughout
            ## the calib program we have situations like that below. However it seems to
            ## me, at least based off of the theory that this is developed off of that the
            ## weights should be strictly greater than zero. However by having 2 * theta 
            ## like below allows for the possibility of non-positive weights. Changing
            ## from 1/pred^(2*theta) to what is now below. DVS
            wts <- 1/((pred^2)^theta)
		if(mmod == "lin")
			regls <- lsfit(x, y, wt = wts, int = intercept)
               
		if(mmod == "quad")
			regls <- lsfit(cbind(x, x * x), y, wt = wts, 
				int = intercept)
		reglsd <- ls.diag(regls)
		coef <- regls$coef
		std.dev <- reglsd$std.dev
		cov.mat <- reglsd$cov.unscaled
		df.error <- length(x) - length(coef)
		## 03-01-07: Taking out the following line. DVS
		#if(print) glsvfe.print(k, coef, reglsd, theta)	
		### now check for convergence or exceeding max number of iterations
		if(theta.calc) {
			if(abs(theta - theta.old) < 0.0001) {
				status <- "converged"
				#cat("\n\n Theta converged\n\n")
				break
			}
			else if(k >= kmax) {
				status <- "failed to converge"
				#cat("\n\n Theta failed to converge - try increasing kmax\n\n")
				break
			}
			k <- k + 1
			theta.old <- theta
		}
		else {
			if(abs(std.dev - std.dev.old) < 0.0001 * 
				std.dev.old) {
				status <- "converged"
				#cat("\n\n Standard Deviation converged\n\n")
				break
			}
		else if(k >= kmax) {
			status <- "failed to converge"
			#cat("\n\n Standard Deviation failed to converge - 
			#	try increasing kmax\n\n")
				break
			}
			k <- k + 1         
			std.dev.old <- std.dev
		}
	}
	if(mmod == "quad") {
		if(sign(coef[2]) != sign(coef[3])) {
			tp <- ( - coef[2]/(2 * coef[3]))
			if(tp <= max(x))
				stop("calibration not valid with this quadratic model")
		}
	}
	if(theta.calc) {
		if(method == "absolute")
			vfe.method <- "pl absolute residuals"
		if(method == "reml")
			vfe.method <- "reml"
		if(method == "squared")
			vfe.method <- "pl squared residuals"
	}
	## originally the next line was else but kept getting a syntax error message
	## so changed to the next line  MM 3/28/03
	if(!theta.calc) {
		if(theta == 0)
			vfe.method <- "OLS"
		else vfe.method <- "theta fixed"
	}

	## 9/26/03 need output of if OLS or POM was used
    if(theta.calc){var.model="POM"}
    if(theta==0){var.model="OLS"}


	## 9/26/03 add usual linear model statistics (SSR, SSE, SST, Rsq, ...)
    SST <- sum((y - mean(y))^2)
    SSE <- sum((regls$residuals)^2)
    SSR <- SST - SSE
    Rsq <- 100*(SSR/SST)

	pom <- new("calib.fit",coefficients = coef, se.coefficients = as.vector(
		std.dev * sqrt(diag(cov.mat))), sigma = std.dev, theta
		 = theta, cov.unscaled = cov.mat, df.residual = 
		df.error, fitted.values = y - regls$residuals, 
		residuals = regls$residuals, vfe.method = vfe.method, 
		kused = k, status = status, m = m, cv = cv, x = x, y = 
		y, var.model=var.model, SST=SST, SSE=SSE, SSR=SSR, Rsq=Rsq,type="lin")

        pom@mmod <- mmod
        #pom@dos <- dos
        
	if(ap.calc) {
		pom@mdc <- mdc.calc(pom, conf = conf,type="lin")$mdc
		pom@rdl <- rdl.calc(pom, conf = conf,type="lin")$rdl
		
		if(pom@rdl > max(x))
			pom@rdlwarn <- "RDL outside range of values"
		else
			pom@rdlwarn <- "RDL within range of values"
		
		pom@cv <- cv
		pom@loq <- loq.calc(pom,type="lin")$loq
	}
	if(lof.calc) {
		if(length(x) > length(unique(x)) & length(unique(x)) != 
			length(coef))
			pom@lof.test <- lof.test(pom)
			
	## 01-25-07 Had made this change some while back. It does not seem
	## necessary with the S4 structure to specify that the lof.test
	## be set to NULL since this slot can just as easily remain empty. DVS
	
	#	else attr(pom, "lof.test") <- NULL
	}
	pom
}

fplpom <- function(x, y, wts, theta, b3start, b4start, mx, alt = F)
{
#########################################################################
## FPL POM function.                                                   ##
## This function does the weighted nls fit with the POM model.         ##
## This function gets called iteratively in the fpl.pom function       ##
##  in order to receive the final estimates of theta,b1,b2,b3,b4.      ##
##                                                                     ##
## modified by Matthew Mitchell on 3/11/03                             ##
##                                                                     ##
##  x: the independent variable (usually dosage here)                  ##
##  y: the response variable (example: OD, optical density)            ##
##  wts: the weights used in the weighted nls                          ##
##  theta: the current value of the POM parameter                      ##
##  mx: the maximum number of iterations in the nls procedure          ## 
##  alt:  If alt=F the nonlinear function fitted is                    ##
##        y = b2 + (trb1-b2)/(1 + (x/b3)^b4)                           ##
##        If alt=T the nonlinear function fitted is                    ##
##        y = b2 + (b1-b2)/(1 + exp(b4(logx - b3)),                    ##       
##        which is reparametrization of the other form                 ##
##                                                                     ##
## output variables:                                                   ##
##   name: "Output from fplpom"                                        ##
##  coefficients: the estimated values for b1,b2,b3,b4                 ##
##  se.coefficients: the standard errors of these                      ##
##  sigma: the residual standard error                                 ##
##  cov.unscaled: this is the variance-covariance matrix               ##
##   of the estimates of b1,b2,b3,b4 without the sigma^2               ##
##  df.residual: the residual degrees of freedom                       ##
##  residual: all the residuals from the fit                           ##
##  fitted: the fitted values from the fit                             ##
##  x,y: x and y values that were input                                ##
##                                                                     ##
## Original comments are ###, added are ##                             ##
#########################################################################
#browser()
## line added by MM on 4/15/03 to prevent taking log of 0
   x[x == 0] <- 0.001
   
   if(missing(mx))
   mx<-50
   else
   mx<-mx
   
    wy <- as.vector(sqrt(wts) * y)
    fpl.frame <- data.frame(x = x, y = y, w = sqrt(wts), wy = wy)

    if(missing(b4start))
        b4start <- 1
    if(missing(b3start)) {
        ymid <- (mean(y[x == min(x)]) + mean(y[x == max(x)]))/2
        y.diff <- abs(y - ymid)
        min.diff <- min(y.diff)
        ymid <- (y[y.diff == min.diff])[1]
        xmid <- (x[y.diff == min.diff])[1]
        b3start <- xmid/(((mean(y[x == min(x)]) - mean(y[x == 
            max(x)]))/(ymid - mean(y[x == max(x)])) - 1)^(1/
            b4start))
        if(alt)
            b3start <- log(b3start)
    }

    attr(fpl.frame, "param") <- list(b1=mean(y[x == min(x)]), b2=mean(y[x == max(x)]),
                                         b3=b3start, b4=b4start)

	## 01-23-07 Note that this was changed from fpl.frame$ptype <- ... to below. Having some
	## troubles getting this working at the moment.
	attr(fpl.frame, "ptype") <- ifelse(alt == T, 2, 1)
	
	## 01-23-07 changed above to be an attribute as otherwise in the logical statement above
	## ptype becomes a vector in order to conform to the dimension of the data frame
	## fpl.frame and error message is produced, however still want to keep fpl.frame
	## an object of class data.frame for the purposes of the non-linear fitting.
	
    nls.out <- nls(wy ~ fpl.model(x, b1, b2, b3, b4, w, parm = attr(fpl.frame, "ptype")),
                 fpl.frame, control = nls.control(maxiter = mx),
                 start=attr(fpl.frame, "param"))
                 
    b <- coef(nls.out)
    n <- length(x)
    fitted <- fpl.model(x, b, parm = attr(fpl.frame,"ptype"))
    gradient <- attributes(fitted)$gradient
    ## 04-03-07: So making a few changes here. It seems as though that throughout
    ## the calib program we have situations like that below. However it seems to
    ## me, at least based off of the theory that this is developed off of that the
    ## weights should be strictly greater than zero. However by having 2 * theta 
    ## like below allows for the possibility of non-positive weights. Changing
    ## from 1/pred^(2*theta) to what is now below. DVS
    wts <- 1/((as.vector(fitted)^2)^theta)
# Was this 04-03-07 wts <- (1/as.vector(fitted)^theta) * (1/as.vector(fitted)^theta)
    resid <- y - as.vector(fitted)
## The next set of lines are computing the std errors of the coefficients
    z <- resid * resid * wts
    sumz <- sum(z)
    sigma <- sqrt(sumz/(n - 4))
    g.mat <- diag(wts)
    sg <- (t(attributes(fitted)$gradient)) %*% g.mat %*% (
        attributes(fitted)$gradient)
    cov.un <- solve(sg)
## These are the output values used in fpl.pom for the next iteration
    fpl.out <- new("calib.fit", coefficients = b, 
        se.coefficients = as.vector(sigma * sqrt(diag(cov.un))),
        sigma = sigma, cov.unscaled = cov.un, df.residual = n - 4, 
        residuals = resid, fitted.values = as.vector(fitted), x = x, y
         = y)
    fpl.out@gradient <- gradient
    fpl.out
}

thplpom <- function(x, y, wts, theta, b3start, mx, alt = F)
{
#########################################################################
## THPL POM function.                                                  ##
## This function does the weighted nls fit with the POM model.         ##
## This function gets called iteratively in the thpl.pom function      ##
##  in order to receive the final estimates of theta,b1,b2,b3.         ##
##                                                                     ##
## written  by Matthew Mitchell on 5/20/03                             ##
##                                                                     ##
##  x: the independent variable (usually glucose conc)                 ##
##  y: the response variable (example: intensity)                      ##
##  wts: the weights used in the weighted nls                          ##
##  theta: the current value of the POM parameter                      ##
##  mx: the maximum number of iterations in the nls procedure          ## 
##                                                                     ##
## output variables:                                                   ##
##   name: "Output from thplpom"                                       ##
##  coefficients: the estimated values for b1,b2,b3                    ##
##  se.coefficients: the standard errors of these                      ##
##  sigma: the residual standard error                                 ##
##  cov.unscaled: this is the variance-covariance matrix               ##
##   of the estimates of b1,b2,b3    without the sigma^2               ##
##  df.residual: the residual degrees of freedom                       ##
##  residual: all the residuals from the fit                           ##
##  fitted: the fitted values from the fit                             ##
##  x,y: x and y values that were input                                ##
##                                                                     ##
## Original comments are ###, added are ##                             ##
#########################################################################
	
	if(missing(mx))
	mx<-25
	else
	mx<-mx
	
    wy <- as.vector(sqrt(wts) * y)
    thpl.frame <- data.frame(x = x, y = y, w = sqrt(wts), wy = wy)

    if(missing(b3start)) {
        ymid <- (mean(y[x == min(x)]) + mean(y[x == max(x)]))/2
        y.diff <- abs(y - ymid)
        min.diff <- min(y.diff)
        ymid <- (y[y.diff == min.diff])[1]
        xmid <- (x[y.diff == min.diff])[1]
        b3start <- xmid/(((mean(y[x == min(x)]) - mean(y[x == 
            max(x)]))/(ymid - mean(y[x == max(x)])) - 1))
		if(alt)
			b3start <- log(b3start)
          }

    attr(thpl.frame, "param") <- list(b1=mean(y[x == min(x)]), b2=mean(y[x == max(x)]),
                                         b3=b3start)
	
	attr(thpl.frame, "ptype") <- ifelse(alt == T, 2, 1)

    nls.out <- nls(wy ~ thpl.model(x, b1, b2, b3, w, parm = attr(thpl.frame, "ptype")),
                 thpl.frame, control = nls.control(maxiter = mx),
                 start=attr(thpl.frame, "param"))
                 
    b <- coef(nls.out)
    n <- length(x)
    fitted <- thpl.model(x, b, parm = attr(thpl.frame, "ptype"))
    gradient <- attributes(fitted)$gradient
    ## 04-03-07: So making a few changes here. It seems as though that throughout
    ## the calib program we have situations like that below. However it seems to
    ## me, at least based off of the theory that this is developed off of that the
    ## weights should be strictly greater than zero. However by having 2 * theta 
    ## like below allows for the possibility of non-positive weights. Changing
    ## from 1/pred^(2*theta) to what is now below. DVS
    wts <- 1/((as.vector(fitted)^2)^theta)
# Was this 04-03-07 DVS    wts <- (1/as.vector(fitted)^theta) * (1/as.vector(fitted)^theta)
    resid <- y - as.vector(fitted)
## The next set of lines are computing the std errors of the coefficients
    z <- resid * resid * wts
    sumz <- sum(z)
    sigma <- sqrt(sumz/(n - 3))
    g.mat <- diag(wts)
    sg <- (t(attributes(fitted)$gradient)) %*% g.mat %*% (
        attributes(fitted)$gradient)
    cov.un <- solve(sg)
## These are the output values used in thpl.pom for the next iteration
    thpl.out <- new("calib.fit", coefficients = b, 
        se.coefficients = as.vector(sigma * sqrt(diag(cov.un))),
        sigma = sigma, cov.unscaled = cov.un, df.residual = n - 3, 
        residuals = resid, fitted.values = as.vector(fitted), x = x, y
         = y)
    thpl.out@gradient <- gradient
    thpl.out
}



fpl.model <- function(x, b1, b2, b3, b4, w = 1, parm = 1){
####################################################################
## The 4 Parameter Logistic Model used in Calibration Analysis    ##
## modified by Matthew Mitchell on 3/11/03                        ##
##                                                                ##
## input variables:                                               ##
##  x: the independent variable (usually dosage here)             ##
##  b1,b2,b3,b4: values used in the FPL fit                       ##
##  w: weights.  In OLS these are 1.                              ##
##  parm:  If parm=1 the nonlinear function fitted is             ##
##        y = b2 + (b1-b2)/(1 + (x/b3)^b4)                        ##
##        If parm=2 the nonlinear function fitted is              ##
##        y = b2 + (b1-b2)/(1 + exp(b4(logx - b3)),               ##       
##        which is reparametrization of the other form            ##
##        since b3 from alt=T is log(b3 from alt=F).              ##
##                                                                ##
## This function is called when doing the nls fitting with        ##
##  the fpl an dfpl.pom functions.                                ##
##                                                                ## 
## Original comments are ###, added are ##                        ##
####################################################################
  
	if(length(b1) == 4) {
		b2 <- b1[2]
		b3 <- b1[3]
		b4 <- b1[4]
		b1 <- b1[1]
	}
	if(parm == 1) {
		.a <- ifelse(x <= 0, 0, (x/b3)^b4)
		.den <- 1 + .a

		.value <- w * ((b1 - b2)/.den + b2)
	}
	else {
		.a <- ifelse(x <= 0, 0, exp(b4 * (log(x) - b3)))
		.den <- 1 + .a
		.value <- w * ((b1 - b2)/.den + b2)
	}
	.grad <- array(.value, c(length(.value), 4), list(NULL, c("b1", 
		"b2", "b3", "b4"))) ## rownames are NULL, columnames are
                                    ## are the parameter names
	.grad[, "b1"] <- w * (1/.den)
	.grad[, "b2"] <- w * (1 - 1/.den)
	if(parm == 1) {

		.grad[, "b3"] <- w * (((b1 - b2) * (b4/b3) * .a)/.den^2
			)
		.grad[, "b4"] <- ifelse(x == 0, 0, (w * ( - (b1 - b2) * 
			.a * log(ifelse(x > 0, x/b3, NA))))/.den^2)
	}
	else {
		.grad[, "b3"] <- w * (((b1 - b2) * b4 * .a)/.den^2)
		.grad[, "b4"] <- ifelse(x == 0, 0, w * ((( - (b1 - b2) * 
			.a * (log(ifelse(x > 0, x, NA)) - b3))/.den^2))
			)
	}
        attr(.value, "gradient") <- .grad
	.value
}


## This is the old thpl.model function which does not have the log parameterization
#thpl.model <- function(x, b1, b2, b3, w = 1)
#{
#####################################################################
### The 3 Parameter Logistic Model used in Calibration Analysis    ##
### written  by Matthew Mitchell on 5/20/03                        ##
###                                                                ##
### input variables:                                               ##
###  x: the independent variable (usually dosage here)             ##
###  b1,b2,b3: values used in the THPL fit                         ##
###  w: weights.  In OLS these are 1.                              ##
###                                                                ##
### This function is called when doing the nls fitting with        ##
###  the thpl an  thpl.pom functions.                              ##
###                                                                ## 
#####################################################################
#  
#	if(length(b1) == 3) {
#		b2 <- b1[2]
#		b3 <- b1[3]
#	        b1 <- b1[1] }
#
#		.a <- ifelse(x <= 0, 0, (x/b3))
#		.den <- 1 + .a
#		.value <- w * ((b1 - b2)/.den + b2)
#
#	.grad <- array(.value, c(length(.value), 3), list(NULL, c("b1", 
#		"b2", "b3"))) 
#	## rownames are NULL, columnames are are the parameter names
#	.grad[, "b1"] <- w * (1/.den)
#	.grad[, "b2"] <- w * (1 - 1/.den)
#	.grad[, "b3"] <- w * ((b1 - b2) * (.a/b3))/.den^2
#		
#	attr(.value, "gradient") <- .grad
#	.value
#}

## This is the new thpl.model that does have the log parameterization

thpl.model <- function(x, b1, b2, b3, w = 1, parm = 1)
{
  
	if(length(b1) == 3) {
		b2 <- b1[2]
		b3 <- b1[3]
	    b1 <- b1[1] 
	}
	if(parm == 1) {
		.a <- ifelse(x <= 0, 0, (x/b3))
		.den <- 1 + .a
		.value <- w * ((b1 - b2)/.den + b2)
	}
	else {
		.a <- ifelse(x <= 0, 0, exp(log(x) - b3))
		.den <- 1 + .a
		.value <- w * ((b1 - b2)/.den + b2)
	}
	.grad <- array(.value, c(length(.value), 3), list(NULL, c("b1", 
		"b2", "b3"))) 
	## rownames are NULL, columnames are are the parameter names
	.grad[, "b1"] <- w * (1/.den)
	.grad[, "b2"] <- w * (1 - 1/.den)
	.grad[, "b3"] <- w * ((b1 - b2) * (.a/b3))/.den^2
		
	attr(.value, "gradient") <- .grad
	.value
}

lin.model <- function(x, beta, w = 1, mmod)
{
##################################################
## called for the linear model fit in lin.pom   ##
##                                              ##
## header added by Matthew Mitchell on 3/18/03  ##
## actually I don't this is used? 3/28/03 MM    ##
##################################################
  
	if(mmod == "lin") {
		.value <- beta[1] + beta[2] * x
		.grad <- array(.value, c(length(.value), 2), list(NULL, 
			c("b1", "b2")))
		.grad[, "b1"] <- w
		.grad[, "b2"] <- w * x
	}
        
	if(mmod == "quad") {
		.value <- beta[1] + beta[2] * x + beta[3] * x * x
		.grad <- array(.value, c(length(.value), 3), list(NULL, 
			c("b1", "b2", "b3")))
		.grad[, "b1"] <- w
		.grad[, "b2"] <- w * x
		.grad[, "b3"] <- w * x * x
	}
	attr(.value, "gradient") <- .grad
	.value
}



vfemin.abs.pom <- function(theta,pred,res)
{
### This function gets called as part of the fpl.pom function
### pdh: 4/23/92 - This function is minimized to find the best value for
### theta as part of the pseudolikelihood approach
### See Eq. (18) of Davidian and Haaland (1990)
### Note that is function uses absolute values and should be more robust
### against possible outliers

	g <- exp(theta * log(pred))
	gdot <- exp(mean(log(g)))
	gdot * sum(abs(res)/g)
}



## This is used to estimate the POM parameter theta in fpl.pom
 
vfemin.pom <- function(theta,pred,res)
{
    g <- sqrt((pred^2)^theta)
    gdot <- exp(mean(log(g)))
    (gdot^2) * sum((res/g)^2)
}

## This function is used to estimate theta

reml.pom <- function(theta,pred)
{
	g <- pred^theta
	gmat.inv <- diag(1/(g * g))
	hmat <- t(grad) %*% gmat.inv %*% grad
	lh <- length(hmat[, 1])
	determ.h <- prod(svd(hmat)$d)
	gdot <- exp(mean(log(g)))
	n <- length(g)
	out <- (gdot^((2 * n)/(n - lh))) * exp((1/(n - lh)) * log(
		determ.h)) * sum((res/g)^2)
	out
}

