glsvfe.print <- function(k, b, se.b, sigma, theta, digits, mod=4, bios=T)
{
################################################################
## Generalized Least Squares and Variance Function Estimation ##
## modified by Matthew Mitchell on 3/12/03                    ##
##                                                            ##
## This function is called when running fpl.pom (or thol.pom) ##
##  and prints                                                ##
##  out the results from the gls and vfe from each iteration. ##
##                                                            ##
## input:                                                     ##
##  k: the kth iteration                                      ##
##  b: the vector of current estimates for b1,b2,b3,b4        ##
##  se.b: the standard error of the estimates                 ##
##  sigma: the residual standard error                        ##
##  theta: current estimate of the POM parameter              ##
##  digits: used to determine the number of significant       ##
##          digits.  It is set to 3 in fpl.pom.               ##
##  mod: 4 if fpl is being used, 3 for thpl                   ##
##                                                            ##
##  The comments I added have ##, original ###                ##
##                                                            ##
## 5/20/03 MM modified so that this prints with new           ##
## three parameter logistic models                            ##
## 5/21/03 added bios (Biosense) option to print F0, Finf,Kd  ##
## for the 3PL if desired                                     ##
################################################################

  
## for 1 < num <=10, vs=1 10 < num <= 100, vs =2 , etc.
## vec.signif is used to determine the number of digits
##  for the residual standard error
	vec.signif <- function(x, digits)
	{
		num <- max(abs(x))
		vs <- ceiling(log(num)/log(10))
		max(digits - vs, 0)
	}
	std.err <- as.vector(se.b)
	coefb <- as.vector(b)

## Do these two lines get used?
	rdigits1 <- vec.signif(coefb, digits)
	rdigits2 <- vec.signif(std.err, digits)

## This is the start of the printout
	cat(paste("\n k=", k, "\n"))
## second argument of round is 3 here, whihc is the thousandth place
##  if it were 1, it would be to the nearest tenth,
##  0 to the nearest one, +1 to the nearest 10
        if(mod == 4){
	print(data.frame(Estimate = round(coefb, digits), Std.Error = 
		round(std.err, digits), row.names = c("B1", "B2", "B3", 
		"B4")))}

        if(mod == 3){
           if(!bios) {
           print(data.frame(Estimate = round(coefb, digits), Std.Error = 
		round(std.err, digits), row.names = c("B1", "B2", "B3")))}
           if(bios) {
           print(data.frame(Estimate = round(coefb, digits), Std.Error = 
		round(std.err, digits), row.names = c("F0", "Finf", "Kd")))}
         }

        
	cat(paste("\tsigma=", format(round(sigma, vec.signif(sigma, 
		digits))), "\n"))
	if(!missing(theta))
		cat(paste("\ttheta=", format(round(theta, 2)), "\n"))
	cat(paste("\n")) 
       
      }


summary.calib <- function(calout, digits = 2)#, invprint=T)
{
##################################################################################
## prints calib objects                                                         ##
## modified by Matthew Mitchell on 3/24/03, 3/27                                ##
##                                                                              ##
## The invprint argument was added by MM because when m, the number of a        ##
## reps is a vector with unequal entries the inverse limits are computed using  ##
## the first element of m only.                                                 ##
##                                                                              ##
## On 3/27/03 code was added to flag those samples whose estimated              ##
## calibrations are less than the mdc.                                          ##
## modified on 8/30/04           												##
## 06-08-07: Made a bunch of changes to this function, now serves as a summary  ##
## method.
##################################################################################  


	x <- calout@Estimated.x
	se.calib <- calout@PredStdErr
	n.x <- length(x)

	cname <- calout@labels$name
	## Next 2 lines were added on 3/27/03 by MM so that "avg rname" (e.g., "avg OD") 
	## will show on the printout rather than "response"
    rname <- calout@rname

	#eval(parse(text=paste("calout@","avg.",rname,sep="")))

    #rname <- paste("avg",rname)
    oor <- ifelse(calout@oor == "x", "x", "")
    oor.ind <- ifelse(oor == "x", 1, 0)

 	## These two lines added by MM on 3/27/03, which are used to
	## indicate those estimated calibrations less than the mdc
    lmdc <- ifelse(calout@lmdc == "X", "x", "")
    lmdc.ind <- ifelse(lmdc == "x", 1, 0)

#    cat(paste("\n\n\ncalibration output\n\n approx. confidence level for limits: ",
#        format(round(calout@conf.level, 2)), "\n\n"
#        ))
#    if(invprint){
#    if(length(calout@samp.names)>0) {
#        snames <- calout@samp.names
#        cal.mat <- matrix(cbind(snames, round(x, 
#            digits), round(se.calib, digits), round(
#            calout@inver.low, digits), round(calout@
#            inver.up, digits), round(calout@wald.low, 
#            digits), round(calout@wald.up, digits), round(
#            calout@response, digits), oor, lmdc), nrow = length(
#            x), ncol = 10, dimnames = list(rep("",
#            length(calout@calib)), c(cname, "calib", 
#            "se.calib", "inver.low", "inver.up", "wald.low",
#            "wald.up", rname, "oor", "lmdc")))
#    }
#    else {
    if(length(calout@samp.names)>0)
        samp.names <- calout@samp.names
    else
    	samp.names <- 1:n.x
    	
    cal.mat <- data.frame(round(x, digits),
    	round(se.calib, digits),
    	round(calout@inver.low, digits),
    	round(calout@inver.up, digits),
    	round(calout@wald.low, digits),
    	round(calout@wald.up, digits),
    	round(calout@avg.response, digits))
    
    	names(cal.mat) <- c("Predicted X", "Std Err X", "Low CI: Inv", 
            "Up CI: Inv", "Low CI: Wald", "Up CI: Wald", "Y")
        rownames(cal.mat) <- samp.names
        
#    		
#        cal.mat <- data.frame(cbind(round(x, digits), 
#            round(se.calib, digits), round(calout@
#            inver.low, digits), round(calout@inver.up, 
#            digits), round(calout@wald.low, digits), round(
#            calout@wald.up, digits), round(calout@avg.response, 
#            digits), oor, lmdc), nrow = length(x), 
#            ncol = 9, dimnames = list(rep("", length(x
#            )), c("calib", "se.calib", "inver.low", 
#            "inver.up", "wald.low", "wald.up", rname, 
#            "oor", "lmdc")))
#    }
#    }
#     if(!invprint){
#    if(length(calout@samp.names)>0) {
#        snames <- calout@samp.names
#        cal.mat <- matrix(cbind(snames, round(x, 
#            digits), round(se.calib, digits), round(calout@wald.low, 
#            digits), round(calout@wald.up, digits),
#            round(calout@avg.response, digits), oor, lmdc), nrow =
#            length(x), ncol = 8, dimnames = list(rep("",
#            length(x)), c(cname, "calib", 
#            "se.calib",  "wald.low",
#            "wald.up", rname, "oor", "lmdc")))
#    }
#    else {
#        cal.mat <- matrix(cbind(round(x, digits), 
#            round(se.calib, digits), round(calout@wald.low, digits),
#            round(calout@wald.up, digits), round(x, 
#            digits), oor, lmdc), nrow = length(x), 
#            ncol = 7, dimnames = list(rep("", length(x
#            )), c("calib", "se.calib", "wald.low", "wald.up", rname, 
#            "oor", "lmdc")))
#    }}
#    if(sum(oor.ind) == 0) {
#        cal.mat <- cal.mat[,  - (ncol(cal.mat))+ 1]
# ## MM 3/27/03 changed this line to print a comment rather than blank and +1 added above
#        post.mess <- "\nNo responses are out of range."
#    }
#    else post.mess <- "\noor='x' indicates response out of range"

 ## next loop added by MM on 3/27/03
#    if(sum(lmdc.ind, na.rm=T) == 0){
#        cal.mat <- cal.mat[ , - (ncol(cal.mat))] 
#        post.stuff <- "\nAll estimated calibrations are above the mdc.\n\n"
#      }
#    else post.stuff <- "\nlmdc='x' indicates that the estimated calibration is less than the mdc\n\n"
    
    #nameHold <- dimnames(cal.mat)
    cal.mat <- format(cal.mat, justify="centre")
    #dimnames(cal.mat) <- nameHold
    #colnames(cal.mat)[2:3] <- slotNames(calout)[1:2]   ## MWM 8/30/04
    #print(cal.mat, quote = F)
    #cat(post.mess)
## added by MM on 3/27/03
    #cat(post.stuff)
    #invisible()
    cal.mat
}

print.fpl.pom <- function(fpl.out, digits = 3)
{
	vec.signif <- function(x, digs)
	{
		num <- max(abs(x))
		vs <- ceiling(log(num)/log(10))
		max(digs - vs, 0)
	}
#  cat(paste('\n\t\t\tEstimate\tStd.Error\n'))
#  for(i in seq(along=fpl.out$coef)) 
#    cat(paste('\tb',i,':', 
#	      '\t\t',format(round(fpl.out$coef[i],digits)), 
#	      '\t\t',format(round(fpl.out$se.coef[i],digits)), 
#	      '\n',sep=''))
#  
	coef.df <- data.frame(Parameter = c("B1", "B2", "B3", "B4"), 
		Estimate = format(round(fpl.out@coefficients, digits)), 
		Std.Error = format(round(fpl.out@se.coefficients, digits)))
	row.names(coef.df) <- c("", " ", "   ", "    ")
	cat("\n")
	cat(paste("  Regression output:\n\n"))
	print(coef.df, right = F)	
	#  cat(paste("\n\n\tsigma=",format(round(fpl.out$sigma,vec.signif(fpl.out$sigma,digits)))))
	cat(paste("\n\n  sigma=", format(round(fpl.out@sigma, 
		vec.signif(fpl.out@sigma, digits)))))
	if(fpl.out@vfe.method != "OLS") {
		cat(paste("\ttheta=", format(round(fpl.out@theta, 2)), 
			"\n"))
		cat(paste("  Theta estimated by method of:", fpl.out@
			vfe.method, "\n"))
		if(fpl.out@status == "converged")
			cat(paste("  Iteration converged at k =", 
				fpl.out@kused))
		else cat(paste("  Iteration failed to converge for kmax =",
				fpl.out@kused))
	}
	cat("\n\n")	# note: i have take first tabs out below
	if(length(fpl.out@mdc)>0) {
		ap.vec <- c(fpl.out@mdc, fpl.out@rdl, fpl.out@loq)
		ap.vec <- ap.vec[is.finite(ap.vec)]
		cat(paste("\n  m=", format(round(fpl.out@m, 1))))
		rdigits3 <- vec.signif(ap.vec, 4)
		cat(paste("\tMDC=", format(round(fpl.out@mdc, rdigits3)
			), "\n"))
		cat(paste("\tRDL=", format(round(fpl.out@rdl, rdigits3)
			), "\n"))
		if(!is.na(fpl.out@loq))
			cat(paste("\tLOQ (c.v.= ", format(round(fpl.out@
				cv, 2)), ")= ", format(round(fpl.out@
				loq, rdigits3)), "\n", sep = ""))
		else cat(paste("\tLOQ not attainable for c.v= ", format(
				round(fpl.out@cv, 2)), "\n", sep = ""))
		cat("\n\n")
	}
        ## Need to think about this logical statement. Previous just did a check for
        ## whether or not the attribute was NULL, however when using slots it does not
        ## come up as being NULL even when there is not information there.
	if(length(fpl.out@lof.test@Fstat)>0) {
		lof.out <- fpl.out@lof.test	
	#      rdigits4 <- vec.signif(c(lof.out$lofss,lof.out$pure.error),digits)
#      cat(paste(" Source\t\tSS\tdf\tMS\tF\tp-value\n"))
#      cat(paste(" LOF\t\t",format(round(lof.out$lofss,rdigits4)), 
#		"",format(round(lof.out$df.lof,0)), 
#		"",format(round(lof.out$lofss/lof.out$df.lof,rdigits4)), 
#		"",format(round(lof.out$F,rdigits4)), 
#		"",format(round(lof.out$p.value,4)),"\n",sep=""))
#      cat(paste(" Pure Error\t",format(round(lof.out$pure.error,rdigits4)), 
#		"",format(round(lof.out$df.pure.error,0)), 
#		"",format(round(lof.out$pure.error/lof.out$df.pure.error,rdigits4)), 
#		"\n\n\n",sep=""))
		lof.df <- data.frame(Source = c("LOF", "Pure Error"), 
			SS = round(c(lof.out@lofss, lof.out@pure.error),
			digits), df = round(c(lof.out@df.lof, lof.out@
			df.pure.error), 0), MS = round(c(lof.out@lofss/
			lof.out@df.lof, lof.out@pure.error/lof.out@
			df.pure.error), digits), Fval = c(round(lof.out@
			Fstat, digits), ""), pvalue = c(round(lof.out@
			p.value, 4), ""))
		row.names(lof.df) <- c("", " ")
		cat(paste("  Lack of Fit test:\n\n"))
		print(lof.df, quote = F, right = F)
	}
	cat("\n")
	return()
}

#############################################################
##                prints thpl objects                      ##
##       written by Matthew Mitchell on 5/20/03            ##
##                                                         ##
##       added bios option on 5/21/03 so that F0, Finf, Kd ##
##       will be printed if desired                        ##
#############################################################

print.thpl.pom <- function(thpl.out, digits = 3, bios=T)
{
	vec.signif <- function(x, digs)
	{
		num <- max(abs(x))
		vs <- ceiling(log(num)/log(10))
		max(digs - vs, 0)
	}

	coef.df <- data.frame(Parameter = c("F0", "Finf", "Kd"), 
		Estimate = format(round(thpl.out@coefficients, digits)), 
		Std.Error = format(round(thpl.out@se.coefficients, digits)))
	row.names(coef.df) <- c("", " ", "  ")
	cat("\n")
	cat(paste("  Regression output:\n\n"))
	print(coef.df, right = F)	

	cat(paste("\n\n  sigma=", format(round(thpl.out@sigma, 
		vec.signif(thpl.out@sigma, digits)))))
	if(thpl.out@vfe.method != "OLS") {
		cat(paste("\ttheta=", format(round(thpl.out@theta, 2)), 
			"\n"))
		cat(paste("  Theta estimated by method of:", thpl.out@
			vfe.method, "\n"))
		if(thpl.out@status == "converged")
			cat(paste("  Iteration converged at k =", 
				thpl.out@kused))
		else cat(paste("  Iteration failed to converge for kmax =",
				thpl.out@kused))
	}
	cat("\n\n")	# note: i have take first tabs out below
	if(!is.null(thpl.out@mdc)) {
		ap.vec <- c(thpl.out@mdc, thpl.out@rdl, thpl.out@loq)
		ap.vec <- ap.vec[is.finite(ap.vec)]
		cat(paste("\n  m=", format(round(thpl.out@m, 1))))
		rdigits3 <- vec.signif(ap.vec, 4)
		cat(paste("\tMDC=", format(round(thpl.out@mdc, rdigits3)
			), "\n"))
		cat(paste("\tRDL=", format(round(thpl.out@rdl, rdigits3)
			), "\n"))
		if(!is.na(thpl.out@loq))
			cat(paste("\tLOQ (c.v.= ", format(round(thpl.out@
				cv, 2)), ")= ", format(round(thpl.out@
				loq, rdigits3)), "\n", sep = ""))
		else cat(paste("\tLOQ not attainable for c.v= ", format(
				round(thpl.out@cv, 2)), "\n", sep = ""))
		cat("\n\n")
	}
	if(length(thpl.out@lof.test@Fstat)>0) {
		lof.out <- thpl.out@lof.test	

		lof.df <- data.frame(Source = c("LOF", "Pure Error"), 
			SS = round(c(lof.out@lofss, lof.out@pure.error),
			digits), df = round(c(lof.out@df.lof, lof.out@
			df.pure.error), 0), MS = round(c(lof.out@lofss/
			lof.out@df.lof, lof.out@pure.error/lof.out@
			df.pure.error), digits), Fval = c(round(lof.out@
			Fstat, digits), ""), pvalue = c(round(lof.out@
			p.value, 4), ""))
		row.names(lof.df) <- c("", " ")
		cat(paste("  Lack of Fit test:\n\n"))
		print(lof.df, quote = F, right = F)
	}
	cat("\n")
	invisible()
}

print.lin.pom <- function(lin.out, digits = c(2, 2, 2, 2, 2))
{
### prints the results from glsvfe.pom() and mdc.pom() in an attractive manner
### mod
### pdh 4/22/92: corrected standard deviation of the coefficients
### pdh 4/23/92: prints method used to estimate theta
### pdh 5/4/92: prints status of the convergence for the pom model
### pdh 5/5/92: include the amplification factor if present
### pdh 7/21/92: takes new pom structure
#	cat(paste("\n", lin.out$type))
#	cat(paste("\n", lin.out$comment))
	cat(paste("\n\t\t\tEstimate\tStd.Error\n"))
	for(i in seq(along = lin.out@coefficients))
		cat(paste("\tb", i, ":", "\t\t", format(signif(lin.out@
			coefficients[i], digits[1])), "\t\t", format(signif(
			lin.out@sigma * sqrt(lin.out@cov.unscaled[i, i]), 
			digits[2])), "\n", sep = ""))
	cat(paste("\n\tsigma=", format(signif(lin.out@sigma, digits[3])
		)))
	cat(paste("\t\ttheta=", format(signif(lin.out@theta, digits[4])
		), "\n"))
	cat(paste("\tTheta estimated by method of:", lin.out@var.model,
		"\n"))
	if(lin.out@status == "converged")
		cat(paste("\tIteration converged at k =", lin.out@kused,
			"\n\n"))
	else cat(paste("\tIteration failed to converge for kmax =", 
			lin.out@kused, "\n\n"))
	cat(paste("\tm=", format(signif(lin.out@m, 1))))
	cat(paste("\tMDC=", format(signif(lin.out@mdc, digits[5])), 
		"\n"))
	cat(paste("\t\tRDL=", format(signif(lin.out@rdl, digits[5])), 
		"\n"))
	cat(paste("\t\tLOQ=", format(signif(lin.out@loq, digits[5])), 
		"\n\n"))
	invisible()
}

