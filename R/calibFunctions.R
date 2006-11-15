calib <- function(calib.fit.out,y0,conf=.9,dilution=1,samp.names,
	truth,times,samp.units="", dose.units="",dose.name="",maxit=1000,
	toler=1e-05,rname="response", extrap=F,xname=x){
	## 06-08-07: Removing the maxit and toler options. Both of these were used
	## previously in the source call to C but since this is not being used
	## anymore makes sense to take it out. DVS

	if(calib.fit.out@type=="lin"){
		calib.out <- calib.lin.pom(calib.fit.out, y0=y0, conf = conf, 
    		m = calib.fit.out@m, dilution = dilution, samp.names, truth, 
		    times, samp.units = samp.units, dose.units = dose.units, 
		    dose.name = dose.name, 
		    #toler = toler, 
		    extrap = F, 
		    #dos = calib.fit.out@dos, 
		    rname = rname, xname = xname)
	}
	else {if(calib.fit.out@type=="thpl"){
		calib.out <- calib.thpl.pom(calib.fit.out, y0=y0, conf = conf, 
			m = calib.fit.out@m, dilution = dilution, samp.names, truth, 
			times, samp.units = samp.units, dose.units = dose.units, 
			dose.name = dose.name, 
			#toler = toler, 
			extrap = F,
		    #dos = calib.fit.out@dos, 
		    rname = rname, xname = xname 
		    #maxit=maxit
		    )
	}
	else {if(calib.fit.out@type=="fpl"){
    	calib.out <- calib.fpl.pom(calib.fit.out, y0=y0, conf = conf, 
    		m = calib.fit.out@m, dilution = dilution, samp.names, truth, 
    		times, samp.units = samp.units, dose.units = dose.units, 
    		dose.name = dose.name, 
    		#toler = toler, 
    		extrap = F,
      		#dos = calib.fit.out@dos, 
      		rname = rname, xname = xname 
      		#maxit=maxit
      		)
	}
    }
    }
 calib.out
}

calib.fpl.pom <- function(fpl.out, y0, conf = 0.9, m = fpl.out@m, dilution = 1, 
	samp.names, truth, times, samp.units = "", dose.units = "", 
    dose.name = "", 
    #maxit = 1000, toler = 1e-005, 
    #dos = fpl.out@dos,
    rname="response", extrap=F, xname="x")
{
###############################################################################
## Computes the Calibration Statistics                                       ##
##                                                                           ##
## input variables:                                                          ##
##  fpl.out: output from the fpl or fpl.pom functions                        ##
##  y0: a vector of the MEAN response values for the unknown x               ##
##  conf: confidence level for calibration estimates                         ##
##  m: number of terms in the average for y0                                 ##
##     BE CAREFUL: by default this will read in the same number as           ##
##     in fpl.out - the means for y0 may be based on a different             ##
##     number of reps for each sample type                                   ##
##  dilution: dilution factor                                                ##
##  samp.names: names of the unknowns                                        ##
##  truth: names of the knowns?                                              ##
##  samp.units, dose.units, dose.name: self-explanatory and blank by default ##
##  maxit, toler, dos: used for dos                                          ##
##  rname: this argument was added by MM on 3/27/03.  This is the name of    ##
##         of the response variable, which was get used when printing        ##  
## comments added by Matthew Mitchell on 3/14/03 are indicated with ##       ##
## the original comments have # or ###                                       ##
##                                                                           ##
## added extrap option as in calib.lin.pom to give user option of            ##
## extrapolating out of range values if desired: MM 9/26/03                  ##
###############################################################################

    n <- length(y0)
    y0 <- array(y0, dim = n)
    if(length(dilution) == 1)
        dilution <- rep(dilution, n)
    dcorr <- 1/dilution
    b <- fpl.out@coefficients
    cov.un <- fpl.out@cov.unscaled    ### Min modified this
    min.x <- min(fpl.out@x)
    max.x <- max(fpl.out@x)     ## MM moved the line from the top to here
    max.range <- fpl.model(max.x, b, parm=fpl.out@parm) ## MM added parm argument on 4/10/03
    min.range <- fpl.model(min.x, b, parm=fpl.out@parm)
    if(sign(b[2] - b[1]) < 0) {
#    oor <- ifelse(y0 > b[1] | y0 < b[2], "x", "o")
        oor <- ifelse(y0 > min.range[1] | y0 < max.range[1], "x", "o")
        too.low <- ifelse(y0 > min.range[1], 1, 0)
        too.high <- ifelse(y0 < max.range[1], 1, 0)
    }
    if(sign(b[2] - b[1]) > 0) {
        oor <- ifelse(y0 < min.range[1] | y0 > max.range[1], "x", "o")
        too.low <- ifelse(y0 < min.range[1], 1, 0)
        too.high <- ifelse(y0 > max.range[1], 1, 0)
    }
## MM added these next 2 lines on 3/27/03 because m may be different here
## than in fpl.out and its these values of m that will determine the mdc
## This is added an output variable 

    cmdc <- mdc.calc(fpl.out, m=m, type="fpl")
    cmdc <- cmdc$mdc*dcorr
    
    cal.val <- fpl.inverse(fpl.out, y0)
    cal.val <- as.array(cal.val)    
    ### correct calibrated value for dilution  ##next line does nothing??
    # cal.val <- cal.val 04-19-07: Removing this line. DVS
    tcrit <- qt((conf + 1)/2, fpl.out@df.residual)
 ## if dos is False
## 02-08-07: This call to C does not seem relavent. DVS
#    if(!dos) {
#        xl <- rep(0, n)
#        xu <- rep(max.x, n)
#        svec <- c(cov.un[, 1], cov.un[, 2], cov.un[, 3], cov.un[, 4])
#        if(!is.loaded(symbol.C("fplpomcal")))
#            dyn.load("/local/lib/calib/sun4/src/fplpom.o")
#        clims.lpred.out <- .C("fplpomcal",
#            xl = as.double(xl),
#            xu = as.double(xu),
#            as.double(b),
#            as.double(svec),
#            as.integer(maxit),
#            as.double(toler),
#            as.double(fpl.out@theta),
#            as.double(tcrit),
#            as.double(fpl.out@sigma),
#            as.double(y0),
#            as.integer(sign(b[2] - b[1])),
#            as.integer(n),
#            as.integer(fpl.out@m),
#            as.integer(fpl.out@parm),
#            as.double(-1))
#        xl <- rep(0, n)
#        xu <- rep(max.x, n)
#        clims.upred.out <- .C("fplpomcal",
#            xl = as.double(xl),
#            xu = as.double(xu),
#            as.double(b),
#            as.double(svec),
#            as.integer(maxit),
#            as.double(toler),
#            as.double(fpl.out@theta),
#            as.double(tcrit),
#            as.double(fpl.out@sigma),
#            as.double(y0),
#            as.integer(sign(b[2] - b[1])),
#            as.integer(n),
#            as.integer(fpl.out@m),
#            as.integer(fpl.out@parm),
#            as.double(1))
#        if(sign(b[2] - b[1]) < 0) {
#            xup.inver <- ifelse(clims.upred.out$xl >= max.x,
#                NA, clims.upred.out$xl)
#            xlow.inver <- ifelse(clims.lpred.out$xl >= 
#                max.x, NA, clims.lpred.out$xl)
#            if(!extrap){                               ## MM 9/26
#            cal.val <- ifelse(y0 > b[1], 0, cal.val)
#            cal.val <- ifelse(y0 < b[2], NA, cal.val)}
#        }
#        else {
#            xup.inver <- ifelse(clims.lpred.out$xl >= max.x,
#                NA, clims.lpred.out$xl)
#            xlow.inver <- ifelse(clims.upred.out$xl >= 
#                max.x, NA, clims.upred.out$xl)
#            if(!extrap){                             ## MM 9/26
#            cal.val <- ifelse(y0 < b[1], 0, cal.val)
#            cal.val <- ifelse(y0 > b[2], NA, cal.val)}
#        }
#    }
 ## is dos=T (default)
 ######################################
 ## computes the inverse limits here ##
 ######################################
## 02-08-07: This call to C does not seem relavent. DVS 
#    else {
        seqy <- seq(as.vector(fpl.model(min(fpl.out@x), fpl.out@
            coefficients, parm = fpl.out@parm)), as.vector(
            fpl.model(max(fpl.out@x), fpl.out@coefficients, parm = 
            fpl.out@parm)), len = 6)
        boundx <- c(fpl.inverse(fpl.out, seqy[1:5]), max.x)[-1]
        xp <- c(seq(0, boundx[1], len = 170), seq(boundx[1], 
            boundx[2], len = 50), seq(boundx[2], boundx[3], 
            len = 30), seq(boundx[3], boundx[4], len = 50), 
            seq(boundx[4], boundx[5], len = 170))
        ypred <- fpl.model(xp, b, parm = fpl.out@parm)
        qn1 <- ((as.vector(ypred))^(2 * fpl.out@theta))/m
        qn2 <- diag(attributes(ypred)$gradient %*% cov.un %*% t(
            attributes(ypred)$gradient))
        sigma <- fpl.out@sigma
        ypl <- as.vector(ypred) + sign(b[2] - b[1]) * tcrit * 
            sigma * sqrt(qn1 + qn2)
        ypu <- as.vector(ypred) - sign(b[2] - b[1]) * tcrit * 
            sigma * sqrt(qn1 + qn2)
        xl <- approx(ypl, xp, y0,ties="ordered")$y
        xu <- approx(ypu, xp, y0,ties="ordered")$y
        if(sign(b[2] - b[1]) < 0) {
            if(!extrap){                                 ## MM 9/26/03
            cal.val <- ifelse(y0 > b[1], 0, cal.val)
            cal.val <- ifelse(y0 < b[2], NA, cal.val)
            cal.val <- ifelse(cal.val > max.x, NA, cal.val)}
            xup <- ifelse(is.na(xu) | xu >= max.x, NA, xu)
            xup <- ifelse(is.na(xup) & cal.val == 0, 0, xup)
            xup.inver <- ifelse(y0 < b[2], NA, xup)
            xlow <- ifelse(y0 > b[1], 0, xl)
            xlow.inver <- ifelse(is.na(xlow) & is.finite(
                cal.val), 0, xlow)
          }
        else {
            if(!extrap){                               ## MM 9/26/03
            cal.val <- ifelse(y0 < b[1], 0, cal.val)
            cal.val <- ifelse(y0 > b[2], NA, cal.val)
            cal.val <- ifelse(cal.val > max.x, NA, cal.val)
            }
            xup <- ifelse(xu >= max.x, NA, xu)
            xup.inver <- ifelse(is.na(xup) & cal.val == 0, 0, xup)
            xlow <- ifelse(y0 < b[1], 0, xl)    
    ### xlow.inver <-  ifelse(is.na(xlow) & y0<b[2],0,xlow)
            xlow.inver <- ifelse(is.na(xlow) & is.finite(cal.val), 0, xlow)
        }
#    }
  ##############################################
  ## computes the Wald intervals              ##
  ## which are computed with the delta method ##
  ##############################################
    if(fpl.out@parm == 1) {
        if(b[1] > b[2]){
        	hyb <- rep(NA,length(y0))
        	hyb[b[1] < y0] <- 0
        	hyb[b[1] >= y0] <- b[3] * exp((1/b[4]) *
                log((b[1] - y0[b[1] >= y0])/(y0[b[1] >= y0] - b[2])))
            }
        ## 05-14-07: Removed the following two lines since errors were
        ## being produced in the ifelse statement
        #hyb <- ifelse(b[1] < y0, 0, b[3] * exp((1/b[4]) *
        #   log((b[1] - y0)/(y0 - b[2]))))
        else {
        	hyb <- rep(NA,length(y0))
        	hyb[b[1] > y0] <- 0
        	hyb[b[1] <= y0] <- b[3] * exp((1/b[4]) * 
                log((b[1] - y0[b[1] <= y0])/(y0[b[1] <= y0] - b[2])))
        ## 05-14-07: Removed the following two lines since errors were
        ## being produced in the ifelse statement
        #hyb <- ifelse(b[1] > y0, 0, b[3] * exp((1/b[4]) * 
        #        log((b[1] - y0)/(y0 - b[2]))))
	        dh.dy <- hyb * ((b[2] - b[1])/(b[4] * (b[1] - y0) * (y0 -  b[2])))
    	    dh.db1 <- hyb/(b[4] * (b[1] - y0))
        	dh.db2 <- hyb/(b[4] * (y0 - b[2]))
	        dh.db3 <- hyb/b[3]
			dh.db4 <- rep(NA,length(y0))
			dh.db4[hyb==0] <- 0
			dh.db4[hyb!=0] <- ( - hyb[hyb!=0]/(b[4] * b[4])) * 
				log((b[1] - y0[hyb!=0])/(y0[hyb!=0] - b[2]))
			}
        ## 05-14-07: Removed the following two lines since errors were
        ## being produced in the ifelse statement
	    #    dh.db4 <- ifelse(hyb == 0, 0, ( - hyb/(b[4] * b[4])) * 
    	#        log((b[1] - y0)/(y0 - b[2])))
    }
    if(fpl.out@parm == 2) {
        if(b[1] > b[2]){ 
	        hyb <- rep(NA,length(y0))
    	    hyb[b[1] < y0] <- 0
        	hyb[b[1] >= y0] <- exp((1/b[4]) * 
        		log((b[1] - y0[b[1] >= y0])/
        		(y0[b[1] >= y0] - b[2])) + b[3])
        	}
        ## 05-14-07: Removed the following two lines since errors were
        ## being produced in the ifelse statement
            #hyb <- ifelse(b[1] < y0, 0, exp((1/b[4]) * 
            #	log((b[1] - y0)/(y0 - b[2])) + b[3]))
        else{ 
	        hyb <- rep(NA,length(y0))
    	    hyb[b[1] > y0] <- 0
        	hyb[b[1] <= y0] <- exp((1/b[4]) * 
        		log((b[1] - y0[b[1] <= y0])/
        		(y0[b[1] <= y0] - b[2])) + b[3])
        	}
        ## 05-14-07: Removed the following two lines since errors were
        ## being produced in the ifelse statement
        #hyb <- ifelse(b[1] > y0, 0, exp((1/b[4]) * 
        #		log((b[1] - y0)/(y0 - b[2])) + b[3]))
        dh.dy <- (hyb * (b[2] - b[1]))/(b[4] * (b[1] - y0) * (y0 - b[2]))
        dh.db1 <- hyb/(b[4] * (b[1] - y0))
        dh.db2 <- hyb/(b[4] * (y0 - b[2]))
        dh.db3 <- hyb
        dh.db4 <- rep(NA,length(y0))
        dh.db4[hyb == 0] <- 0
        dh.db4[hyb != 0] <- ( - hyb[hyb != 0]/(b[4] * b[4])) * 
            log((b[1] - y0[hyb != 0])/(y0[hyb != 0] - b[2]))
		## 05-14-07: Removed the following two lines since errors were
        ## being produced in the ifelse statement        
        #ifelse(hyb == 0, 0, ( - hyb/(b[4] * b[4])) * 
        #    log((b[1] - y0)/(y0 - b[2])))
    }
    xnot.hat <- hyb
    sigma2 <- fpl.out@sigma * fpl.out@sigma
    var.xnot.hat <- ((dh.dy * dh.dy) * sigma2 * (y0^(2 * fpl.out@theta)))/m +
        sigma2 * (dh.db1 * (dh.db1 * cov.un[1, 1] + 
        dh.db2 * cov.un[2, 1] + dh.db3 * cov.un[3, 1] + dh.db4 * 
        cov.un[4, 1]) + dh.db2 * (dh.db1 * cov.un[1, 2] + 
        dh.db2 * cov.un[2, 2] + dh.db3 * cov.un[3, 2] + dh.db4 * 
        cov.un[4, 2]) + dh.db3 * (dh.db1 * cov.un[1, 3] + 
        dh.db2 * cov.un[2, 3] + dh.db3 * cov.un[3, 3] + dh.db4 * 
        cov.un[4, 3]) + dh.db4 * (dh.db1 * cov.un[1, 4] + 
        dh.db2 * cov.un[2, 4] + dh.db3 * cov.un[3, 4] + dh.db4 * 
        cov.un[4, 4]))
    tcrit <- qt((1 + conf)/2, fpl.out@df.residual)
    se.xhat <- sqrt(var.xnot.hat)
    xup <- xnot.hat + tcrit * se.xhat
    xlow <- xnot.hat - tcrit * se.xhat
    if(!extrap)                              ## MM 9/26
		cal.val <- ifelse(hyb < max.x, hyb, NA)
     ## it seems like this was computed already
      ## then recomputed as hyp and then set equal?
    
    if(extrap)
		cal.val <- hyb
    ## will not have NAs for high values but have 0s for ones too low
    xlow.wald <- ifelse(is.finite(xlow), ifelse(xlow < 0, 0, xlow), NA)
    xup.wald <- ifelse(xup <= max.x, xup, NA)
    xup.waldl <- NULL
    xlo.waldl <- NULL
  ####################################################
  ## Another set of confidence intervals for parm=1 ##
  ####################################################
    if(fpl.out@parm == 1) {
		if(b[1] > b[2]){ 
	        hybl <- rep(NA,length(y0))
    	    hybl[b[1] < y0] <- 0
        	hybl[b[1] >= y0] <- log(b[3]) + (1/b[4]) * 
        		log((b[1] - y0[b[1] >= y0])/(y0[b[1] >= y0] - b[2]))
        	}
        else{ 
	        hybl <- rep(NA,length(y0))
    	    hybl[b[1] > y0] <- 0
        	hybl[b[1] <= y0] <- log(b[3]) + (1/b[4]) * 
        		log((b[1] - y0[b[1] <= y0])/(y0[b[1] <= y0] - b[2]))
        	}		
        ## 05-14-07: Was getting errors from the next line, added above
        ## to account for shape of curve and avoid taking the log of
        ## negative values. DVS
        #hybl <- log(b[3]) + (1/b[4]) * log((b[1] - y0)/(y0 - b[2]))
        dh.dyl <-  - (1/(b[4] * (b[1] - y0)) + 1/(b[4] * (y0 - b[2])))
        dh.db1l <- 1/(b[4] * (b[1] - y0))
        dh.db2l <- 1/(b[4] * (y0 - b[2]))
        dh.db3l <- 1/b[3]
        if(b[1] > b[2]){ 
	        dh.db4l <- rep(NA,length(y0))
    	    dh.db4l[b[1] < y0] <- 0
        	dh.db4l[b[1] >= y0] <- (1/(b[4] * b[4])) * 
        		log((y0[b[1] >= y0] - b[2])/(b[1] - y0[b[1] >= y0]))
        	}
        else{ 
	        dh.db4l <- rep(NA,length(y0))
    	    dh.db4l[b[1] > y0] <- 0
        	dh.db4l[b[1] <= y0] <- (1/(b[4] * b[4])) * 
        		log((y0[b[1] <= y0] - b[2])/(b[1] - y0[b[1] <= y0]))
        	}		
#        dh.db4l <- (1/(b[4] * b[4])) * log((y0 - b[2])/(b[1] - y0))
        var.hybl <- (dh.dyl * dh.dyl) * sigma2 * (y0^(2 * 
            fpl.out@theta)) + sigma2 * (dh.db1l * (dh.db1l * 
            cov.un[1, 1] + dh.db2l * cov.un[2, 1] + dh.db3l *
            cov.un[3, 1] + dh.db4l * cov.un[4, 1]) + 
            dh.db2l * (dh.db1l * cov.un[1, 2] + dh.db2l * 
            cov.un[2, 2] + dh.db3l * cov.un[3, 2] + dh.db4l *
            cov.un[4, 2]) + dh.db3l * (dh.db1l * cov.un[1, 
            3] + dh.db2l * cov.un[2, 3] + dh.db3l * cov.un[
            3, 3] + dh.db4l * cov.un[4, 3]) + dh.db4l * (
            dh.db1l * cov.un[1, 4] + dh.db2l * cov.un[2, 4] +
            dh.db3l * cov.un[3, 4] + dh.db4l * cov.un[4, 4]))
        xup.waldl <- hybl + tcrit * sqrt(var.hybl)
        xlo.waldl <- hybl - tcrit * sqrt(var.hybl)
        log.up <- log(max.x)
        xup.waldl <- ifelse(xup.waldl > log.up, max.x, exp(xup.waldl))
        xlo.waldl <- ifelse(xlo.waldl > log.up, NA, exp(xlo.waldl))
        xup.waldl <- ifelse(is.na(cal.val), NA, xup.waldl)
        xlo.waldl <- ifelse(is.na(cal.val), NA, xlo.waldl)
    }
    xlow.inver <- dcorr * xlow.inver
    xup.inver <- dcorr * xup.inver
    xlow.wald <- dcorr * xlow.wald
    xup.wald <- dcorr * xup.wald
    cal.val <- dcorr * cal.val
    se.xhat <- dcorr * se.xhat

    ## next line added by MM on 3/27/03 to compute when the
    ## estimated calibration is less than the mdc
    lmdc <- ifelse(as.vector(cal.val) < as.vector(cmdc), "X", "O")

 ############################################################
 ## creates a data set with all the calibration statistics ##
 ############################################################

    if(any(is.null(xup.waldl))) { ## 05-14-07: Note, made a change from
    ## is.na to is.null. DVS
        clims.df <- new("Calib",Estimated.x = as.vector(cal.val), 
            PredStdErr = as.vector(se.xhat), inver.low = 
            as.vector(xlow.inver), inver.up = as.vector(
            xup.inver), wald.low = as.vector(xlow.wald), 
            wald.up = as.vector(xup.wald), avg.response = 
            as.vector(y0), dilution = as.vector(dilution), 
            oor = as.vector(oor), mdc = as.vector(cmdc),
            lmdc = as.vector(lmdc))
    }
    else {
        clims.df <- new("Calib",Estimated.x = as.vector(cal.val), 
            PredStdErr = as.vector(se.xhat), inver.low = 
            as.vector(xlow.inver), inver.up = as.vector(
            xup.inver), wald.low = as.vector(xlow.wald), 
            wald.up = as.vector(xup.wald), avg.response = 
            as.vector(y0), waldl.low = as.vector(xlo.waldl),
            waldl.up = as.vector(xup.waldl), dilution = 
            as.vector(dilution), oor = as.vector(oor), mdc = as.vector(cmdc),
            lmdc = as.vector(lmdc))
    }
    ## Don't think that these are necessary anymore
    #names(clims.df)[1:2] <- c(paste("Estimated",xname,sep=""),"PredStdErr")  ## MWM 8/31/04 added
    #names(clims.df)[match("response",names(clims.df))] <- paste("avg",rname,sep=".") ## MWM 9/2/04 added

clims.df@wald.up <- ifelse(is.na(clims.df@wald.up) & !is.na(clims.df@wald.low),
                           clims.df@Estimated.x + (clims.df@Estimated.x-clims.df@wald.low),
                           clims.df@wald.up)   ## MWM 9/2/04
    
    
    lab.given <- F
 #########################################################
 ## used if set value for truth, known controls         ##
 #########################################################
    if(!missing(truth)) {
        ord.truth <- order(truth)
        y0 <- y0[ord.truth]
        cal.val <- cal.val[ord.truth]
        xlow.inver <- xlow.inver[ord.truth]
        xup.inver <- xup.inver[ord.truth]
        xlow.wald <- xlow.wald[ord.truth]
        xup.wald <- xup.wald[ord.truth]
        truth <- truth[ord.truth]
        labs <- list(name = "Control", plot.ind = 1:length(
            truth),    # ## too.low=too.low, 
###    too.high=too.high, 
        too.low = too.low[ord.truth], too.high = too.high[
            ord.truth],     
    ### i think its here that the truth (| samp names ?) get set as 
### xlabs which turn into tick labels on axis in plot.calib
### order truth is screwing things up for plot
### subscripting too.low and too.high as above fixed it !! (moc 5/8/95)
        xlabs = as.character(truth), cont.pl = F, xunits = 
            samp.units, dose.units = dose.units, dose.name
             = dose.name)
        lab.given <- T
        clims.df@truth <- truth
        if(missing(samp.names))
            clims.df@samp.names <- truth
        else clims.df@samp.names <- samp.names
        ## This line does not make sense in the context of
        ## rewriting things to S4 format

## Note sure what to do with this line yet
        #clims.df <- clims.df[ord.truth,  ]
    }
    
    ## This line is not necessary as the slot will simply be empty
#    else attributes(clims.df)$truth <- NULL
    
##########################################################
## used if times are provided                           ##
##########################################################
    if(!missing(times)) {
      times = as.character(times)  ## Added by DVS 10/9/06
        clims.df@times <- times
        if(missing(samp.names))
            clims.df@samp.names <- times
        else clims.df@samp.names <- samp.names
        labs <- list(name = "Time", plot.ind = times, too.low
             = too.low, too.high = too.high, xlabs = 
            as.character(times), cont.pl = T, xunits = 
            samp.units, dose.units = dose.units, dose.name
             = dose.name)
        lab.given <- T
    }

    ## Not necessary, slot will simply be empty
#    else attributes(clims.df)$times <- NULL
    
##########################################################
## will be executed is truth and times are missing      ##
##########################################################
    if(!lab.given) {
        labs <- list(name = "Sample", plot.ind = 1:length(y0), 
            too.low = too.low, too.high = too.high, xlabs
             = as.character(1:length(y0)), cont.pl = F, 
            xunits = samp.units, dose.units = dose.units, 
            dose.name = dose.name)
        if(!missing(samp.names)) {
            clims.df@samp.names <- samp.names
            labs$xlabs <- samp.names
            if(is.numeric(samp.names))
                labs$plot.ind <- samp.names
        }
       # else attributes(clims.df)$samp.names <- NULL
    }
    clims.df@labels <- labs
    clims.df@max.x <- max.x

## MM 9/29/03
    clims.df@extrap <- extrap
     

## MM added these 2 lines on 3/27/03 to indicate whether m here is the same
## as in the fpl print out (this affects the MDC)
   ifelse(all(m == fpl.out@m), {repeq <- T}, {repeq <- F})
   clims.df@repeq <- repeq
    
## MM added this line on 3/27/03 that the actual name of the response
## variable will be shown on the printout
    clims.df@rname <- rname
  
## MM 3/27/03 see change earlier
    if(dcorr[1] == 1){
    	clims.df@mdc <- fpl.out@mdc
    }
    clims.df@conf.level <- conf
    clims.df@calib.fit <- fpl.out
    clims.df
}

calib.thpl.pom <- function(thpl.out, y0, conf = 0.9, m = thpl.out@m, dilution = 1, 
    samp.names, truth, times, samp.units = "", dose.units = "", 
    dose.name = "", 
    #maxit = 1000, toler = 1e-005, dos = thpl.out@dos,
    rname="response",extrap=F,xname="x")
{
###############################################################################
## Computes the Calibration Statistics                                       ##
##                                                                           ##
## input variables:                                                          ##
##  thpl.out: output from the thpl or thpl.pom functions                     ##
##  y0: a vector of the MEAN response values for the unknown x               ##
##  conf: confidence level for calibration estimates                         ##
##  m: number of terms in the average for y0                                 ##
##     BE CAREFUL: by default this will read in the same number as           ##
##     in thpl.out - the means for y0 may be based on a different            ##
##     number of reps for each sample type                                   ##
##  dilution: dilution factor                                                ##
##  samp.names: names of the unknowns                                        ##
##  truth: names of the knowns?                                              ##
##  samp.units, dose.units, dose.name: self-explanatory and blank by default ##
##  maxit, toler, dos: used for dos                                          ##
##  rname: this argument was added by MM on 3/27/03.  This is the name of    ##
##         of the response variable, which was get used when printing        ##  
## written  by Matthew Mitchell on 5/20/03                                   ##
## the original comments have # or ###                                       ##
## add extrapolation choice on 8/30/04                                       ##
###############################################################################
  
    n <- length(y0)
    y0 <- array(y0, dim = n)
    if(length(dilution) == 1)
        dilution <- rep(dilution, n)
    dcorr <- 1/dilution
    b <- thpl.out@coefficients
    cov.un <- thpl.out@cov.unscaled    ### Min modified this
    min.x <- min(thpl.out@x)
    max.x <- max(thpl.out@x)     ## MM moved the line from the top to here
    max.range <- thpl.model(max.x, b) 
    min.range <- thpl.model(min.x, b)
    if(sign(b[2] - b[1]) < 0) {
#    oor <- ifelse(y0 > b[1] | y0 < b[2], "x", "o")
        oor <- ifelse(y0 > min.range[1] | y0 < max.range[1], "x", "o")
        too.low <- ifelse(y0 > min.range[1], 1, 0)
        too.high <- ifelse(y0 < max.range[1], 1, 0)
    }
    if(sign(b[2] - b[1] > 0)) {
#    oor <- ifelse(y0 < b[1] | y0 > b[2], "x", "o")
        oor <- ifelse(y0 < min.range[1] | y0 > max.range[1], "x", "o")
        too.low <- ifelse(y0 < min.range[1], 1, 0)
        too.high <- ifelse(y0 > max.range[1], 1, 0)
    }

## MM added these next 2 lines on 3/27/03 because m may be different here
## than in thpl.out and its these values of m that will determine the mdc
## This is added an output variable 

    cmdc <-  mdc.calc(thpl.out, m=m,type="thpl")
    cmdc <-  cmdc$mdc*dcorr
    
    cal.val <- thpl.inverse(thpl.out, y0)
    cal.val <- as.array(cal.val)    
    ### correct calibrated value for dilution
    cal.val <- cal.val
    tcrit <- qt((conf + 1)/2, thpl.out@df.residual)
 ## if dos is False
## 02-08-07: This call to C does not seem relavent. DVS
#    if(!dos) {
#        xl <- rep(0, n)
#        xu <- rep(max.x, n)
#        svec <- c(cov.un[, 1], cov.un[, 2], cov.un[, 3])
#        if(!is.loaded(symbol.C("thplpomcal")))
#            dyn.load("/local/lib/calib/sun4/src/thplpom.o")
#        clims.lpred.out <- .C("thplpomcal",
#            xl = as.double(xl),
#            xu = as.double(xu),
#            as.double(b),
#            as.double(svec),
#            as.integer(maxit),
#            as.double(toler),
#            as.double(thpl.out@theta),
#            as.double(tcrit),
#            as.double(thpl.out@sigma),
#            as.double(y0),
#            as.integer(sign(b[2] - b[1])),
#            as.integer(n),
#            as.integer(thpl.out@m),
#            as.double(-1))
#        xl <- rep(0, n)
#        xu <- rep(max.x, n)
#        clims.upred.out <- .C("thplpomcal",
#            xl = as.double(xl),
#            xu = as.double(xu),
#            as.double(b),
#            as.double(svec),
#            as.integer(maxit),
#            as.double(toler),
#            as.double(thpl.out@theta),
#            as.double(tcrit),
#            as.double(thpl.out@sigma),
#            as.double(y0),
#            as.integer(sign(b[2] - b[1])),
#            as.integer(n),
#            as.integer(thpl.out@m),
#            as.double(1))
#        if(sign(b[2] - b[1]) < 0) {
#            xup.inver <- ifelse(clims.upred.out$xl >= max.x,
#                NA, clims.upred.out$xl)
#            xlow.inver <- ifelse(clims.lpred.out$xl >= 
#                max.x, NA, clims.lpred.out$xl)
#            if(!extrap){                                   ## MWM 8/30/04
#            cal.val <- ifelse(y0 > b[1], 0, cal.val) 
#            cal.val <- ifelse(y0 < b[2], NA, cal.val)}
#        }
#        else {
#            xup.inver <- ifelse(clims.lpred.out$xl >= max.x,
#                NA, clims.lpred.out$xl)
#            xlow.inver <- ifelse(clims.upred.out$xl >= 
#                max.x, NA, clims.upred.out$xl)
#            if(!extrap){                                 ## MWM 8/30/04
#            cal.val <- ifelse(y0 < b[1], 0, cal.val)
#            cal.val <- ifelse(y0 > b[2], NA, cal.val)}
#        }
#    }
 ## is dos=T (default)
 ######################################
 ## computes the inverse limits here ##
 ######################################
#    else {
        seqy <- seq(as.vector(thpl.model(min(thpl.out@x), thpl.out@coefficients)),
                as.vector(thpl.model(max(thpl.out@x), thpl.out@coefficients)), len = 6)
        boundx <- c(thpl.inverse(thpl.out, seqy[1:5]), max.x)[-1]
        xp <- c(seq(0, boundx[1], len = 170), seq(boundx[1], 
            boundx[2], len = 50), seq(boundx[2], boundx[3], 
            len = 30), seq(boundx[3], boundx[4], len = 50), 
            seq(boundx[4], boundx[5], len = 170))
        ypred <- thpl.model(xp, b)
        qn1 <- ((as.vector(ypred))^(2 * thpl.out@theta))/m
        qn2 <- diag(attributes(ypred)$gradient %*% cov.un %*% t(
            attributes(ypred)$gradient))
        sigma <- thpl.out@sigma
        ypl <- as.vector(ypred) + sign(b[2] - b[1]) * tcrit * 
            sigma * sqrt(qn1 + qn2)
        ypu <- as.vector(ypred) - sign(b[2] - b[1]) * tcrit * 
            sigma * sqrt(qn1 + qn2)
        xl <- approx(ypl, xp, y0, ties = "ordered")$y
        xu <- approx(ypu, xp, y0, ties = "ordered")$y
        if(sign(b[2] - b[1]) < 0) {
            if(!extrap){                                   ## MWM 8/30/04
            cal.val <- ifelse(y0 > b[1], 0, cal.val)
            cal.val <- ifelse(y0 < b[2], NA, cal.val)
            cal.val <- ifelse(cal.val > max.x, NA, cal.val)}
            xup <- ifelse(is.na(xu) | xu >= max.x, NA, xu)
            xup <- ifelse(is.na(xup) & cal.val == 0, 0, xup)
            xup.inver <- ifelse(y0 < b[2], NA, xup)
            xlow <- ifelse(y0 > b[1], 0, xl)
            xlow.inver <- ifelse(is.na(xlow) & is.finite(
                cal.val), 0, xlow)
        }
        else {
            if(!extrap){                                   ## MWM 8/30/04
            cal.val <- ifelse(y0 < b[1], 0, cal.val)
            cal.val <- ifelse(y0 > b[2], NA, cal.val)
            cal.val <- ifelse(cal.val > max.x, NA, cal.val)}
            xup <- ifelse(xu >= max.x, NA, xu)
            xup.inver <- ifelse(is.na(xup) & cal.val == 0, 0, xup)
            xlow <- ifelse(y0 < b[1], 0, xl)    
    ### xlow.inver <-  ifelse(is.na(xlow) & y0<b[2],0,xlow)
            xlow.inver <- ifelse(is.na(xlow) & is.finite(cal.val), 0, xlow)}
        
#    }
  ##############################################
  ## computes the Wald intervals              ##
  ## which are computed with the delta method ##
  ##############################################
   
        if(b[1] > b[2])
            hyb <- ifelse(b[1] < y0, 0, b[3]*((b[1] - y0)/(y0 - b[2])))
        else hyb <- ifelse(b[1] > y0,0, b[3]*((b[1] - y0)/(y0 - b[2]))) 
        dh.dy <- hyb * ((b[2] - b[1])/( (b[1] - y0) * (y0 -  b[2])))
        dh.db1 <- hyb/(b[1] - y0)
        dh.db2 <- hyb/(y0 - b[2])
        dh.db3 <- hyb/b[3]
           
   
    xnot.hat <- hyb
    sigma2 <- thpl.out@sigma * thpl.out@sigma
    var.xnot.hat <- ((dh.dy * dh.dy) * sigma2 * (y0^(2 * thpl.out@
        theta)))/m + sigma2 * (dh.db1 * (dh.db1 * cov.un[1, 1] + 
        dh.db2 * cov.un[2, 1] + dh.db3 * cov.un[3, 1]) +
        dh.db2 * (dh.db1 * cov.un[1, 2] + 
        dh.db2 * cov.un[2, 2] + dh.db3 * cov.un[3, 2]) +
        dh.db3 * (dh.db1 * cov.un[1, 3] + 
        dh.db2 * cov.un[2, 3] + dh.db3 * cov.un[3, 3]))
      
    tcrit <- qt((1 + conf)/2, thpl.out@df.residual)
    se.xhat <- sqrt(var.xnot.hat)
    xup <- xnot.hat + tcrit * se.xhat
    xlow <- xnot.hat - tcrit * se.xhat
    if(!extrap){                                 ## MWM 9/30/04
    cal.val <- ifelse(hyb < max.x, hyb, NA)}
    if(extrap){                                ## MWM 9/30/04
    cal.val <- hyb}
    xlow.wald <- ifelse(is.finite(xlow), ifelse(xlow < 0, 0, xlow), NA)
    xup.wald <- ifelse(xup <= max.x, xup, NA)
    xup.waldl <- NULL
    xlo.waldl <- NULL
  ####################################################
  ## Another set of confidence intervals            ##
  ####################################################
  
        hybl <- log(b[3]) + log((b[1] - y0)/(y0 - b[2]))
        dh.dyl <-  - (b[1] - y0) + 1/(y0 - b[2])
        dh.db1l <- 1/ (b[1] - y0)
        dh.db2l <- 1/ (y0 - b[2])
        dh.db3l <- 1/b[3]
      
        var.hybl <- (dh.dyl * dh.dyl) * sigma2 * (y0^(2 * 
            thpl.out@theta)) + sigma2 * (dh.db1l * (dh.db1l * 
            cov.un[1, 1] + dh.db2l * cov.un[2, 1] + dh.db3l *
            cov.un[3, 1]) + 
            dh.db2l * (dh.db1l * cov.un[1, 2] + dh.db2l * 
            cov.un[2, 2] + dh.db3l * cov.un[3, 2])
            + dh.db3l * (dh.db1l * cov.un[1, 
            3] + dh.db2l * cov.un[2, 3] + dh.db3l * cov.un[
            3, 3]))
        xup.waldl <- hybl + tcrit * sqrt(var.hybl)
        xlo.waldl <- hybl - tcrit * sqrt(var.hybl)
        log.up <- log(max.x)
        xup.waldl <- ifelse(xup.waldl > log.up, max.x, exp(xup.waldl))
        xlo.waldl <- ifelse(xlo.waldl > log.up, NA, exp(xlo.waldl))
        xup.waldl <- ifelse(is.na(cal.val), NA, xup.waldl)
        xlo.waldl <- ifelse(is.na(cal.val), NA, xlo.waldl)
   
    xlow.inver <- dcorr * xlow.inver
    xup.inver <- dcorr * xup.inver
    xlow.wald <- dcorr * xlow.wald
    xup.wald <- dcorr * xup.wald
    cal.val <- dcorr * cal.val
    se.xhat <- dcorr * se.xhat

    ## next line added by MM on 3/27/03 to compute when the
    ## estimated calibration is less than the mdc
    lmdc <- ifelse(as.vector(cal.val) < as.vector(cmdc), "X", "O")


 ############################################################
 ## creates a data set with all the calibration statistics ##
 ############################################################
    if(is.null(xup.waldl)) {
        clims.df <- new("Calib", Estimated.x = as.vector(cal.val), 
            PredStdErr = as.vector(se.xhat), inver.low = 
            as.vector(xlow.inver), inver.up = as.vector(
            xup.inver), wald.low = as.vector(xlow.wald), 
            wald.up = as.vector(xup.wald), avg.response = 
            as.vector(y0), dilution = as.vector(dilution), 
            oor = as.vector(oor), mdc = as.vector(cmdc),
            lmdc = as.vector(lmdc))
    }
    else {
        clims.df <- new("Calib", Estimated.x  = as.vector(cal.val), 
            PredStdErr = as.vector(se.xhat), inver.low = 
            as.vector(xlow.inver), inver.up = as.vector(
            xup.inver), wald.low = as.vector(xlow.wald), 
            wald.up = as.vector(xup.wald), avg.response = 
            as.vector(y0), waldl.low = as.vector(xlo.waldl),
            waldl.up = as.vector(xup.waldl), dilution = 
            as.vector(dilution), oor = as.vector(oor), mdc = as.vector(cmdc),
            lmdc = as.vector(lmdc))
    }
#    names(clims.df)[1:2] <- c(paste("Estimated",xname),"PredStdErr")  ## MWM 8/31/04 added
#     names(clims.df)[match("response",names(clims.df))] <- paste("avg",rname,sep=".") ## MWM 9/2/04
   clims.df@wald.up <- ifelse(is.na(clims.df@wald.up) & !is.na(clims.df@wald.low),
                           clims.df@Estimated.x + (clims.df@Estimated.x-clims.df@wald.low),
                           clims.df@wald.up)   ## MWM 9/2/04
    
    
    lab.given <- F
 #########################################################
 ## used if set value for truth, known controls         ##
 #########################################################
    if(!missing(truth)) {
        ord.truth <- order(truth)
        y0 <- y0[ord.truth]
        cal.val <- cal.val[ord.truth]
        xlow.inver <- xlow.inver[ord.truth]
        xup.inver <- xup.inver[ord.truth]
        xlow.wald <- xlow.wald[ord.truth]
        xup.wald <- xup.wald[ord.truth]
        truth <- truth[ord.truth]
        labs <- list(name = "Control", plot.ind = 1:length(
            truth),    # ## too.low=too.low, 
###    too.high=too.high, 
        too.low = too.low[ord.truth], too.high = too.high[
            ord.truth],     
    ### i think its here that the truth (| samp names ?) get set as 
### xlabs which turn into tick labels on axis in plot.calib
### order truth is screwing things up for plot
### subscripting too.low and too.high as above fixed it !! (moc 5/8/95)
        xlabs = as.character(truth), cont.pl = F, xunits = 
            samp.units, dose.units = dose.units, dose.name
             = dose.name)
        lab.given <- T
        clims.df@truth <- truth
        if(missing(samp.names))
            clims.df@samp.names <- truth
        else clims.df@samp.names <- samp.names
        clims.df <- clims.df[ord.truth,  ]
    }
    else clims.df@truth <- NULL
##########################################################
## used if times are provided                           ##
##########################################################
    if(!missing(times)) {
        clims.df@times <- times
        if(missing(samp.names))
            clims.df@samp.names <- times
        else clims.df@samp.names <- samp.names
        labs <- list(name = "Time", plot.ind = times, too.low
             = too.low, too.high = too.high, xlabs = 
            as.character(times), cont.pl = T, xunits = 
            samp.units, dose.units = dose.units, dose.name
             = dose.name)
        lab.given <- T
    }
    else clims.df@times <- NULL
##########################################################
## will be executed is truth and times are missing      ##
##########################################################

    if(!lab.given) {
        labs <- list(name = "Sample", plot.ind = 1:length(y0), 
            too.low = too.low, too.high = too.high, xlabs
             = as.character(1:length(y0)), cont.pl = F, 
            xunits = samp.units, dose.units = dose.units, 
            dose.name = dose.name)
        if(!missing(samp.names)) {
            clims.df@samp.names <- samp.names
            labs$xlabs <- samp.names
            if(is.numeric(samp.names))
                labs$plot.ind <- samp.names
        }
        else clims.df@samp.names <- NULL
    }
    clims.df@labels <- labs
    clims.df@max.x <- max.x

## MWM 9/30/04
    clims.df@extrap <- extrap


## MM added these 2 lines on 3/27/03 to indicate whether m here is the same
## as in the thpl print out (this affects the MDC)
   ifelse(all(m == thpl.out@m), {repeq <- T}, {repeq <- F})
   clims.df@repeq <- repeq
    
## MM added this line on 3/27/03 that the actual name of the response
## variable will be shown on the printout
    clims.df@rname <- rname
  
## MM 3/27/03 see change earlier
    if(dcorr[1] == 1)
    {clims.df@mdc <- thpl.out@mdc}
    clims.df@conf.level <- conf
    #class(clims.df) <- c("calib", "data.frame")
    clims.df
}

#######################################################################
## Computes the calibration statistics for the linear standard curve ##
## modifications made by Matthew Mitchell on 3/27/03                 ##
##                                                                   ##
## derivatives for quadratic are incorrect - see loq.lin for more    ##
## details - MM 4/15/03                                              ##
## added xname argument on 8/30/04 MWM                               ##
#######################################################################

calib.lin.pom <- function(lin.out, y0, conf = 0.9, m = lin.out@m, dilution = 1, 
	samp.names, truth, times, samp.units = "", dose.units = "", 
	dose.name = "", 
	#toler = 0.001, 
	extrap = F, 
	#dos = lin.out@dos,
    rname="response",xname="x")
{
#	if(!dos) {
#		if(!is.loaded(symbol.C("linlow"))) {
#			cat("Loading linpom.o\n")
#			dyn.load("/local/lib/calib/sun4/src/linpom.o")
#		}
#	}
	if(missing(m))
		m <- lin.out@m
	b <- lin.out@coefficients
	n <- length(y0)
	alpha <- 1 - (1 - conf)/2
	tcrit <- qt(alpha, lin.out@df.residual)
	if(length(dilution) == 1)
		dilution <- rep(dilution, n)
	dcorr <- 1/dilution
	max.x <- max(lin.out@x)
	min.x <- min(lin.out@x, 0)
        
## MM added these next 2 lines on 3/27/03 because m may be different here
## than in lin.out and its these values of m that will determine the mdc
## This is added an output variable 

    cmdc <-  mdc.calc(lin.out, m=m,type="lin")
    cmdc <-  cmdc$mdc*dcorr
        
	if(lin.out@mmod == "lin") {
		cal.val <- (y0 - b[1])/b[2]
		resp.range <- c(b[1], b[1] + b[2] * max.x)
	}
     
        
	if(lin.out@mmod == "quad") {
		resp.range <- c(b[1], b[1] + b[2] * max.x + b[3] * 
			max.x * max.x)
		cal1 <- ( - b[2] - sqrt(b[2]^2 - 4 * b[3] * (b[1] - y0)
			))/(2 * b[3])
		cal2 <- ( - b[2] + sqrt(b[2]^2 - 4 * b[3] * (b[1] - y0)
			))/(2 * b[3])
		cal.val <- ifelse(cal1 > min(lin.out@x) & cal1 < max(
			lin.out@x), cal1, cal2)
	}
# now get the calibrated values
	xlow <- rep(min.x, n)
	xup <- rep(max.x, n)
	cov.un <- lin.out@cov.unscaled
	svec <- c(cov.un[, 1], cov.un[, 2])	
	# calculate  limits on lower pred. line
## 02-08-07: This call to C does not seem relavent. DVS
#	if(!dos) {
#		clims.lpred.out <- .C("linlow",
#			xl = as.double(xlow),
#			xu = as.double(xup),
#			as.double(b),
#			as.double(svec),
#			as.integer(5000),
#			as.double(toler),
#			as.double(lin.out@theta),
#			as.double(tcrit),
#			as.double(lin.out@sigma),
#			as.double(y0),
#			as.integer(n),
#			as.integer(m),
#			as.integer(length(b)))
#		if(b[2] < 0)
#			xlow <- clims.lpred.out$xl
#		if(b[2] > 0)
#			xup <- clims.lpred.out$xl
#		xsl <- rep(min.x, n)
#		xsu <- rep(max.x, n)
#		clims.upred.out <- .C("linup",
#			xl = as.double(xsl),
#			xu = as.double(xsu),
#			as.double(b),
#			as.double(svec),
#			as.integer(5000),
#			as.double(toler),
#			as.double(lin.out@theta),
#			as.double(tcrit),
#			as.double(lin.out@sigma),
#			as.double(y0),
#			as.integer(n),
#			as.integer(m),
#			as.integer(length(b)))
#		if(b[2] < 0)
#			xup <- clims.upred.out$xu
#		if(b[2] > 0)
#			xlow <- clims.upred.out$xu
#	}
#	else {
		xp <- seq(min.x, max.x, len = 300)
		if(lin.out@mmod == "quad") {
			ypred <- b[1] + b[2] * xp + b[3] * xp * xp
			grad <- cbind(1, xp, xp * xp)
		}
		if(lin.out@mmod == "lin") {
			ypred <- b[1] + b[2] * xp
			grad <- cbind(1, xp)
			resp.range <- c(b[1], b[1] + b[2] * max.x)
		}
    
		qn1 <- (ypred^(2 * lin.out@theta))/m
		qn2 <- diag(grad %*% cov.un %*% t(grad))
		ypl <- ypred + sign(b[2]) * tcrit * lin.out@sigma * 
			sqrt(qn1 + qn2)
		ypu <- ypred - sign(b[2]) * tcrit * lin.out@sigma * 
			sqrt(qn1 + qn2)
		xlow <- approx(ypl, xp, y0,ties = "ordered")$y
		xup <- approx(ypu, xp, y0,ties = "ordered")$y
#	}
## 05-14-07: Making some changes here, the old logic statement, below
## didn't seem to make much sense. DVS
#	if(sign(b[2] < 0)) {
	if(sign(b[2]) < 0) {
		xup <- ifelse(is.na(xup) | xup >= max.x, NA, xup)
		xup.inver <- ifelse(y0 < resp.range[2], NA, xup)
		xlow <- ifelse(y0 > resp.range[1], 0, xlow)
		xlow.inver <- ifelse(is.na(xlow) & is.finite(cal.val), 
			0, xlow)
		oor <- ifelse(y0 > resp.range[1] | y0 < resp.range[2], 
			"x", "o")
		too.low <- ifelse(y0 > resp.range[1], 1, 0)
		too.high <- ifelse(y0 < resp.range[2], 1, 0)
	}
	else {
		xup <- ifelse(xup >= max.x, NA, xup)
		xup.inver <- ifelse(is.na(xup) & y0 < resp.range[1], 0, 
			xup)
		xlow <- ifelse(y0 < resp.range[1], 0, xlow)
		xlow.inver <- ifelse(is.na(xlow) & y0 < resp.range[2], 
			0, xlow)
		oor <- ifelse(y0 < resp.range[1] | y0 > resp.range[2], 
			"x", "o")
		too.low <- ifelse(y0 < resp.range[1], 1, 0)
		too.high <- ifelse(y0 > resp.range[2], 1, 0)
	}
	hyb <- cal.val
	sigma2 <- lin.out@sigma * lin.out@sigma
##  This is incorrect - correct tcrit is at the beginning of the program
##  MM 4/14/03
##	tcrit <- qt((1 + alpha)/2, lin.out$df.res)
	if(lin.out@mmod == "lin") {
		dh.dy <- rep(1/b[2], n)
		dh.da <- rep(-1/b[2], n)
		dh.db <-  - (y0 - b[1])/(b[2] * b[2])
		var.xnot.hat <- ((dh.dy * dh.dy) * sigma2 * (y0^(2 * 
			lin.out@theta)))/m + sigma2 * (dh.da * dh.da * 
			cov.un[1, 1] + 2 * dh.da * dh.db * cov.un[1, 2] +
			dh.db * dh.db * cov.un[2, 2])
	}
        
	if(lin.out@mmod == "quad") {
##		dh.dy <- 1/(hyb + (b[2]/(2 * b[1])))
                dh.dy <- 1/(b[2] + 2*b[3]*hyb)
##		dh.da <- b[2]/(2 * b[1] * b[1]) + ((y0 - b[3])/b[1]) * (
##			1/(hyb + (b[2]/(2 * b[1])))) - (hyb + (b[2]/(2 * 
##			b[1])))/(2 * b[1] * b[1])
                dh.da <- -1/(b[2] + 2*b[3]*hyb)
##		dh.db <- -1/(2 * b[1]) + (b[2]/(2 * b[1])) * (1/(hyb + (
##			b[2]/(2 * b[1]))))
                dh.db <- -hyb/(b[2] + 2*b[3]*hyb)
##		dh.dc <- -1/(hyb + (b[2]/(2 * b[1])))
                dh.dc <- -(hyb*hyb)/(b[2] + 2*b[3]*hyb)
		var.xnot.hat <- ((dh.dy * dh.dy) * sigma2 * (y0^(2 * 
			lin.out@theta)))/m + sigma2 * (dh.da * (dh.da * 
			cov.un[1, 1] + dh.db * cov.un[1, 2] + dh.dc * 
			cov.un[1, 3]) + dh.db * (dh.da * cov.un[1, 2] + 
			dh.db * cov.un[2, 2] + dh.dc * cov.un[2, 3]) + 
			dh.dc * (dh.da * cov.un[1, 3] + dh.db * cov.un[
			2, 3] + dh.dc * cov.un[3, 3]))
	}
	se.xhat <- sqrt(var.xnot.hat)
	xup <- hyb + tcrit * se.xhat
	xlow <- hyb - tcrit * se.xhat
	if(!extrap){                                                     ## MM 9/27/03
	    cal.val <- ifelse(hyb < max.x, hyb, NA)
          xlow.wald <- ifelse(is.finite(xlow), ifelse(xlow < 0, 0, xlow), NA)
	    xup.wald <- ifelse(xup <= max.x, xup, NA)
	    }
		if(extrap){
		cal.val <- hyb
		xlow.wald <- xlow
		xup.wald <- xup
		}
       cal.val <- ifelse(cal.val > 0, cal.val, 0)  ## MM 8/30/04 - change negatives to 0
	cal.val <- dcorr * cal.val
	xlow.inver <- dcorr * xlow.inver
	xup.inver <- dcorr * xup.inver
	xlow.wald <- dcorr * xlow.wald
	xup.wald <- dcorr * xup.wald
	clims.df <- new("Calib",Estimated.x = as.vector(cal.val), PredStdErr = 
		as.vector(se.xhat), inver.low = as.vector(xlow.inver), 
		inver.up = as.vector(xup.inver), wald.low = as.vector(
		xlow.wald), wald.up = as.vector(xup.wald), avg.response = 
		as.vector(y0), dilution = as.vector(dilution), oor = 
		as.vector(oor))
#        names(clims.df)[1:2] <- c(paste("Estimated",xname),"PredStdErr")  ## MWM 8/30/04 added
#        names(clims.df)[match("response",names(clims.df))] <- paste("avg",rname,sep=".") ## MWM 9/2/04

       clims.df@wald.up <- ifelse(is.na(clims.df@wald.up) & !is.na(clims.df@wald.low),
                           clims.df@Estimated.x + (clims.df@Estimated.x-clims.df@wald.low),
                           clims.df@wald.up)   ## MWM 9/2/04
    
        
	lab.given <- F
	if(!missing(truth)) {
		ord.truth <- order(truth)
		y0 <- y0[ord.truth]
		cal.val <- cal.val[ord.truth]
		xlow <- xlow[ord.truth]
		xup <- xup[ord.truth]
		truth <- truth[ord.truth]
		labs <- list(name = "Control", plot.ind = truth, 
			too.low = too.low, too.high = too.high, xlabs
			 = as.character(truth), cont.pl = F, xunits = 
			samp.units, dose.units = dose.units, dose.name
			 = dose.name)
		lab.given <- T
		clims.df@truth <- truth
		if(missing(samp.names))
			clims.df@samp.names <- sapply(truth,
				function(a)
			paste("Cntl-", a, sep = ""))
		else clims.df@samp.names <- samp.names
#		clims.df <- clims.df[ord.truth,  ]
	}
	else clims.df@truth <- NULL
	if(!missing(times)) {
		clims.df@times <- times
		if(missing(samp.names))
			clims.df@samp.names <- times
		else clims.df@samp.names <- samp.names
		labs <- list(name = "Time", plot.ind = times, too.low
			 = too.low, too.high = too.high, xlabs = 
			as.character(times), cont.pl = T, xunits = 
			samp.units, dose.units = dose.units, dose.name
			 = dose.name)
		lab.given <- T
	}
	else clims.df@times <- NULL
	if(!lab.given) {
		labs <- list(name = "Sample", plot.ind = 1:length(y0), 
			too.low = too.low, too.high = too.high, xlabs
			 = as.character(1:length(y0)), cont.pl = F, 
			xunits = samp.units, dose.units = dose.units, 
			dose.name = dose.name)
		if(!missing(samp.names)) {
			clims.df@samp.names <- samp.names
			labs$xlabs <- samp.names
		}
		else clims.df@samp.names <- NULL
	}
	clims.df@labels <- labs
	clims.df@max.x <- max.x
## MM 9/29/03
      clims.df@extrap <- extrap
        
## MM added these 2 lines on 3/27/03 to indicate whether m here is the same
## as in the lin print out (this affects the MDC)
   ifelse(all(m == lin.out@m), {repeq <- T}, {repeq <- F})
   clims.df@repeq <- repeq
    
## MM added this line on 3/27/03 that the actual name of the response
## variable will be shown on the printout
   clims.df@rname <- rname
  
## MM 3/27/03 see change earlier
    if(dcorr[1] == 1)
    	clims.df@mdc <- lin.out@mdc
        
	clims.df@conf.level <- conf
#	class(clims.df) <- c("calib", "data.frame")
	clims.df
}

fpl.inverse <- function(fpl.out, y)
{
####################################################################
## Inverting the 4 Parameter Logistic Model to get calibrations   ##
## modified by Matthew Mitchell on 3/11/03                        ##
##                                                                ##
## input variables:                                               ##
##  fpl.out: The output from the nls fit                          ##
##  y: the y-coordinates                                          ##
##                                                                ##
## This function is called when performing the assay performance  ##
##  summary measure is the fpl function                           ##
####################################################################

  
  
###   	function returns x for b vector and y value for logistic eqn.
### mod: pdh - 8/22/91: temporary fix due to error in ^ function 
### 		for V3.0 Beta 3
	b <- fpl.out@coefficients
	a <- (b[1] - y)/(y - b[2])
	a <- ifelse(a < 0, NA, a)
	c <- 1/b[4]
	if(fpl.out@parm == 1)
		inver <- as.vector(b[3] * (a^c))
	else inver <- as.vector(exp(c * log(a) + b[3]))
	inver
}

thpl.inverse <- function(thpl.out, y)
{
####################################################################
## Inverting the 3 Parameter Logistic Model to get calibrations   ##
## written by Matthew Mitchell on 5/20/03                         ##
##                                                                ##
## input variables:                                               ##
##  thpl.out: The output from the nls fit                         ##
##  y: the y-coordinates                                          ##
##                                                                ##
## This function is called when performing the assay performance  ##
##  summary measure is the thpl function                          ##
####################################################################

  
  
###   	function returns x for b vector and y value for logistic eqn.
### mod: pdh - 8/22/91: temporary fix due to error in ^ function 
### 		for V3.0 Beta 3
	b <- thpl.out@coefficients
	a <- (b[1] - y)/(y - b[2])
	a <- ifelse(a < 0, NA, a)
        inver <- as.vector(b[3]*a)
	inver
}

