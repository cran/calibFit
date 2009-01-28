calib <- function(x,y0,conf=.9, dilution=1, samp.names, m = x@m,
		truth, times, samp.units="", dose.units = "", dose.name = "", maxit = 1000,
		toler = 1e-05, rname="response", extrap = FALSE, xname = x){
	## 06-08-07: Removing the maxit and toler options. Both of these were used
	## previously in the source call to C but since this is not being used
	## anymore makes sense to take it out. DVS
	
	if((x@type == "lin") | (x@type == "quad")){
		calib.out <- calib.lin.pom(x, y0 = y0, conf = conf, 
				m = m, dilution = dilution, samp.names, truth, 
				times, samp.units = samp.units, dose.units = dose.units, 
				dose.name = dose.name, 
				extrap = F, 
				rname = rname, xname = xname)
	}
	else if(x@type == "thpl"){
		calib.out <- calib.thpl.pom(x, y0 = y0, conf = conf, 
				m = m, dilution = dilution, samp.names, truth, 
				times, samp.units = samp.units, dose.units = dose.units, 
				dose.name = dose.name, 
				extrap = F,
				rname = rname, xname = xname)
	}
	else if(x@type == "fpl"){
		calib.out <- calib.fpl.pom(x, y0 = y0, conf = conf, 
				m = m, dilution = dilution, samp.names, truth, 
				times, samp.units = samp.units, dose.units = dose.units, 
				dose.name = dose.name, 
				extrap = F,
				rname = rname, xname = xname)
	}
	else if(x@type == "tpl"){
		calib.out <- calib.tpl.pom(x, y0 = y0, conf = conf,
				m = m, dilution = dilution, samp.names, truth,
				times, samp.units = samp.units, dose.units = dose.units,
				dose.name = dose.name,
				extrap = F,
				rname = rname, xname = xname)
	}
	else
		stop("'type' not valid.")
	calib.out
}

calib.fpl.pom <- function(calibFit.out, y0, conf = calibFit.out@conf.level, 
		m = calibFit.out@m, dilution = 1, 
		samp.names, truth, times, samp.units = "", dose.units = "", 
		dose.name = "", rname="response", extrap = F, xname = "x")
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
	
	if(length(dilution) == 1)
		dilution <- rep(dilution, n)
	
	dcorr <- 1/dilution
	b <- calibFit.out@coefficients
	cov.un <- calibFit.out@cov.unscaled    
	min.x <- min(calibFit.out@x)
	max.x <- max(calibFit.out@x)     
	max.range <- fpl.model(max.x, b, logParm = calibFit.out@logParm) 
	min.range <- fpl.model(min.x, b, logParm = calibFit.out@logParm)
	if(sign(b[2] - b[1]) < 0) {
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
	
	cmdc <- mdc(calibFit.out, m = m) * dcorr
	cal.val <- fpl.inverse(calibFit.out, y0)
	tcrit <- qt((conf + 1)/2, calibFit.out@df.residual)
	
	seqy <- seq(fpl.model(min(calibFit.out@x), calibFit.out@coefficients, 
					logParm = calibFit.out@logParm), 
			fpl.model(max(calibFit.out@x), calibFit.out@coefficients, 
					logParm = calibFit.out@logParm), len = 6)
	boundx <- c(fpl.inverse(calibFit.out, seqy[1:5]), max.x)[-1]
	xp <- c(seq(0, boundx[1], len = 170), seq(boundx[1], 
					boundx[2], len = 50), seq(boundx[2], boundx[3], 
					len = 30), seq(boundx[3], boundx[4], len = 50), 
			seq(boundx[4], boundx[5], len = 170))
	ypred <- fpl.model(xp, b, logParm = calibFit.out@logParm)
	qn1 <- ((as.vector(ypred))^(2 * calibFit.out@theta))/m
	
	## 1-07-08: Note that the term '1/length(out@x)' was missing from below, i.e. the 1/n term
	## in the confidence/prediction interval calculation. DVS
	qn2 <- (1/length(calibFit.out@x)) * diag(attributes(ypred)$gradient %*% cov.un %*% t(attributes(ypred)$gradient))
	sigma <- calibFit.out@sigma
	ypl <- as.vector(ypred) + sign(b[2] - b[1]) * tcrit * sigma * sqrt(qn1 + qn2)
	ypu <- as.vector(ypred) - sign(b[2] - b[1]) * tcrit * sigma * sqrt(qn1 + qn2)
	xl <- try(approx(ypl, xp, y0, ties = "mean"), TRUE)
	if(class(xl) == "try-error")
		xl <- approx(ypl, xp, y0, ties = "ordered")
	xu <- try(approx(ypu, xp, y0, ties = "mean"), TRUE)
	if(class(xu) == "try-error")
		xu <- approx(ypu, xp, y0, ties = "ordered")
	xl <- xl$y
	xu <- xu$y
	
	if(sign(b[2] - b[1]) < 0) {
		if(!extrap){                                 
			cal.val <- ifelse(y0 > b[1], 0, cal.val)
			cal.val <- ifelse(y0 < b[2], NA, cal.val)
			cal.val <- ifelse(cal.val > max.x, NA, cal.val)
		}
		xup <- ifelse(is.na(xu) | xu >= max.x, NA, xu)
		xup <- ifelse(is.na(xup) & cal.val == 0, 0, xup)
		xup.inver <- ifelse(y0 < b[2], NA, xup)
		xlow <- ifelse(y0 > b[1], 0, xl)
		xlow.inver <- ifelse(is.na(xlow) & is.finite(cal.val), 0, xlow)
	}
	else {
		if(!extrap){                               
			cal.val <- ifelse(y0 < b[1], 0, cal.val)
			cal.val <- ifelse(y0 > b[2], NA, cal.val)
			cal.val <- ifelse(cal.val > max.x, NA, cal.val)
		}
		xup <- ifelse(xu >= max.x, NA, xu)
		xup.inver <- ifelse(is.na(xup) & cal.val == 0, 0, xup)
		xlow <- ifelse(y0 < b[1], 0, xl)    
		xlow.inver <- ifelse(is.na(xlow) & is.finite(cal.val), 0, xlow)
	}
	
	###########################################################################################
	## computes the Wald intervals which are computed with the delta method ##
	###########################################################################################
	if(!calibFit.out@logParm) {
		if(b[1] > b[2]){
			hyb <- rep(NA,length(y0))
			hyb[b[1] < y0] <- 0
			hyb[b[1] >= y0] <- b[3] * exp((1/b[4]) * log((b[1] - y0[b[1] >= y0])/(y0[b[1] >= y0] - b[2])))
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
		}
		## 05-14-07: Removed the following two lines since errors were
		## being produced in the ifelse statement
		#hyb <- ifelse(b[1] > y0, 0, b[3] * exp((1/b[4]) * 
		#        log((b[1] - y0)/(y0 - b[2]))))
		dh.dy <- hyb * ((b[2] - b[1])/(b[4] * (b[1] - y0) * (y0 -  b[2])))
		dh.db1 <- hyb/(b[4] * (b[1] - y0))
		dh.db2 <- hyb/(b[4] * (y0 - b[2]))
		dh.db3 <- hyb/b[3]
		dh.db4 <- rep(NA,length(y0))
		dh.db4[hyb == 0] <- 0
		dh.db4[hyb != 0] <- ( - hyb[hyb!=0]/(b[4] * b[4])) * 
				log((b[1] - y0[hyb!=0])/(y0[hyb!=0] - b[2]))
		
		## 05-14-07: Removed the following two lines since errors were
		## being produced in the ifelse statement
		#    dh.db4 <- ifelse(hyb == 0, 0, ( - hyb/(b[4] * b[4])) * 
		#        log((b[1] - y0)/(y0 - b[2])))
	}
	if(calibFit.out@logParm) {
		if(b[1] > b[2]){ 
			hyb <- rep(NA,length(y0))
			hyb[b[1] < y0] <- 0
			hyb[b[1] >= y0] <- exp((1/b[4]) * log((b[1] - y0[b[1] >= y0])/
									(y0[b[1] >= y0] - b[2])) + b[3])
		}
		## 05-14-07: Removed the following two lines since errors were
		## being produced in the ifelse statement
		#hyb <- ifelse(b[1] < y0, 0, exp((1/b[4]) * 
		#	log((b[1] - y0)/(y0 - b[2])) + b[3]))
		else{ 
			hyb <- rep(NA,length(y0))
			hyb[b[1] > y0] <- 0
			hyb[b[1] <= y0] <- exp((1/b[4]) * log((b[1] - y0[b[1] <= y0])/
									(y0[b[1] <= y0] - b[2])) + b[3])
		}
		
		dh.dy <- (hyb * (b[2] - b[1]))/(b[4] * (b[1] - y0) * (y0 - b[2]))
		dh.db1 <- hyb/(b[4] * (b[1] - y0))
		dh.db2 <- hyb/(b[4] * (y0 - b[2]))
		dh.db3 <- hyb
		dh.db4 <- rep(NA,length(y0))
		dh.db4[hyb == 0] <- 0
		dh.db4[hyb != 0] <- ( - hyb[hyb != 0]/(b[4] * b[4])) * 
				log((b[1] - y0[hyb != 0])/(y0[hyb != 0] - b[2]))
	}
	dh.df <- cbind(dh.db1, dh.db2, dh.db3, dh.db4)
	xnot.hat <- hyb
	sigma2 <- sigma^2
	
	var.xnot.hat <- sigma2 * (dh.dy^2 * (y0^2)^calibFit.out@theta/m  + diag(dh.df %*% cov.un %*% t(dh.df)))
	tcrit <- qt((1 + conf)/2, calibFit.out@df.residual)
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
	if(!calibFit.out@logParm) {
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
		
		var.hybl <- (dh.dyl * dh.dyl) * sigma2 * (y0^(2 * 
						calibFit.out@theta)) + sigma2 * (dh.db1l * (dh.db1l * 
						cov.un[1, 1] + dh.db2l * cov.un[2, 1] + dh.db3l *
						cov.un[3, 1] + dh.db4l * cov.un[4, 1]) + 
					dh.db2l * (dh.db1l * cov.un[1, 2] + dh.db2l * 
						cov.un[2, 2] + dh.db3l * cov.un[3, 2] + dh.db4l *
						cov.un[4, 2]) + dh.db3l * (dh.db1l * cov.un[1, 3] + dh.db2l * 
						cov.un[2, 3] + dh.db3l * cov.un[3, 3] + dh.db4l * cov.un[4, 3]) + 
					dh.db4l * (dh.db1l * cov.un[1, 4] + dh.db2l * cov.un[2, 4] +
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
#	browser()
	## next line added by MM on 3/27/03 to compute when the
	## estimated calibration is less than the mdc
	lmdc <- ifelse(as.vector(cal.val) < as.vector(cmdc), "X", "O")
	
	############################################################
	## creates a data set with all the calibration statistics ##
	############################################################
	
	if(any(is.null(xup.waldl))) { ## 05-14-07: Note, made a change from
		## is.na to is.null. DVS
		clims.df <- new("calib",Estimated.x = as.vector(cal.val), 
				PredStdErr = as.vector(se.xhat), 
				inver.low = as.vector(xlow.inver), 
				inver.up = as.vector(xup.inver), 
				wald.low = as.vector(xlow.wald), 
				wald.up = as.vector(xup.wald), 
				avg.response = as.vector(y0), 
				dilution = as.vector(dilution), 
				oor = as.vector(oor), mdc = as.vector(cmdc),
				lmdc = as.vector(lmdc))
	}
	else {
		clims.df <- new("calib",Estimated.x = as.vector(cal.val), 
				PredStdErr = as.vector(se.xhat), 
				inver.low = as.vector(xlow.inver), 
				inver.up = as.vector(xup.inver), 
				wald.low = as.vector(xlow.wald), 
				wald.up = as.vector(xup.wald), 
				avg.response = as.vector(y0), 
				waldl.low = as.vector(xlo.waldl),
				waldl.up = as.vector(xup.waldl), 
				dilution = as.vector(dilution), 
				oor = as.vector(oor), mdc = as.vector(cmdc),
				lmdc = as.vector(lmdc))
	}
	
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
	}
	
	##########################################################################################################
	## used if times are provided                           
	##########################################################################################################
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
	
	##########################################################################################################
	## will be executed is truth and times are missing      
	##########################################################################################################
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
	ifelse(all(m == calibFit.out@m), {repeq <- T}, {repeq <- F})
	clims.df@repeq <- repeq
	
	## MM added this line on 3/27/03 that the actual name of the response
	## variable will be shown on the printout
	clims.df@rname <- rname
	
	## MM 3/27/03 see change earlier
	if(dcorr[1] == 1){
		clims.df@mdc <- calibFit.out@mdc
	}
	clims.df@conf.level <- conf
#	clims.df@calib.fit <- calibFit.out
	clims.df
}

calib.thpl.pom <- function(calibFit.out, y0, conf = calibFit.out@conf.level, 
		m = calibFit.out@m, dilution = 1, 
		samp.names, truth, times, samp.units = "", dose.units = "", 
		dose.name = "", rname="response", extrap = FALSE, xname = "x"){
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
	
	if(length(dilution) == 1)
		dilution <- rep(dilution, n)
	dcorr <- 1/dilution
	b <- calibFit.out@coefficients
	cov.un <- calibFit.out@cov.unscaled
	min.x <- min(calibFit.out@x)
	max.x <- max(calibFit.out@x)
	max.range <- thpl.model(max.x, b) 
	min.range <- thpl.model(min.x, b)
	if(sign(b[2] - b[1]) < 0) {
		oor <- ifelse(y0 > min.range[1] | y0 < max.range[1], "x", "o")
		too.low <- ifelse(y0 > min.range[1], 1, 0)
		too.high <- ifelse(y0 < max.range[1], 1, 0)
	}
	if(sign(b[2] - b[1] > 0)) {
		oor <- ifelse(y0 < min.range[1] | y0 > max.range[1], "x", "o")
		too.low <- ifelse(y0 < min.range[1], 1, 0)
		too.high <- ifelse(y0 > max.range[1], 1, 0)
	}
	
	## MM added these next 2 lines on 3/27/03 because m may be different here
	## than in thpl.out and its these values of m that will determine the mdc
	## This is added an output variable 
	
	cmdc <-  mdc(calibFit.out, m = m, conf = conf)
	cmdc <-  cmdc * dcorr
	
	cal.val <- thpl.inverse(calibFit.out, y0)
	cal.val <- as.array(cal.val)    
	### correct calibrated value for dilution
	cal.val <- cal.val
	tcrit <- qt((conf + 1)/2, calibFit.out@df.residual)
	
	####################################################################################
	## computes the inverse limits here 
	####################################################################################
	
	seqy <- seq(thpl.model(min(calibFit.out@x), calibFit.out@coefficients, 
					logParm = calibFit.out@logParm),
			thpl.model(max(calibFit.out@x), calibFit.out@coefficients, 
					logParm = calibFit.out@logParm), len = 6)
	boundx <- c(thpl.inverse(calibFit.out, seqy[1:5]), max.x)[-1]
	xp <- c(seq(0, boundx[1], len = 170), 
			seq(boundx[1], boundx[2], len = 50), 
			seq(boundx[2], boundx[3], len = 30), 
			seq(boundx[3], boundx[4], len = 50), 
			seq(boundx[4], boundx[5], len = 170))
	ypred <- thpl.model(xp, b)
	qn1 <- ((as.vector(ypred))^(2 * calibFit.out@theta))/m
	
	## 1-07-08: Note that the term '1/length(out@x)' was missing from below, i.e. the 1/n term
	## in the confidence/prediction interval calculation. DVS
	qn2 <- (1/length(calibFit.out@x)) * diag(attributes(ypred)$gradient %*% cov.un %*% t(
					attributes(ypred)$gradient))
	sigma <- calibFit.out@sigma
	ypl <- as.vector(ypred) + sign(b[2] - b[1]) * tcrit * 
			sigma * sqrt(qn1 + qn2)
	ypu <- as.vector(ypred) - sign(b[2] - b[1]) * tcrit * 
			sigma * sqrt(qn1 + qn2)
	xl <- try(approx(ypl, xp, y0, ties = "mean"), TRUE)
	if(class(xl) == "try-error")
		xl <- approx(ypl, xp, y0, ties = "ordered")
	xu <- try(approx(ypu, xp, y0, ties = "mean"), TRUE)
	if(class(xu) == "try-error")
		xu <- approx(ypu, xp, y0, ties = "ordered")
	xl <- xl$y
	xu <- xu$y
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
			cal.val <- ifelse(cal.val > max.x, NA, cal.val)
		}
		xup <- ifelse(xu >= max.x, NA, xu)
		xup.inver <- ifelse(is.na(xup) & cal.val == 0, 0, xup)
		xlow <- ifelse(y0 < b[1], 0, xl)    
		xlow.inver <- ifelse(is.na(xlow) & is.finite(cal.val), 0, xlow)
	}
	##############################################
	## computes the Wald intervals              ##
	## which are computed with the delta method ##
	##############################################
	
	if(b[1] > b[2])
		hyb <- ifelse(b[1] < y0, 0, b[3]*((b[1] - y0)/(y0 - b[2])))
	else 
		hyb <- ifelse(b[1] > y0,0, b[3]*((b[1] - y0)/(y0 - b[2]))) 
	dh.dy <- hyb * ((b[2] - b[1])/( (b[1] - y0) * (y0 -  b[2])))
	dh.db1 <- hyb/(b[1] - y0)
	dh.db2 <- hyb/(y0 - b[2])
	dh.db3 <- hyb/b[3]
	
	xnot.hat <- hyb
	sigma2 <- calibFit.out@sigma^2
	var.xnot.hat <- ((dh.dy^2) * sigma2 * (y0^(2 * calibFit.out@theta)))/m + 
			sigma2 * (dh.db1 * (dh.db1 * cov.un[1, 1] + 
					dh.db2 * cov.un[2, 1] + dh.db3 * cov.un[3, 1]) +
				dh.db2 * (dh.db1 * cov.un[1, 2] + 
					dh.db2 * cov.un[2, 2] + dh.db3 * cov.un[3, 2]) +
				dh.db3 * (dh.db1 * cov.un[1, 3] + 
					dh.db2 * cov.un[2, 3] + dh.db3 * cov.un[3, 3]))
	
	tcrit <- qt((1 + conf)/2, calibFit.out@df.residual)
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
					calibFit.out@theta)) + sigma2 * (dh.db1l * (dh.db1l * 
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
		clims.df <- new("calib", Estimated.x = as.vector(cal.val), 
				PredStdErr = as.vector(se.xhat), inver.low = 
						as.vector(xlow.inver), inver.up = as.vector(
						xup.inver), wald.low = as.vector(xlow.wald), 
				wald.up = as.vector(xup.wald), avg.response = 
						as.vector(y0), dilution = as.vector(dilution), 
				oor = as.vector(oor), mdc = as.vector(cmdc),
				lmdc = as.vector(lmdc))
	}
	else {
		clims.df <- new("calib", Estimated.x  = as.vector(cal.val), 
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
		labs <- list(name = "Control", plot.ind = 1:length(truth), 
				too.low = too.low[ord.truth], too.high = too.high[ord.truth],     
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
	else 
		clims.df@truth <- NULL
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
	else 
		clims.df@times <- NULL
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
		else 
			clims.df@samp.names <- NULL
	}
	clims.df@labels <- labs
	clims.df@max.x <- max.x
	
	## MWM 9/30/04
	clims.df@extrap <- extrap
	
	## MM added these 2 lines on 3/27/03 to indicate whether m here is the same
	## as in the thpl print out (this affects the MDC)
	ifelse(all(m == calibFit.out@m), {repeq <- T}, {repeq <- F})
	clims.df@repeq <- repeq
	
	## MM added this line on 3/27/03 that the actual name of the response
	## variable will be shown on the printout
	clims.df@rname <- rname
	
	## MM 3/27/03 see change earlier
	if(dcorr[1] == 1)
		clims.df@mdc <- calibFit.out@mdc
	clims.df@conf.level <- conf
	
	clims.df
}



calib.tpl.pom <- function(calibFit.out, y0, conf = 0.9, m = calibFit.out@m, dilution = 1, 
		samp.names, truth, times, samp.units = "", dose.units = "", 
		dose.name = "", 
		rname="response", extrap=FALSE, xname = "x")
{
	
	###############################################################################
	## Computes the Calibration Statistics                                       ##
	##                                                                           ##
	## input variables:                                                          ##
	##  tpl.out: output from the tpl or tpl.pom functions                        ##
	##  y0: a vector of the MEAN response values for the unknown x               ##
	##  conf: confidence level for calibration estimates                         ##
	##  m: number of terms in the average for y0                                 ##
	##     BE CAREFUL: by default this will read in the same number as           ##
	##     in tpl.out - the means for y0 may be based on a different             ##
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
	
	## Extracting necessary variables for calibration calculation
	n <- length(y0)
	theta <- calibFit.out@theta
	logParm <- calibFit.out@logParm
	sigma <- calibFit.out@sigma
	b <- calibFit.out@coefficients
	cov.un <- calibFit.out@cov.unscaled    
	min.x <- min(calibFit.out@x)
	max.x <- max(calibFit.out@x)     
	max.range <- tpl.model(max.x, b, logParm = logParm) 
	min.range <- tpl.model(min.x, b, logParm = logParm)
	tcrit <- qt((conf + 1)/2, calibFit.out@df.residual)
	
	if(length(dilution) == 1)
		dilution <- rep(dilution, n)
	dcorr <- 1/dilution
	
	## This checks to see if the values we are looking to get the calibrated values for are outside the
	## range of x values, i.e. are we trying to extrapolate outside the sample.
	if(b[2] < b[1]) {
		oor <- ifelse(y0 > min.range[1] | y0 < max.range[1], "x", "o")
		too.low <- ifelse(y0 > min.range[1], 1, 0)
		too.high <- ifelse(y0 < max.range[1], 1, 0)
	}
	if(b[2] > b[1]) {
		oor <- ifelse(y0 < min.range[1] | y0 > max.range[1], "x", "o")
		too.low <- ifelse(y0 < min.range[1], 1, 0)
		too.high <- ifelse(y0 > max.range[1], 1, 0)
	}
	## MM added these next 2 lines on 3/27/03 because m may be different here
	## than in tpl.out and its these values of m that will determine the mdc
	## This is added an output variable 
	## Getting the MDC
	
	mdc.val <- mdc(calibFit.out, m = m, conf = conf) * dcorr
	cal.val <- tpl.inverse(calibFit.out, y0)
	
	### correct calibrated value for dilution  ##next line does nothing??
	seqy <- seq(tpl.model(min.x, b, logParm = logParm), 
			tpl.model(max.x, b, logParm = logParm), len = 6)
	boundx <- c(tpl.inverse(calibFit.out, seqy[1:5]), max.x)[-1]
	xp <- c(seq(0, boundx[1], len = 170), 
			seq(boundx[1], boundx[2], len = 50), 
			seq(boundx[2], boundx[3], len = 30), 
			seq(boundx[3], boundx[4], len = 50), 
			seq(boundx[4], boundx[5], len = 170))
	ypred <- tpl.model(xp, b, logParm = logParm)
	qn1 <- ((ypred^2)^theta)/m
	## 1-07-08: Note that the term '1/length(out@x)' was missing from below, i.e. the 1/n term
	## in the confidence/prediction interval calculation. DVS
	qn2 <- (1/length(calibFit.out@x)) * diag(attributes(ypred)$gradient %*% cov.un %*% t(attributes(ypred)$gradient))
	
	ypl <- ypred + sign(b[2] - b[1]) * tcrit * sigma * sqrt(qn1 + qn2)
	ypu <- ypred - sign(b[2] - b[1]) * tcrit * sigma * sqrt(qn1 + qn2)
	xl <- try(approx(ypl, xp, y0, ties = "mean"), TRUE)
	if(class(xl) == "try-error")
		xl <- approx(ypl, xp, y0, ties = "ordered")
	xu <- try(approx(ypu, xp, y0, ties = "mean"), TRUE)
	if(class(xu) == "try-error")
		xu <- approx(ypu, xp, y0, ties = "ordered")
	xl <- xl$y
	xu <- xu$y
	
	if(b[2] < b[1]){
		## Whether calibration should be extrapolated to x's outside the range
		if(!extrap){                                 
			cal.val <- ifelse(y0 > b[1], 0, cal.val)
			cal.val <- ifelse(y0 < b[2], NA, cal.val)
			cal.val <- ifelse(cal.val > max.x, NA, cal.val)
		}
		xup <- ifelse(is.na(xu) | xu >= max.x, NA, xu)
		xup <- ifelse(is.na(xup) & cal.val == 0, 0, xup)
		xup.inver <- ifelse(y0 < b[2], NA, xup)
		xlow <- ifelse(y0 > b[1], 0, xl)
		xlow.inver <- ifelse(is.na(xlow) & is.finite(cal.val), 0, xlow)
	}
	else {
		if(!extrap){                               
			cal.val <- ifelse(y0 < b[1], 0, cal.val)
			cal.val <- ifelse(y0 > b[2], NA, cal.val)
			cal.val <- ifelse(cal.val > max.x, NA, cal.val)
		}
		xup <- ifelse(xu >= max.x, NA, xu)
		xup.inver <- ifelse(is.na(xup) & cal.val == 0, 0, xup)
		xlow <- ifelse(y0 < b[1], 0, xl)    
		xlow.inver <- ifelse(is.na(xlow) & is.finite(cal.val), 0, xlow)
	}
	
	############################################################################################
	## Computes the Wald intervals which are computed with the delta method 
	############################################################################################
	if(!logParm) {
		if(b[1] > b[2]){
			hyb <- rep(NA,length(y0))
			hyb[(b[1] < y0)|(b[2] >= y0)] <- 0
			hyb[(b[1] >= y0)&(b[2] < y0)] <- b[3] * exp((1/b[4]) *
							log((b[1] - y0[(b[1] >= y0) & (b[2] < y0)])/(y0[(b[1] >= y0) & (b[2] < y0)] - b[2])))
		}
		else {
			hyb <- rep(NA,length(y0))
			hyb[(b[1] > y0)|(b[2] <= y0)] <- 0
			hyb[(b[1] <= y0)&(b[2] > y0)] <- b[3] * exp((1/b[4]) * 
							log((b[1] - y0[(b[1] <= y0)&(b[2] > y0)])/(y0[(b[1] <= y0)&(b[2] > y0)] - b[2])))
		}
		dh.dy <- hyb * ((b[2] - b[1])/(b[4] * (b[1] - y0) * (y0 -  b[2])))
		dh.db3 <- hyb/b[3]
		dh.db4 <- rep(NA,length(y0))
		dh.db4[hyb==0] <- 0
		dh.db4[hyb!=0] <- ( - hyb[hyb!=0]/(b[4] * b[4])) * 
				log((b[1] - y0[hyb!=0])/(y0[hyb!=0] - b[2]))
	}
	if(logParm){
		if(b[1] > b[2]){ 
			hyb <- rep(NA,length(y0))
			hyb[(b[1] < y0)|(b[2] >= y0)] <- 0
			hyb[(b[1] >= y0)&(b[2] < y0)] <- exp((1/b[4]) * log((b[1] - y0[(b[1] >= y0)&(b[2] < y0)])/
									(y0[(b[1] >= y0)&(b[2] < y0)] - b[2])) + b[3])
		}
		else{ 
			hyb <- rep(NA,length(y0))
			by12 <- (b[1] >= y0)|(b[2] <= y0)
			hyb[by12] <- 0
			hyb[!by12] <- exp((1/b[4]) * log((b[1] - y0[!by12])/(y0[!by12] - b[2])) + b[3])
		}
		dh.dy <- (hyb * (b[2] - b[1]))/(b[4] * (b[1] - y0) * (y0 - b[2]))
		dh.db3 <- hyb
		dh.db4 <- rep(NA,length(y0))
		dh.db4[hyb == 0] <- 0
		dh.db4[hyb != 0] <- ( - hyb[hyb != 0]/(b[4] * b[4])) * log((b[1] - y0[hyb != 0])/(y0[hyb != 0] - b[2]))
	}
	xnot.hat <- hyb
	dh.db <- cbind(dh.db3, dh.db4)
	
	var.xnot.hat <- sigma^2 * (dh.dy^2 * ((y0^2)^theta)/m + apply(dh.db, 1, function(z) t(z) %*% cov.un %*% z))
	
	se.xhat <- sqrt(var.xnot.hat)
	xup <- xnot.hat + tcrit * se.xhat
	xlow <- xnot.hat - tcrit * se.xhat
	if(!extrap) 
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
	##################################################################################################
	## Another set of confidence intervals for when the model is not log parameterized
	##################################################################################################
	if(!logParm) {
		if(b[1] > b[2]){ 
			hybl <- rep(NA,length(y0))
			hybl[(b[1] < y0)|(b[2] >= y0)] <- 0
			hybl[(b[1] >= y0)&(b[2] < y0)] <- log(b[3]) + (1/b[4]) * 
					log((b[1] - y0[(b[1] >= y0)&(b[2] < y0)])/(y0[(b[1] >= y0)&(b[2] < y0)] - b[2]))
		}
		else{ 
			hybl <- rep(NA,length(y0))
			hybl[(b[1] > y0) | (b[2] <= y0)] <- 0
			hybl[(b[1] <= y0) & (b[2] > y0)] <- log(b[3]) + (1/b[4]) * 
					log((b[1] - y0[(b[1] <= y0) & (b[2] > y0)])/(y0[(b[1] <= y0)&(b[2] > y0)] - b[2]))
		}		
		## 05-14-07: Was getting errors from the next line, added above
		## to account for shape of curve and avoid taking the log of
		## negative values. DVS
		
		dh.dyl <-  - (1/(b[4] * (b[1] - y0)) + 1/(b[4] * (y0 - b[2])))
		dh.db3l <- 1/b[3]
		
		if(b[1] > b[2]){ 
			dh.db4l <- rep(NA,length(y0))
			dh.db4l[(b[1] < y0)|(b[2] >= y0)] <- 0
			dh.db4l[(b[1] >= y0)&(b[2] < y0)] <- (1/(b[4] * b[4])) * 
					log((y0[(b[1] >= y0)&(b[2] < y0)] - b[2])/(b[1] - y0[(b[1] >= y0)&(b[2] < y0)]))
		}
		else{ 
			dh.db4l <- rep(NA,length(y0))
			dh.db4l[(b[1] >= y0) | (b[2] <= y0)] <- 0
			dh.db4l[(b[1] < y0)&(b[2] > y0)] <- (1/(b[4] * b[4])) * 
					log((y0[(b[1] < y0)&(b[2] > y0)] - b[2])/(b[1] - y0[(b[1] < y0)&(b[2] > y0)]))
		}		
		
		var.hybl <- (dh.dyl * dh.dyl) * sigma^2 * (y0^2)^theta + sigma^2 * (dh.db3l * (dh.db3l * cov.un[
								1, 1] + dh.db4l * cov.un[2, 1]) + dh.db4l * (
						dh.db3l * cov.un[1, 2] + dh.db4l * cov.un[2, 2]))
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
	lmdc <- ifelse(as.vector(cal.val) < as.vector(mdc.val), "X", "O")
	
	############################################################
	## creates a data set with all the calibration statistics ##
	############################################################
	
	if(any(is.null(xup.waldl))) {
		calib.out <- new("calib",
				Estimated.x = cal.val,
				PredStdErr = se.xhat, 
				inver.low = xlow.inver, 
				inver.up = xup.inver, 
				wald.low = xlow.wald, 
				wald.up = xup.wald, 
				avg.response = y0, 
				dilution = dilution, 
				oor = oor, 
				mdc = mdc.val,
				lmdc = lmdc)
	}
	else {
		calib.out <- new("calib",
				Estimated.x = cal.val, 
				PredStdErr = se.xhat, 
				inver.low = xlow.inver, 
				inver.up = xup.inver, 
				wald.low = xlow.wald, 
				wald.up = xup.wald, 
				avg.response = y0, 
				waldl.low = xlo.waldl,
				waldl.up = xup.waldl, 
				dilution = dilution, 
				oor = oor, 
				mdc = mdc.val,
				lmdc = lmdc)
	}
	
	calib.out@wald.up <- ifelse(is.na(calib.out@wald.up) & !is.na(calib.out@wald.low),
			calib.out@Estimated.x + (calib.out@Estimated.x - calib.out@wald.low),
			calib.out@wald.up)  
	
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
		labs <- list(name = "Control", plot.ind = 1:length(truth),
				too.low = too.low[ord.truth], too.high = too.high[ord.truth],     
				### i think its here that the truth (| samp names ?) get set as 
				### xlabs which turn into tick labels on axis in plot.calib
				### order truth is screwing things up for plot
				### subscripting too.low and too.high as above fixed it !! (moc 5/8/95)
				xlabs = as.character(truth), cont.pl = F, xunits = samp.units, 
				dose.units = dose.units, dose.name = dose.name)
		lab.given <- T
		calib.out@truth <- truth
		if(missing(samp.names))
			calib.out@samp.names <- truth
		else 
			calib.out@samp.names <- samp.names
		
	}
	
	## This line is not necessary as the slot will simply be empty
	
	##########################################################
	## used if times are provided                           ##
	##########################################################
	if(!missing(times)) {
		times <- as.character(times)  ## Added by DVS 10/9/06
		calib.out@times <- times
		if(missing(samp.names))
			calib.out@samp.names <- times
		else 
			calib.out@samp.names <- samp.names
		labs <- list(name = "Time", plot.ind = times, too.low = too.low, 
				too.high = too.high, xlabs = as.character(times), 
				cont.pl = T, xunits = samp.units, dose.units = dose.units, 
				dose.name = dose.name)
		lab.given <- T
	}
	
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
			calib.out@samp.names <- samp.names
			labs$xlabs <- samp.names
			if(is.numeric(samp.names))
				labs$plot.ind <- samp.names
		}
		# else attributes(clims.df)$samp.names <- NULL
	}
	calib.out@labels <- labs
	calib.out@max.x <- max.x
	
	## MM 9/29/03
	calib.out@extrap <- extrap
	
	## MM added these 2 lines on 3/27/03 to indicate whether m here is the same
	## as in the tpl print out (this affects the MDC)
	ifelse(all(m == calibFit.out@m), {repeq <- T}, {repeq <- F})
	calib.out@repeq <- repeq
	
	## MM added this line on 3/27/03 that the actual name of the response
	## variable will be shown on the printout
	calib.out@rname <- rname
	
	## MM 3/27/03 see change earlier
	if(dcorr[1] == 1){
		calib.out@mdc <- calibFit.out@mdc
	}
	calib.out@conf.level <- conf
	calib.out
}


#######################################################################
## Computes the calibration statistics for the linear standard curve ##
## modifications made by Matthew Mitchell on 3/27/03                 ##
##                                                                   ##
## derivatives for quadratic are incorrect - see loq.lin for more    ##
## details - MM 4/15/03                                              ##
## added xname argument on 8/30/04 MWM                               ##
#######################################################################

calib.lin.pom <- function(out, y0, conf = 0.95, m = out@m, dilution = 1, 
		samp.names, truth, times, samp.units = "", dose.units = "", 
		dose.name = "", 
		extrap = F, 
		rname="response",xname="x")
{
	
	if(missing(m))
		m <- out@m
	b <- out@coefficients
	n <- length(y0)
	alpha <- 1 - (1 - conf)/2
	tcrit <- qt(alpha, out@df.residual)
	if(length(dilution) == 1)
		dilution <- rep(dilution, n)
	dcorr <- 1/dilution
	max.x <- max(out@x)
	min.x <- min(out@x, 0)
	
	## MM added these next 2 lines on 3/27/03 because m may be different here
	## than in lin.out and its these values of m that will determine the mdc
	## This is added an output variable 
	
	cmdc <-  mdc(out, m = m, conf = conf)
	cmdc <-  cmdc * dcorr
	
	if(out@type == "lin") {
		cal.val <- (y0 - b[1])/b[2]
		resp.range <- c(b[1], b[1] + b[2] * max.x)
	}
	
	
	if(out@type == "quad") {
		resp.range <- c(b[1], b[1] + b[2] * max.x + b[3] * max.x * max.x)
		cal1 <- ( - b[2] - sqrt(b[2]^2 - 4 * b[3] * (b[1] - y0)))/(2 * b[3])
		cal2 <- ( - b[2] + sqrt(b[2]^2 - 4 * b[3] * (b[1] - y0)))/(2 * b[3])
		cal.val <- ifelse(cal1 > min(out@x) & cal1 < max(out@x), cal1, cal2)
	}
# now get the calibrated values
	xlow <- rep(min.x, n)
	xup <- rep(max.x, n)
	cov.un <- out@cov.unscaled
	svec <- c(cov.un[, 1], cov.un[, 2])	
	
	xp <- seq(min.x, max.x, len = 300)
	if(out@type == "quad") {
		ypred <- b[1] + b[2] * xp + b[3] * xp * xp
		grad <- cbind(1, xp, xp * xp)
	}
	if(out@type == "lin") {
		ypred <- b[1] + b[2] * xp
		grad <- cbind(1, xp)
		resp.range <- c(b[1], b[1] + b[2] * max.x)
	}
	
	qn1 <- ((ypred^2)^out@theta)/m
	## 1-07-08: Note that the term '1/length(out@x)' was missing from below, i.e. the 1/n term
	## in the confidence/prediction interval calculation. DVS
	qn2 <- (1/length(out@x)) * diag(grad %*% cov.un %*% t(grad))
	ypl <- ypred + sign(b[2]) * tcrit * out@sigma * 
			sqrt(qn1 + qn2)
	ypu <- ypred - sign(b[2]) * tcrit * out@sigma * 
			sqrt(qn1 + qn2)
	xlow <- try(approx(ypl, xp, y0, ties = "mean"), TRUE)
	if(class(xlow) == "try-error")
		xlow <- approx(ypl, xp, y0, ties = "ordered")
	xup <- try(approx(ypu, xp, y0, ties = "mean"), TRUE)
	if(class(xup) == "try-error")
		xup <- approx(ypu, xp, y0, ties = "ordered")
	xlow <- xlow$y
	xup <- xup$y
	
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
	sigma2 <- out@sigma * out@sigma
	##  This is incorrect - correct tcrit is at the beginning of the program
	##  MM 4/14/03
	##	tcrit <- qt((1 + alpha)/2, lin.out$df.res)
	if(out@type == "lin") {
		dh.dy <- rep(1/b[2], n)
		dh.da <- rep(-1/b[2], n)
		dh.db <-  - (y0 - b[1])/(b[2] * b[2])
		var.xnot.hat <- ((dh.dy * dh.dy) * sigma2 * (y0^(2 * 
							out@theta)))/m + sigma2 * (dh.da * dh.da * 
					cov.un[1, 1] + 2 * dh.da * dh.db * cov.un[1, 2] +
					dh.db * dh.db * cov.un[2, 2])
	}
	
	if(out@type == "quad") {
		
		dh.dy <- 1/(b[2] + 2 * b[3] * hyb)
		dh.da <- -1/(b[2] + 2 * b[3] * hyb)
		dh.db <- -hyb/(b[2] + 2 * b[3] * hyb)
		dh.dc <- -(hyb^2)/(b[2] + 2 * b[3] * hyb)
		var.xnot.hat <- ((dh.dy * dh.dy) * sigma2 * (y0^(2 * 
							out@theta)))/m + sigma2 * (dh.da * (dh.da * 
						cov.un[1, 1] + dh.db * cov.un[1, 2] + dh.dc * 
						cov.un[1, 3]) + dh.db * (dh.da * cov.un[1, 2] + 
						dh.db * cov.un[2, 2] + dh.dc * cov.un[2, 3]) + 
					dh.dc * (dh.da * cov.un[1, 3] + dh.db * cov.un[2, 3] + 
						dh.dc * cov.un[3, 3]))
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
	
	clims.df <- new("calib",Estimated.x = as.vector(cal.val), PredStdErr = 
					as.vector(se.xhat), inver.low = as.vector(xlow.inver), 
			inver.up = as.vector(xup.inver), wald.low = as.vector(
					xlow.wald), wald.up = as.vector(xup.wald), avg.response = 
					as.vector(y0), dilution = as.vector(dilution), oor = 
					as.vector(oor))
	
	clims.df@wald.up <- ifelse(is.na(clims.df@wald.up) & !is.na(clims.df@wald.low),
			clims.df@Estimated.x + (clims.df@Estimated.x-clims.df@wald.low),
			clims.df@wald.up)  
	
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
		else 
			clims.df@samp.names <- NULL
	}
	clims.df@labels <- labs
	clims.df@max.x <- max.x
	## MM 9/29/03
	clims.df@extrap <- extrap
#	browser()
	## MM added these 2 lines on 3/27/03 to indicate whether m here is the same
	## as in the lin print out (this affects the MDC)
	ifelse(all(m == out@m), {repeq <- T}, {repeq <- F})
	clims.df@repeq <- repeq
	
	## MM added this line on 3/27/03 that the actual name of the response
	## variable will be shown on the printout
	clims.df@rname <- rname
	
	## MM 3/27/03 see change earlier
	if(dcorr[1] == 1)
		clims.df@mdc <- out@mdc
	
	clims.df@conf.level <- conf
#	class(clims.df) <- c("calib", "data.frame")
	clims.df
}

fpl.inverse <- function(out, y)
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
	b <- out@coefficients
	a <- (b[1] - y)/(y - b[2])
	a <- ifelse(a < 0, NA, a)
	c <- 1/b[4]
	if(out@logParm)
		## No long used
#	if(out@parm == 1)
		inver <- as.vector(exp(c * log(a) + b[3]))
	else 
		inver <- as.vector(b[3] * (a^c))
	
	return(inver)
}

thpl.inverse <- function(out, y)
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
	b <- out@coefficients
	a <- (b[1] - y)/(y - b[2])
	a <- ifelse(a < 0, NA, a)
	
	## 07-01-08: For some reason there is no log parameterized inverse. DVS
	if(out@logParm)
		inver <- as.vector(exp(log(a) + b[3]))
	else 
		inver <- as.vector(b[3] * a)
#	inver <- as.vector(b[3] * a)
	return(inver)
}


tpl.inverse <- function(out, y)
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
	b <- out@coefficients
	a <- (b[1] - y)/(y - b[2])
	a <- ifelse(a < 0, NA, a)
	c <- 1/b[4]
#	if(out@parm == 1)
	if(out@logParm)
		inver <- as.vector(exp(c * log(a) + b[3]))
	else
		inver <- as.vector(b[3] * (a^c)) 
	
	return(inver)
}
