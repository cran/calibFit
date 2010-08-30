## calibfit-functions.R
## 
## All class definitions for the package calib.
##
## Author: samarov
###############################################################################

b3start <- function(x, y, b1, b2, b4, logParm = TRUE){
	##===============================================================
	## Function for finding starting values for beta3
	## x: Typically concentration or dose.
	## y: OD of some sort.
	## b1: Supplied by estimate from startVals function.
	## b3: Supplied by estimate from startVals function.
	## b4: Supplied by estimate from startVals function.
	## logParm: Is x on the log scale.
	## ==============================================================
	ymid <- (b1 + b2)/2
	ydiff <- abs(y-ymid)
	
	mindiff <- min(ydiff, na.rm = TRUE)
	
	ymid <- (y[ydiff == mindiff])[1]
	xmid <- (x[ydiff == mindiff])[1]
	
	b3 <- xmid/(((b1 - b2)/(ymid - b2) - 1)^(1/b4))
	
	if(logParm)
		b3 <- log(b3)
	b3
}

startVals <- function(x, y, b1, b2, b3, b4, logParm = TRUE){
	##=================================================================
	## Function for finding starting values of the four 
	## the and two parameter logistic regression models.
	## x: Typically concentration or dose
	## y: OD of some sort
	## b1: If not provided will be calculated as below. In the
	##	   tpl model this is fixed.
	## b2: If not provided will be calculated as below. In the
	##	   tpl model this is fixed.
	## b3: If not provided will be calculated as below.
	## b4: If not provided will be calculated as below. In the thpl
	##	   model this is fixed.
	## logParm: Is x on the log scale. 
	##==================================================================
	
	## starting value for beta_1
	if(missing(b1))
		b1 <- mean(y[x==min(x)], na.rm = TRUE)
	else
		b1 <- b1
	
	## starting value for beta_2
	if(missing(b2))
		b2 <- mean(y[x==max(x, na.rm = TRUE)], na.rm = TRUE)
	else
		b2 <- b2
	
	## starting value for beta_4
	if(missing(b4))
		b4 <- 1
	else
		b4 <- b4
	
	## starting value for beta_3
	if(missing(b3))
		b3 <- b3start(x, y, b1 = b1, b2 = b2, b4 = b4, logParm = logParm)
	else
		b3 <- b3
	
	out <- list(b1 = b1, b2 = b2, b3 = b3, b4 = b4)
	return(out)
	
}

calib.fit <- function(x, y, 
		b1start, b2start, b3start, b4start, calcDiagnostics = TRUE,
		m, cv = 0.2, conf = 0.95, 
		mx = 50, lof.calc = T, lowLim = 1e-3,
		type=c("log.fpl.pom","fpl.pom","log.fpl","fpl","log.tpl.pom",
				"tpl.pom","log.tpl","tpl","log.thpl.pom",
				"thpl.pom","log.thpl","thpl","quad.pom","quad","lin.pom","lin"))
{
#	browser()
	####################################################################
	## Function for fitting a standard curve to the data. 			  ##
	##                                                                ##
	## input variables:                                               ##
	##  x: the independent variable (usually dosage here)             ##
	##  y: the response variable (example: OD, optical density)       ##
	##  m: This is the number of replicates (number of y's at         ##
	##     each x).                                                   ##
	##  cv: This is the acceptable coefficient of variation.          ##
	##      The limits of quantitation are calculated with            ##
	##      this constraint.                                          ##
	##  conf: The confidence level used for any prediction interval   ##
	##  b1start, b2start, b3start, b4start: the initial value used	  ## 
	##  mx: the maximum number of iterations used in the non-linear   ##
	##     least-squares fit.                                         ##
	##  lof.calc: if T, the lack of fit statistics are computed       ##
	##  lowLim: Any x-values less than zero are set to lowLim.		  ##
	##  type: Selects the type of model to be fit.					  ##
	####################################################################
	
	require(nlme)
	## Setting x values less than zero equal to a small positive constant to avoid
	## fitting issues.
	x[x == 0] <- lowLim
	
	## Available models
	# TODO: Make log, model, pom three separate choices, instead of combos? - parses them anyway
	availableModels <- c("log.fpl.pom","fpl.pom","log.fpl","fpl","log.tpl.pom",
			"tpl.pom","log.tpl","tpl","quad.pom","quad","log.thpl.pom",
			"thpl.pom","log.thpl","thpl","lin.pom","lin")
	## Checking to make sure only valid model have been selected
	if(any(!(type %in% availableModels)))
		stop(paste("'",paste(type[!(type %in% availableModels)], collapse = ", "),"'", 
						"are not valid models please select one or some subset of the following models\n", 
						paste(availableModels, collapse = ",\n")))
	
	## This statement parses out the argument for the model type to be
	## used. The syntax for the input statement of the model needs to
	## be as shown in the 'type' argument above. This is more for convenience purposes
	## regarding the automatic model selection procedure in R. 
	## logParm, modelType, and pom get set here
	modelSelect <- t(sapply(strsplit(type,"\\."), function(x){
						## If a log parameterized model is to be fit 'log' should always come first
						## in the statement, i.e. log.x.x
						if(x[1] == "log"){
							## log paramerization is set to be TRUE
							logParm <- TRUE
							## The model type is always the second x.model.x
							modelType <- x[2]
							## If a power of the mean model is to be fit this should come last
							## in the statement. 
							if(!is.na(x[3]))
								pom <- TRUE
							else
								pom <- FALSE
						}
						## Same as above but for the non-log model.
						else{
							logParm <- FALSE
							modelType <- x[1]
							if(!is.na(x[2]))
								pom <- TRUE
							else
								pom <- FALSE
						}
						
						as.matrix(c(I(logParm), modelType, I(pom)))
						
					}))
	
	## Model information
	modelSelect.df <- data.frame(matrix(NA, nrow = length(type), ncol = 3))
	names(modelSelect.df) <- c("logParm","modelType","pom")
	modelSelect.df$logParm <- as.logical(modelSelect[,1])
	modelSelect.df$modelType <- as.character(modelSelect[,2])
	modelSelect.df$pom <- as.logical(modelSelect[,3])
	
	## Setting up data frame
	dFrame <- data.frame(x = x, y = y)
	
	## The non-linear fit is performed
	
	for(i in 1:length(type)){
		
		## Model information
		logParm <- modelSelect.df[i, 1]
		modelType <- modelSelect.df[i, 2]
		pom <- modelSelect.df[i, 3]
		
		## Starting values
		sVals <- startVals(x, y, b1start, b2start, b3start, b4start, logParm = logParm)
		
		## Should a power of the mean variance model be used
		if(pom) {
			wts <- varPower()
		} else
			wts <- NULL
		
		## TODO: combine common elements of all the if{} sections by model below outside the if segments - reduce duplicate code
		## TODO: to reduce duplicate blocks, use parse/eval to set model statement specifically by type
		## TODO: use logParm directly instead of setting it in the if/else statement within each large if block - set lowLim either way?
		## TODO: check to see if try-error class could be looked at more specifically -- i.e. make sure it fails because it's not possible to fit, not for other reason
		
		## Fit a four parameter logistic model to the data	
		if(modelType == "fpl"){
			
			dFrame.i <- dFrame
			## If the model is log parameterized
			if(logParm){
				## If log parameterization is used x values less than or equal to zero are set to lowLim
				dFrame.i$x[dFrame.i$x <= 0] <- lowLim
				model <- try(gnls(y ~ fpl.model(x = x, b1 = b1, b2 = b2, b3 = b3, b4 = b4, w = 1, logParm = TRUE),
								data = dFrame.i, start = sVals, weights = wts,
								control = nls.control(maxiter = mx)), TRUE)
			}
			## If the model is not log parameterized
			else{
				model <- try(gnls(y ~ fpl.model(x = x, b1 = b1, b2 = b2, b3 = b3, b4 = b4, w = 1, logParm = FALSE),
								data = dFrame.i, start = sVals, weights = wts,
								control = nls.control(maxiter = mx)), TRUE)
			}
			## Check to see if we were able to fit a model to the data. If not keep cycling through.
#			browser()
			if(all(class(model) != "try-error")){
				if(i != 1)
					cat(paste("Warning", type[1], "produced error (", model[1], "), used", type[i], "instead"), "\n")
				break
			}
		}
		## Fit a three parameter logistic model to the data
		if(modelType == "thpl"){
			## Getting the starting values. These are essentially the same the fpl model but we
			## only use the b1, b2 and b3.
			sVals.i <- sVals[1:3]
			dFrame.i <- dFrame
			## If the model is log parameterized
			
			if(logParm){
				## If log parameterization is used x values less than or equal to zero are set to lowLim
				dFrame.i$x[dFrame.i$x <= 0] <- lowLim	
				model <- try(gnls(y ~ thpl.model(x, b1, b2, b3, logParm = TRUE),
								data = dFrame.i, start = sVals.i, weights = wts,
								control = nls.control(maxiter = mx)), TRUE)
			}
			else{
				## If the model is not log parameterized
				model <- try(gnls(y ~ thpl.model(x, b1, b2, b3, logParm = FALSE),
								data = dFrame.i, start = sVals.i, weights = wts,
								control = nls.control(maxiter = mx)), TRUE)				
			}
			## Check to see if we were able to fit a model to the data. If not keep cycling through.
			if(all(class(model) != "try-error")){
				if(i != 1)
					cat(paste("Warning", type[1], "produced error (", model[1], "), used", type[i], "instead"), "\n")
				break
			}
		}
		## Fit a two parameter logistic model to the data
		if(modelType == "tpl"){
			## For the two parameter logistic model b1 and b2 are fixed. If values are not provided 
			## the values will be set to the starting value estimates calculated from the startVals
			## function.
			if(missing(b1start)) {
				b1 <- sVals[[1]] 
			} else {
				b1 <- b1start
			}
			if(missing(b2start)) {
				b2 <- sVals[[2]]
			} else {
				b2 <- b2start
			}
			## The starting values for the tpl model. Only b3 and b4 are estimated.
			sVals.i <- sVals[3:4]
			dFrame.i <- dFrame
			## In order to get the gnls function to consider the parameters b1 and b2 to be fixed
			## need to include them as columns in the data frame.
			dFrame.i$b1 <- b1
			dFrame.i$b2 <- b2
			## If the model is log parameterized
			if(logParm){
				## If log parameterization is used x values less than or equal to zero are set to lowLim
				dFrame.i$x[dFrame.i$x <= 0] <- lowLim
				model <- try(gnls(y ~ tpl.model(x, b1, b2, b3, b4, logParm = TRUE),
								data = dFrame.i, start = sVals.i, weights = wts,
								control = nls.control(maxiter = mx)), TRUE)
			}
			## If the model is not log parameterized
			else{
				model <- try(gnls(y ~ tpl.model(x, b1, b2, b3, b4, logParm = FALSE),
								data = dFrame.i, start = sVals.i, weights = wts,
								control = nls.control(maxiter = mx)), TRUE)
			}
			## Check to see if we were able to fit a model to the data. If not keep cycling through.
			if(all(class(model) != "try-error")){
				if(i != 1)
					cat(paste("Warning", type[1], "produced error (", model[1], "), used", type[i], "instead"), "\n")
				break
			}
		}
		## Fit a quadratic model to the data
		if(modelType == "quad"){
			## For the quadtratic model need to include an additional quadratic term
			dFrame.i <- dFrame
			dFrame.i$x2 <- dFrame$x^2
			model <- try(gls(y ~ x + x2, weights = wts, data = dFrame.i, 
							control = glsControl(maxIter = mx)), TRUE)
			if(all(class(model) != "try-error")){
				if(i != 1)
					cat(paste("Warning", type[1], "produced error (", model[1], "), used", type[i], "instead"), "\n")
				break
			}
		}
		## Fit a linear model to the data
		if(modelType == "lin"){
			model <- try(gls(y ~ x, weights = wts, data = dFrame, 
							control = glsControl(maxIter = mx)), TRUE)
			if(all(class(model) != "try-error")){
				if(i != 1)
					cat(paste("Warning", type[1], "produced error (", model[1], "), used", type[i], "instead"), "\n")
				break
			}
			
		}
		## If no model can be found to fit the data then break.
		if((i == length(type)) & any(class(model) == "try-error"))
			stop("calib.fit is unable to fit a model to the data")
	}	
	
	## Getting the model coefficients
	b <- coef(model)
	
	## Getting theta
	theta <- coef(model$modelStruct)
	if(is.null(theta)) {
		theta <- 0  
	} else {
		theta <- theta
	}
	
	## If the model fit was tpl only b3 and b4 are fit, need to also return b1 and b2.
	if(modelType == "tpl"){
		b <- c(b1, b2, b)
		names(b) <- c("b1", "b2", "b3", "b4")
	}
	## Getting the fitted values
	fitVals <- as.numeric(fitted(model)) 
	
	## Getting weights from the model fit	
	wts <- 1/(fitVals^2)^theta
	
	## this is the gradient matrix
	if(modelType == "fpl")
		grad <- attr(fpl.model(x, b,  w = wts, logParm = logParm), "gradient")
	
	if(modelType == "thpl")
		grad <- attr(thpl.model(x, b, w = wts, logParm = logParm), "gradient")
	
	if(modelType == "tpl")
		grad <- attr(tpl.model(x, b,  w = wts, logParm = logParm), "gradient")
	
	if(modelType == "quad")
		grad <- attr(lin.model(x, b,  w = wts, type = modelType), "gradient")
	
	if(modelType == "lin")
		grad <- attr(lin.model(x, b,  w = wts, type = modelType), "gradient")
	
	## Determing the number of replicates
	if(missing(m))
		m <- length(x[x == min(x)])
	
	## Creating an object of class calib.fit   
	out <- new("calib.fit",
			coefficients = b, 
			se.coefficients = sqrt(diag(model$varBeta)),
			sigma = model$sigma, 
			cov.unscaled = ((1/model$sigma^2) * model$varBeta),
			df.residual = (length(x) - length(b)),
			residuals = model$residuals, 
			fitted.values = fitVals, 
			theta = theta,  
			pom = pom,
			method = model$method, 
			kused = model$numIter, 
			x = x, y = y, 
			m = m, cv = cv, 
			logParm = logParm, 
			conf.level = conf,
			gradient = grad, type = modelType)
	
	## Should model diagnostics be calculated
	if(calcDiagnostics) {
		
		out@mdc <- mdc(out, m = m, conf = conf)
		out@rdl <- rdl(out, m = m, conf = conf)
#		browser()
		## Doing a check on the RDL. If it falls outside the x range then
		## a warning is produced.
		if(!is.null(out@rdl) & !is.na(out@rdl)){
			if(out@rdl > max(x)) {
				out@rdlwarn <- "RDL outside range of values"
			} else {
				out@rdlwarn <- "RDL within range of values"
			}
		} else {
			out@rdlwarn <- "RDL returned NULL"
		}
		
		out@loq <- loq(out)
		out@conf.level <- conf
		
	}
	
	if(lof.calc) {
		if(length(x) > length(unique(x)))
			out@lof.test <- lof.test(out)
		else 
			out@lof.test <- NULL
	}
	
	return(out)
}
