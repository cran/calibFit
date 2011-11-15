# Contains accessor functions for calibfit related objects (i.e mdc, rdl, loq, etc.)
# 
# Author: samarov
###############################################################################

## ---------------------------------------------------------------------------
## Method for mdc
##============================================================================
setMethod("mdc", signature = "calib.fit",
		function(x,...){
			vars <- list(...)
			conf <- vars$conf
			if(is.null(conf)){
				if(!is.null(x@conf.level))
					conf <- x@conf.level
				else
					conf <- 0.95
			}
			m <- vars$m
			if(is.null(m)){
				if(!is.null(x@conf.level))
					m <- x@m
				else
					m <- 2
			}
			if((x@type == "lin") | (x@type == "quad"))
				mdc.out <- mdc.lin.pom(x, m = m, conf = conf)
			else if(x@type == "thpl")
				mdc.out <- mdc.thpl.pom(x, m = m, conf = conf)
			else if(x@type == "fpl")
				mdc.out <- mdc.fpl.pom(x, m = m, conf = conf)
			else if(x@type == "tpl")
				mdc.out <- mdc.tpl.pom(x, m = m, conf=conf)
			else
				stop("'type' not valid")
			return(mdc.out)
		})

## ---------------------------------------------------------------------------
## Method for rdl
##============================================================================
setMethod("rdl", signature = "calib.fit",
		function(x,...){
			vars <- list(...)
			conf <- vars$conf
			if(is.null(conf)){
				if(!is.null(x@conf.level))
					conf <- x@conf.level
				else
					conf <- 0.95
			}
			m <- vars$m
			if(is.null(m)){
				if(!is.null(x@conf.level))
					m <- x@m
				else
					m <- 2
			}
			if((x@type == "lin") | (x@type == "quad"))
				rdl.out <- rdl.lin.pom(x, m = m, conf = conf)
			else if(x@type == "thpl")
				rdl.out <- rdl.thpl.pom(x, m = m, conf = conf)
			else if(x@type == "fpl")
				rdl.out <- rdl.fpl.pom(x, m = m, conf = conf)
			else if(x@type == "tpl")
				rdl.out <- rdl.tpl.pom(x, m = m, conf = conf)
			else
				stop("'type' not valid")
			return(rdl.out)
		})

## ---------------------------------------------------------------------------
## accessor method for slot loq
## ============================================================================	
setMethod("loq",
		signature=c("calib.fit"),
		definition = function(x,...){
			vars <- list(...)
			m <- vars$m
			if(is.null(m)){
				if(!is.null(x@m))
					m <- x@m
				else
					m <- 2
			}
			cv <- vars$cv
			if(is.null(cv)){
				if(!is.null(x@cv))
					cv <- x@cv
				else
					cv <- .2
			}
			vlen <- vars$vlen
			if(is.null(vlen)){
				vlen <- 700
			}
			mit <- vars$mit
			if(is.null(mit)){
				mit <- 10000
			}
			toler <- vars$toler
			if(is.null(toler)){
				toler <- .001
			}
			
			if((x@type == "lin") | (x@type == "quad")){
				loq.out <- loq.lin.pom(x, m = m, cv = cv,
						vlen = vlen)
			}
			else if(x@type == "thpl"){
				loq.out <- loq.thpl.pom(x, m = m, cv = cv,
						vlen = vlen, mit = mit, toler = toler)
			}
			else if(x@type == "fpl"){
				loq.out <- loq.fpl.pom(x, m = m, cv = cv,
						vlen = vlen, mit = mit, toler = toler)
			}
			else if(x@type == "tpl"){
				loq.out <- loq.tpl.pom(x, m = m, cv = cv,
						vlen = vlen, mit = mit, toler = toler)
			}
			else
				stop("Not a valid object.")
			
			return(loq.out)
		})

##=============================================================================
## Print method. For the time being summary, show will be the same
##-----------------------------------------------------------------------------
setMethod("print", signature = "calib.fit", function(x, ...){
			cat("\nObject of class calib.fit\n")
			cat("Model type:",x@type,"\n")
			cat("Power of the mean variance (theta) used:",ifelse(x@pom,x@theta,"NA"),"\n")
			cat("Log parameterized:",x@logParm,"\n")
			cat("Method:",x@method,"\n")
			cat("\n")
			cat("Coefficients:\n")
			se.coef <- x@se.coefficients
			if(x@type == "tpl")
				se.coef <- c(NA,NA,x@se.coefficients)
			coef.df <- data.frame(rbind(x@coefficients,se.coef))
			names(coef.df) <- names(x@coefficients)
			rownames(coef.df) <- c("Estimates","Std. Error")
			print(coef.df)
			cat("\n")
			cat("Residual degrees of freedom:",x@df.residual,"\n")
			cat("Model Standard Deviation:",x@sigma,"\n")
			cat("Number of replicates m:",x@m,"\n")
			if(is.null(x@mdc))
				x@mdc <- NA
			if(is.null(x@rdl))
				x@rdl <- NA
			if(is.null(x@loq))
				x@loq <- NA
			diag.df <- data.frame(cbind(x@mdc,x@rdl,x@loq))
			names(diag.df) <- c("MDC","RDL","LOQ")
			rownames(diag.df) <- "Estimates"
			cat("\n")
			print(diag.df)
			cat("\n")
			invisible()
		})
setMethod("summary",signature="calib.fit",function(object,...){
			print(object,...)
		})
setMethod("show",signature="calib.fit",function(object){
			print(object)
		})

## ============================================================================
## accessor method for slot coefficients
## ============================================================================
setMethod("coefficients",
		signature=c("calib.fit"),
		definition = function(object) object@coefficients)

## ============================================================================
## accessor method for slot coefficients
## ============================================================================
setMethod("coef",
		signature=c("calib.fit"),
		definition = function(object) object@coefficients)


## ============================================================================
## accessor method for slot fitted
## ============================================================================
setMethod("fitted",
		signature=c("calib.fit"),
		definition = function(object) object@fitted.values)

## ============================================================================
## accessor method for slot residuals
## ============================================================================
setMethod("residuals",
		signature=c("calib.fit"),
		definition = function(object) object@residuals)

## ============================================================================
## accessor method for slot residuals
## ============================================================================
setMethod("resid",
		signature=c("calib.fit"),
		definition = function(object) object@residuals)