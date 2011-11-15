## calibplotMethods.R
## 
## Plot methods for calib.fit and calib objects
##
## Author: samarov
###############################################################################

## ============================================================================
## Plot methods for calib objects
## =======================================================================
setMethod("plot", signature(x = "calib", y = "missing"),
		function(x, ...) plot.calib(x, ...))

## ============================================================================
## Plot mehods for calib.fit objects
## =======================================================================
setMethod("plot", signature(x = "calib.fit", y = "missing"),
		function(x, type = NULL, ...){
			if(is.null(type))
				type <- "calib.fit" 
			if(type == "calib.fit"){
				if(x@type == "fpl")
					plot.fpl.pom(x, ...)
				else if(x@type == "thpl")
					plot.thpl.pom(x, ...)
				else if((x@type == "lin") | (x@type == "quad"))
					plot.lin.pom(x, ...)
				else if(x@type == "tpl")
					plot.tpl.pom(x, ...)
				else
					stop("Object type not valid.")
			}
			else if(type == "precprof"){
				if (x@type == "fpl")
					precprof.fpl.pom(x, ...)
				else if((x@type == "lin") | (x@type=="lin"))
					precprof.lin.pom(x, ...)
				else if(x@type == "thpl")
					precprof.thpl.pom(x, ...)
				else if(x@type == "tpl")
					precprof.tpl.pom(x, ...)
				else
					stop("Object type not valid.")
			}
			else if(type == "diagplot")
				diagplot(x,...)
			else
				stop("'type' not valid.")
			
		})
