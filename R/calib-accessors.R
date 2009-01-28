# calib accessor functions
# 
# Author: samarov
###############################################################################

##=============================================================================
## print method for calib
##-----------------------------------------------------------------------------
setMethod("print", signature = "calib", function(x,...){
			cat("\nAn object of class calib\n")
#			cat("Predicted Stadandard Error:",x@PredStdErr,"\n")
#			ci.df <- data.frame(rbind(c(x@inver.low,x@inver.up),
#							c(x@wald.low,x@wald.up)))
#			names(ci.df) <- c("Lower","Upper")
#			rownames(ci.df) <- c("Inverse Intervals", "Wald Intervals")
#			cat("Confidence Intervals:\n")
#			print(ci.df)
			cat("\n")
		})

setMethod("summary",signature="calib",function(object,...){
			print(object,...)
		})
setMethod("show",signature="calib",function(object){
			print(object)
		})