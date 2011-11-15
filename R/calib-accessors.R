# calib accessor functions
# 
# Author: samarov
###############################################################################

##=============================================================================
## print method for calib
##-----------------------------------------------------------------------------
setMethod("print", signature = "calib", function(x,...){
			cat("\nAn object of class calib\n")
			
			df = data.frame(Average.Response = x@avg.response,
					Estimated.Concentration = x@Estimated.x,
					Predicted.Standard.Error=x@PredStdErr,
					Inverse.Lower.CI=x@inver.low,
					Inverse.Upper.CI=x@inver.up,
					Wald.Lower.CI=x@wald.low,
					Wald.Upper.CI=x@wald.up)
			print(df)
#			cat("Predicted Standard Error:",x@PredStdErr,"\n")
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