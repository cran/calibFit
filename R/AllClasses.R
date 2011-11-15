## AllClasses.R
## 
## All class definitions for the package calib.
##
## Author: samarov
###############################################################################

## This class will allow a slot to take class numeric or NULL
setClass("numOrNULL")
setClassUnion("numOrNULL",c("numeric","logical","NULL"))

## This class will allow a slot to take class matrix or numeric
setClass("numOrMat")
setClassUnion("numOrMat",c("numeric","matrix"))

## This class will allow a slot to take class character or logical
setClass("charOrLogic")  
setClassUnion("charOrLogic",c("character","logical","NULL"))

setClass("charOrNULL")
setClassUnion("charOrNULL",c("character","NULL"))

## This class will allow a slot to take class numeric or integer
setClass("intOrNum")
setClassUnion("intOrNum",c("integer","numeric","NULL"))

## This class will allow a slot to take class numeric or integer
setClass("dataFrameOrNull")
setClassUnion("dataFrameOrNull",c("data.frame","NULL"))


## 07-28-08: I don't think that this should be its own class. It would be more appropriate to have this as a 
## data frame object. DVS
#setClass("lof.test",
#		representation(Fstat="numeric",
#				p.value="numeric",
#				lofss="numeric",
#				df.lof="numeric",
#				pure.error="numeric",
#				df.pure.error="numeric",
#				sse="numeric",
#				df.sse="intOrNum"))

setClass("calib.fit",
		representation(coefficients = "numeric",
				se.coefficients = "numeric",
				sigma = "numeric",
				cov.unscaled = "matrix",
				pom = "logical",
				theta = "numeric",
				df.residual = "intOrNum",
				fitted.values = "numeric",
				residuals = "numeric",
				method = "character",
				kused = "numeric",
				status = "character",
				x = "numOrMat",
				y = "numOrMat",
				logParm = "logical",
				m = "numeric",
				cv = "numeric",
				mdc = "numOrNULL",
				rdl = "numOrNULL",
				loq = "numOrNULL",
#				cf="numeric",
				gradient = "matrix",
				lof.test = "dataFrameOrNull",
				var.model = "character",
				conf.level = "numeric",
				## No long planning on using these
#				mmod="character",
#				SST="numeric",
#				SSE="numeric",
#				SSR="numeric",
#				Rsq="numeric",
				type = "character",
				rdlwarn = "character"),
				prototype = list(se.coefficients = numeric(0),
						se.coefficients = numeric(0),
						sigma = numeric(0),
						cov.unscaled = matrix(0),
						pom = logical(0),
						theta = numeric(0),
						df.residual = integer(0),
						fitted.values = numeric(0),
						method = character(0),
						kused = numeric(0),
						status = character(0),
						x = numeric(0),
						y = numeric(0),
						logParm = logical(0),
						m = numeric(0),
						cv = numeric(0),
						mdc = numeric(0),
						rdl = numeric(0),
						loq = numeric(0),
						gradient = matrix(0),
						lof.test = NULL,
						var.model = character(0),
						conf.level = numeric(0),
						type = character(0),
						rdlwarn = character(0)
						))

setClass("calib",
		representation(Estimated.x = "numeric",
				PredStdErr = "numeric",
				inver.low = "numOrNULL",
				inver.up="numOrNULL",
				wald.low="numOrNULL",
				wald.up="numOrNULL",
#				waldl.low="numOrNULL",
#				waldl.up="numOrNULL",
				avg.response="numeric",
				dilution="numeric",
				oor="character",
				lmdc="charOrLogic",
				row.names="character",
				labels="list",
				max.x="numeric",
				extrap="logical",
				repeq="logical",
				rname="character",
				conf.level="numeric",
				type="character",
				mdc="numOrNULL",
				truth="charOrNULL",
				times="charOrNULL",
				samp.names="charOrNULL"
#				lof.test="lof.test",
#				calib.fit="calib.fit"
				))
