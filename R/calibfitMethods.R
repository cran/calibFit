## calibfitMethods.R
## 
## Accessor functions and other methods for calibfit objects
##
## Author: samarov
###############################################################################

## ============================================================================
## accessor method for slot mdc
## ============================================================================
setMethod("mdc",
	signature=c("calib.fit"),
	definition = function(object) object@mdc,
	valueClass="numOrNull")
## ============================================================================
## accessor method for slot rdl
## ============================================================================
setMethod("rdl",
	signature=c("calib.fit"),
	definition = function(object) object@rdl,
	valueClass="numOrNull")
## ============================================================================
## accessor method for slot loq
## ============================================================================	
setMethod("loq",
	signature=c("calib.fit"),
	definition = function(object) object@loq,
	valueClass="numOrNull")
## ============================================================================
## accessor method for slot coefficients
## ============================================================================
setMethod("coefficients",
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