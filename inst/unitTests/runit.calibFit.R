##################################################################################
## unit tests are setup to run with R CMD Check
## results are for calib 2.02
# 
# Author: samarov, haaland
###############################################################################
load("data/RUnitcalibFit.rda")
ELISA <- read.table("../../data/ELISA.txt", header = TRUE)
HPLC <- read.table("../../data/HPLC.txt", header = TRUE)
library(RUnit)

####################################################################
## this is the setup for the RUnit testing and would happen only once
#model.fpl <- calib.fit(ELISA[,1], ELISA[,2], type = "fpl")
#model.fpl.pom <- calib.fit(ELISA[,1], ELISA[,2], type = "fpl.pom")
#model.log.fpl <- calib.fit(ELISA[,1], ELISA[,2], type = "log.fpl")
#model.log.fpl.pom <- calib.fit(ELISA[,1], ELISA[,2], type = "log.fpl.pom")
#model.thpl <- calib.fit(ELISA[,1], ELISA[,2], type = "thpl")
#model.thpl.pom <- calib.fit(ELISA[,1], ELISA[,2], type = "thpl.pom")
#model.tpl <- calib.fit(ELISA[,1], ELISA[,2], type = "tpl")
#model.tpl.pom <- calib.fit(ELISA[,1], ELISA[,2], type = "tpl.pom")
#model.log.tpl <- calib.fit(ELISA[,1], ELISA[,2], type = "log.tpl")
#model.log.tpl.pom <- calib.fit(ELISA[,1], ELISA[,2], type = "log.tpl.pom")
#model.lin <- calib.fit(HPLC[,1], HPLC[,2], type = "lin.pom")
#model.quad <- calib.fit(HPLC[,1], HPLC[,2], type = "quad.pom")
#save(model.fpl,model.fpl.pom,model.log.fpl,model.log.fpl.pom,model.thpl,model.thpl.pom,
#		model.tpl,model.tpl.pom,model.log.tpl,model.log.tpl.pom,model.lin,model.quad,
#		file="data/RUnitcalibFit.rda")


## the test functions use these results so evaluate outside
mod.fpl <- calib.fit(ELISA[,1], ELISA[,2], type = "fpl")
mod.fpl.pom <- calib.fit(ELISA[,1], ELISA[,2], type = "fpl.pom")
mod.log.fpl <- calib.fit(ELISA[,1], ELISA[,2], type = "log.fpl")
mod.log.fpl.pom <- calib.fit(ELISA[,1], ELISA[,2], type = "log.fpl.pom")
mod.thpl <- calib.fit(ELISA[,1], ELISA[,2], type = "thpl")
mod.thpl.pom <- calib.fit(ELISA[,1], ELISA[,2], type = "thpl.pom")
mod.tpl <- calib.fit(ELISA[,1], ELISA[,2], type = "tpl")
mod.tpl.pom <- calib.fit(ELISA[,1], ELISA[,2], type = "tpl.pom")
mod.log.tpl <- calib.fit(ELISA[,1], ELISA[,2], type = "log.tpl")
mod.log.tpl.pom <- calib.fit(ELISA[,1], ELISA[,2], type = "log.tpl.pom")
mod.lin <- calib.fit(HPLC[,1], HPLC[,2], type = "lin.pom")
mod.quad <- calib.fit(HPLC[,1], HPLC[,2], type = "quad.pom")


test.calibFit <- function(){	
	## fpl
	checkEqualsNumeric(mod.fpl@coefficients, model.fpl@coefficients, tolerance = 1e-2)
	checkEqualsNumeric(mod.fpl@se.coefficients, model.fpl@se.coefficients, tolerance = 1e-2)
	checkEqualsNumeric(mod.fpl@cov.unscaled, model.fpl@cov.unscaled, tolerance = 1e-2)
	checkEqualsNumeric(mod.fpl@x, model.fpl@x, tolerance = 1e-2)
	checkEqualsNumeric(mod.fpl@y, model.fpl@y, tolerance = 1e-2)
	checkEqualsNumeric(mod.fpl@mdc, model.fpl@mdc, tolerance = 1e-2)
	checkEqualsNumeric(mod.fpl@rdl, model.fpl@rdl, tolerance = 1e-2)
	checkEqualsNumeric(mod.fpl@fitted.values, model.fpl@fitted.values, tolerance = 1e-2)
	checkEqualsNumeric(mod.fpl@sigma, model.fpl@sigma, tolerance = 1e-2)
	checkEqualsNumeric(mod.fpl@cv, model.fpl@cv, tolerance = 1e-2)
	checkEquals(mod.fpl@theta, model.fpl@theta, tolerance = 1e-2)
	checkEqualsNumeric(mod.fpl@residuals, model.fpl@residuals, tolerance = 1e-2)
	
	## fpl.pom
	checkEqualsNumeric(mod.fpl.pom@coefficients, model.fpl.pom@coefficients, tolerance = 1e-2)
	checkEqualsNumeric(mod.fpl.pom@se.coefficients, model.fpl.pom@se.coefficients, tolerance = 1e-2)
	checkEqualsNumeric(mod.fpl.pom@cov.unscaled, model.fpl.pom@cov.unscaled, tolerance = 1e-2)
	checkEqualsNumeric(mod.fpl.pom@x, model.fpl.pom@x, tolerance = 1e-2)
	checkEqualsNumeric(mod.fpl.pom@y, model.fpl.pom@y, tolerance = 1e-2)
	checkEqualsNumeric(mod.fpl.pom@mdc, model.fpl.pom@mdc, tolerance = 1e-2)
	checkEqualsNumeric(mod.fpl.pom@rdl, model.fpl.pom@rdl, tolerance = 1e-2)
	checkEqualsNumeric(mod.fpl.pom@fitted.values, model.fpl.pom@fitted.values, tolerance = 1e-2)
	checkEqualsNumeric(mod.fpl.pom@sigma, model.fpl.pom@sigma, tolerance = 1e-2)
	checkEqualsNumeric(mod.fpl.pom@cv, model.fpl.pom@cv, tolerance = 1e-2)
	checkEqualsNumeric(mod.fpl.pom@theta, model.fpl.pom@theta, tolerance = 1e-2)
	checkEqualsNumeric(mod.fpl.pom@residuals, model.fpl.pom@residuals, tolerance = 1e-2)
	
	## log.fpl
	checkEqualsNumeric(mod.log.fpl@coefficients, model.log.fpl@coefficients, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.fpl@se.coefficients, model.log.fpl@se.coefficients, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.fpl@cov.unscaled, model.log.fpl@cov.unscaled, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.fpl@x, model.log.fpl@x, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.fpl@y, model.log.fpl@y, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.fpl@mdc, model.log.fpl@mdc, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.fpl@rdl, model.log.fpl@rdl, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.fpl@fitted.values, model.log.fpl@fitted.values, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.fpl@sigma, model.log.fpl@sigma, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.fpl@cv, model.log.fpl@cv, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.fpl@theta, model.log.fpl@theta, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.fpl@residuals, model.log.fpl@residuals, tolerance = 1e-2)
	
	## log.fpl.pom
	checkEqualsNumeric(mod.log.fpl.pom@coefficients, model.log.fpl.pom@coefficients, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.fpl.pom@se.coefficients, model.log.fpl.pom@se.coefficients, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.fpl.pom@cov.unscaled, model.log.fpl.pom@cov.unscaled, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.fpl.pom@x, model.log.fpl.pom@x, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.fpl.pom@y, model.log.fpl.pom@y, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.fpl.pom@mdc, model.log.fpl.pom@mdc, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.fpl.pom@rdl, model.log.fpl.pom@rdl, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.fpl.pom@fitted.values, model.log.fpl.pom@fitted.values, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.fpl.pom@sigma, model.log.fpl.pom@sigma, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.fpl.pom@cv, model.log.fpl.pom@cv, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.fpl.pom@theta, model.log.fpl.pom@theta, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.fpl.pom@residuals, model.log.fpl.pom@residuals, tolerance = 1e-2)
	
	## thpl
	checkEqualsNumeric(mod.thpl@coefficients, model.thpl@coefficients, tolerance = 1e-2)
	checkEqualsNumeric(mod.thpl@se.coefficients, model.thpl@se.coefficients, tolerance = 1e-2)
	checkEqualsNumeric(mod.thpl@cov.unscaled, model.thpl@cov.unscaled, tolerance = 1e-2)
	checkEqualsNumeric(mod.thpl@x, model.thpl@x, tolerance = 1e-2)
	checkEqualsNumeric(mod.thpl@y, model.thpl@y, tolerance = 1e-2)
	checkEqualsNumeric(mod.thpl@mdc, model.thpl@mdc, tolerance = 1e-2)
	checkEqualsNumeric(mod.thpl@rdl, model.thpl@rdl, tolerance = 1e-2)
	checkEqualsNumeric(mod.thpl@fitted.values, model.thpl@fitted.values, tolerance = 1e-2)
	checkEqualsNumeric(mod.thpl@sigma, model.thpl@sigma, tolerance = 1e-2)
	checkEqualsNumeric(mod.thpl@cv, model.thpl@cv, tolerance = 1e-2)
	checkEquals(mod.thpl@theta, model.thpl@theta, tolerance = 1e-2)
	checkEqualsNumeric(mod.thpl@residuals, model.thpl@residuals, tolerance = 1e-2)
	
	## thpl.pom
	checkEqualsNumeric(mod.thpl.pom@coefficients, model.thpl.pom@coefficients, tolerance = 1e-2)
	checkEqualsNumeric(mod.thpl.pom@se.coefficients, model.thpl.pom@se.coefficients, tolerance = 1e-2)
	checkEqualsNumeric(mod.thpl.pom@cov.unscaled, model.thpl.pom@cov.unscaled, tolerance = 1e-2)
	checkEqualsNumeric(mod.thpl.pom@x, model.thpl.pom@x, tolerance = 1e-2)
	checkEqualsNumeric(mod.thpl.pom@y, model.thpl.pom@y, tolerance = 1e-2)
	checkEqualsNumeric(mod.thpl.pom@mdc, model.thpl.pom@mdc, tolerance = 1e-2)
	checkEqualsNumeric(mod.thpl.pom@rdl, model.thpl.pom@rdl, tolerance = 1e-2)
	checkEqualsNumeric(mod.thpl.pom@fitted.values, model.thpl.pom@fitted.values, tolerance = 1e-2)
	checkEqualsNumeric(mod.thpl.pom@sigma, model.thpl.pom@sigma, tolerance = 1e-2)
	checkEqualsNumeric(mod.thpl.pom@cv, model.thpl.pom@cv, tolerance = 1e-2)
	checkEqualsNumeric(mod.thpl.pom@theta, model.thpl.pom@theta, tolerance = 1e-2)
	checkEqualsNumeric(mod.thpl.pom@residuals, model.thpl.pom@residuals, tolerance = 1e-2)
	
	## tpl
	checkEqualsNumeric(mod.tpl@coefficients, model.tpl@coefficients, tolerance = 1e-2)
	checkEqualsNumeric(mod.tpl@se.coefficients, model.tpl@se.coefficients, tolerance = 1e-2)
	checkEqualsNumeric(mod.tpl@cov.unscaled, model.tpl@cov.unscaled, tolerance = 1e-2)
	checkEqualsNumeric(mod.tpl@x, model.tpl@x, tolerance = 1e-2)
	checkEqualsNumeric(mod.tpl@y, model.tpl@y, tolerance = 1e-2)
	checkEqualsNumeric(mod.tpl@mdc, model.tpl@mdc, tolerance = 1e-1)
	checkEqualsNumeric(mod.tpl@rdl, model.tpl@rdl, tolerance = 1e-1)
	checkEqualsNumeric(mod.tpl@fitted.values, model.tpl@fitted.values, tolerance = 1e-2)
	checkEqualsNumeric(mod.tpl@sigma, model.tpl@sigma, tolerance = 1e-2)
	checkEqualsNumeric(mod.tpl@cv, model.tpl@cv, tolerance = 1e-2)
	checkEquals(mod.tpl@theta, model.tpl@theta, tolerance = 1e-2)
	checkEqualsNumeric(mod.tpl@residuals, model.tpl@residuals, tolerance = 1e-2)
	
	## tpl.pom
	checkEqualsNumeric(mod.tpl.pom@coefficients, model.tpl.pom@coefficients, tolerance = 1e-2)
	checkEqualsNumeric(mod.tpl.pom@se.coefficients, model.tpl.pom@se.coefficients, tolerance = 1e-1)
	checkEqualsNumeric(mod.tpl.pom@cov.unscaled, model.tpl.pom@cov.unscaled, tolerance = 1e-2)
	checkEqualsNumeric(mod.tpl.pom@x, model.tpl.pom@x, tolerance = 1e-2)
	checkEqualsNumeric(mod.tpl.pom@y, model.tpl.pom@y, tolerance = 1e-2)
	checkEqualsNumeric(mod.tpl.pom@mdc, model.tpl.pom@mdc, tolerance = 1e-1)
	checkEqualsNumeric(mod.tpl.pom@rdl, model.tpl.pom@rdl, tolerance = 1e-1)
	checkEqualsNumeric(mod.tpl.pom@fitted.values, model.tpl.pom@fitted.values, tolerance = 1e-2)
	checkEqualsNumeric(mod.tpl.pom@sigma, model.tpl.pom@sigma, tolerance = 1e-1)
	checkEqualsNumeric(mod.tpl.pom@cv, model.tpl.pom@cv, tolerance = 1e-2)
	checkEqualsNumeric(mod.tpl.pom@theta, model.tpl.pom@theta, tolerance = 1e-2)
	checkEqualsNumeric(mod.tpl.pom@residuals, model.tpl.pom@residuals, tolerance = 1e-2)
	
	## log.tpl
	checkEqualsNumeric(mod.log.tpl@coefficients, model.log.tpl@coefficients, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.tpl@se.coefficients, model.log.tpl@se.coefficients, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.tpl@cov.unscaled, model.log.tpl@cov.unscaled, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.tpl@x, model.log.tpl@x, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.tpl@y, model.log.tpl@y, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.tpl@mdc, model.log.tpl@mdc, tolerance = 1e-1)
	checkEqualsNumeric(mod.log.tpl@rdl, model.log.tpl@rdl, tolerance = 1e-1)
	checkEqualsNumeric(mod.log.tpl@fitted.values, model.log.tpl@fitted.values, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.tpl@sigma, model.log.tpl@sigma, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.tpl@cv, model.log.tpl@cv, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.tpl@theta, model.log.tpl@theta, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.tpl@residuals, model.log.tpl@residuals, tolerance = 1e-2)
	
	## log.tpl.pom
	checkEqualsNumeric(mod.log.tpl.pom@coefficients, model.log.tpl.pom@coefficients, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.tpl.pom@se.coefficients, model.log.tpl.pom@se.coefficients, tolerance = 1e-1)
	checkEqualsNumeric(mod.log.tpl.pom@cov.unscaled, model.log.tpl.pom@cov.unscaled, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.tpl.pom@x, model.log.tpl.pom@x, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.tpl.pom@y, model.log.tpl.pom@y, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.tpl.pom@mdc, model.log.tpl.pom@mdc, tolerance = 1e-1)
	checkEqualsNumeric(mod.log.tpl.pom@rdl, model.log.tpl.pom@rdl, tolerance = 1e-1)
	checkEqualsNumeric(mod.log.tpl.pom@fitted.values, model.log.tpl.pom@fitted.values, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.tpl.pom@sigma, model.log.tpl.pom@sigma, tolerance = 1e-1)
	checkEqualsNumeric(mod.log.tpl.pom@cv, model.log.tpl.pom@cv, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.tpl.pom@theta, model.log.tpl.pom@theta, tolerance = 1e-2)
	checkEqualsNumeric(mod.log.tpl.pom@residuals, model.log.tpl.pom@residuals, tolerance = 1e-2)
	
}

## Test that the accessors are acting in the way that they should
test.calibFitAccessors <- function(){

	## fpl	
	checkEquals(mod.fpl@mdc, mdc(mod.fpl))
	checkEquals(mod.fpl@rdl, rdl(mod.fpl))
	checkEquals(mod.fpl@loq, loq(mod.fpl))
	
	checkEquals(mod.log.fpl@mdc, mdc(mod.log.fpl))
	checkEquals(mod.log.fpl@rdl, rdl(mod.log.fpl))
	checkEquals(mod.log.fpl@loq, loq(mod.log.fpl))
	
	checkEquals(mod.fpl.pom@mdc, mdc(mod.fpl.pom))
	checkEquals(mod.fpl.pom@rdl, rdl(mod.fpl.pom))
	checkEquals(mod.fpl.pom@loq, loq(mod.fpl.pom))
		
	checkEquals(mod.log.fpl.pom@mdc, mdc(mod.log.fpl.pom))
	checkEquals(mod.log.fpl.pom@rdl, rdl(mod.log.fpl.pom))
	checkEquals(mod.log.fpl.pom@loq, loq(mod.log.fpl.pom))
	
	# tpl
	checkEquals(mod.tpl@mdc, mdc(mod.tpl))
	checkEquals(mod.tpl@rdl, rdl(mod.tpl))
	checkEquals(mod.tpl@loq, loq(mod.tpl))
	
	checkEquals(mod.log.tpl@mdc, mdc(mod.log.tpl))
	checkEquals(mod.log.tpl@rdl, rdl(mod.log.tpl))
	checkEquals(mod.log.tpl@loq, loq(mod.log.tpl))
	
	checkEquals(mod.tpl.pom@mdc, mdc(mod.tpl.pom))
	checkEquals(mod.tpl.pom@rdl, rdl(mod.tpl.pom))
	checkEquals(mod.tpl.pom@loq, loq(mod.tpl.pom))
	
	checkEquals(mod.log.tpl.pom@mdc, mdc(mod.log.tpl.pom))
	checkEquals(mod.log.tpl.pom@rdl, rdl(mod.log.tpl.pom))
	checkEquals(mod.log.tpl.pom@loq, loq(mod.log.tpl.pom))
			
	# thpl
	checkEquals(mod.thpl@mdc, mdc(mod.thpl))
	checkEquals(mod.thpl@rdl, rdl(mod.thpl))
	checkEquals(mod.thpl@loq, loq(mod.thpl))
	
	checkEquals(mod.thpl.pom@mdc, mdc(mod.thpl.pom))
	checkEquals(mod.thpl.pom@rdl, rdl(mod.thpl.pom))
	checkEquals(mod.thpl.pom@loq, loq(mod.thpl.pom))
	
	## lin
	checkEquals(mod.lin@mdc, mdc(mod.lin))
	checkEquals(mod.lin@rdl, rdl(mod.lin))
	checkEquals(mod.lin@loq, loq(mod.lin))
	
	checkEquals(mod.quad@mdc, mdc(mod.quad))
	checkEquals(mod.quad@rdl, rdl(mod.quad))
	checkEquals(mod.quad@loq, loq(mod.quad))
		
}