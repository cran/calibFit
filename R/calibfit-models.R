##############################################################################
## calibfit-models.R
##
## models used by calib.fit
## 
## Author: Haaland
###############################################################################




fpl.model <- function(x, b1, b2, b3, b4, w = 1, logParm = TRUE){
	####################################################################
	## The 4 Parameter Logistic Model used in Calibration Analysis    ##
	## modified by Matthew Mitchell on 3/11/03                        ##
	##                                                                ##
	## input variables:                                               ##
	##  x: the independent variable (usually dosage here)             ##
	##  b1,b2,b3,b4: values used in the FPL fit                       ##
	##  w: weights.  In OLS these are 1.                              ##
	##  parm:  If parm=1 the nonlinear function fitted is             ##
	##        y = b2 + (b1-b2)/(1 + (x/b3)^b4)                        ##
	##        If parm=2 the nonlinear function fitted is              ##
	##        y = b2 + (b1-b2)/(1 + exp(b4(logx - b3)),               ##       
	##        which is reparametrization of the other form            ##
	##        since b3 from alt=T is log(b3 from alt=F).              ##
	##                                                                ##
	## This function is called when doing the nls fitting with        ##
	##  the fpl and fpl.pom functions.                                ##
	##                                                                ## 
	## Original comments are ###, added are ##                        ##
	####################################################################
	if(length(b1) == 4) {
		b2 <- b1[2]
		b3 <- b1[3]
		b4 <- b1[4]
		b1 <- b1[1]
	}
	if(!logParm) {
		.a <- ifelse(x <= 0, 0, (x/b3)^b4)
		.den <- 1 + .a
		
		.value <- w * ((b1 - b2)/.den + b2)
	}
	else {
		.a <- ifelse(x <= 0, 0, exp(b4 * (log(x) - b3)))
		.den <- 1 + .a
		.value <- w * ((b1 - b2)/.den + b2)
	}
	.grad <- array(.value, c(length(.value), 4), list(NULL, c("b1", 
							"b2", "b3", "b4"))) 
	
	.grad[, "b1"] <- w * (1/.den)
	.grad[, "b2"] <- w * (1 - 1/.den)
	if(!logParm) {
		if(any(b3 <= 0))
			stop("Warning, b3 is less than or equal to zero.")
		.grad[, "b3"] <- w * (((b1 - b2) * (b4/b3) * .a)/.den^2)
		.grad[, "b4"] <- ifelse(x == 0, 0, (w * ( - (b1 - b2) * 
								.a * log(ifelse(x/b3 > 0, x/b3, NA))))/.den^2)
	}
	else {
		.grad[, "b3"] <- w * (((b1 - b2) * b4 * .a)/.den^2)
		.grad[, "b4"] <- ifelse(x == 0, 0, w * ((( - (b1 - b2) * 
									.a * (log(ifelse(x > 0, x, NA)) - b3))/.den^2))
		)
	}
	attr(.value, "gradient") <- .grad
	.value
}



thpl.model <- function(x, b1, b2, b3, w = 1, logParm = TRUE)
###################################################################
## The 3 parameter logistic model used in calibration analysis
##
###################################################################
{
	if(length(b1) == 3) {
		b2 <- b1[2]
		b3 <- b1[3]
		b1 <- b1[1] 
	}
	if(!logParm) {
		.a <- ifelse(x <= 0, 0, (x/b3))
		.den <- 1 + .a
		.value <- w * ((b1 - b2)/.den + b2)
	}
	else {
		.a <- ifelse(x <= 0, 0, exp(log(x) - b3))
		.den <- 1 + .a
		.value <- w * ((b1 - b2)/.den + b2)
	}
	.grad <- array(.value, c(length(.value), 3), list(NULL, c("b1", 
							"b2", "b3"))) 
	## rownames are NULL, columnames are are the parameter names
	.grad[, "b1"] <- w * (1/.den)
	.grad[, "b2"] <- w * (1 - 1/.den)
#	.grad[, "b3"] <- w * ((b1 - b2) * (.a/b3))/.den^2
	
	if(!logParm) {
		if(any(b3 <= 0))
			stop("Warning, b3 is less than or equal to zero.")
		.grad[, "b3"] <- w * (((b1 - b2) * (1/b3) * .a)/.den^2)
	}
	else {
		.grad[, "b3"] <- w * (((b1 - b2) * .a)/.den^2)
	}
	
	attr(.value, "gradient") <- .grad
	.value
}


##TODO: Figure out whether this is incorrect -- should be tpl, why is everyt
tpl.model <- function(x, b1, b2, b3, b4, w = 1, logParm = TRUE){
	####################################################################
	## The 4 Parameter Logistic Model used in Calibration Analysis    ##
	## modified by Matthew Mitchell on 3/11/03                        ##
	##                                                                ##
	## input variables:                                               ##
	##  x: the independent variable (usually dosage here)             ##
	##  b1,b2,b3,b4: values used in the FPL fit                       ##
	##  w: weights.  In OLS these are 1.                              ##
	##  parm:  If parm=1 the nonlinear function fitted is             ##
	##        y = b2 + (b1-b2)/(1 + (x/b3)^b4)                        ##
	##        If parm=2 the nonlinear function fitted is              ##
	##        y = b2 + (b1-b2)/(1 + exp(b4(logx - b3)),               ##       
	##        which is reparametrization of the other form            ##
	##        since b3 from alt=T is log(b3 from alt=F).              ##
	##                                                                ##
	## This function is called when doing the nls fitting with        ##
	##  the fpl and fpl.pom functions.                                ##
	##                                                                ## 
	## Original comments are ###, added are ##                        ##
	####################################################################
	
	if(length(b1) == 4) {
		b2 <- b1[2]
		b3 <- b1[3]
		b4 <- b1[4]
		b1 <- b1[1]
	}
	## If we are not working on the log scale
	if(!logParm) {
		.a <- ifelse(x <= 0, 0, (x/b3)^b4)
		.den <- 1 + .a
		
		.value <- w * ((b1 - b2)/.den + b2)
	}
	## If we are working on the log scale
	else {
		.a <- ifelse(x <= 0, 0, exp(b4 * (log(x) - b3)))
		.den <- 1 + .a
		.value <- w * ((b1 - b2)/.den + b2)
	}
	.grad <- array(.value, c(length(.value), 2), list(NULL, 
					c("b3", "b4"))) ## rownames are NULL, columnames are
	## are the parameter names
	if(!logParm) {
		if(any(b3 <= 0))
			stop("Warning, b3 is less than or equal to zero.")
		.grad[, "b3"] <- w * (((b1 - b2) * (b4/b3) * .a)/.den^2)
		.grad[, "b4"] <- ifelse(x == 0, 0, (w * ( - (b1 - b2) * 
								.a * log(ifelse(x > 0, x/b3, NA))))/.den^2)
	}
	else {
		.grad[, "b3"] <- w * (((b1 - b2) * b4 * .a)/.den^2)
		.grad[, "b4"] <- ifelse(x == 0, 0, w * ((( - (b1 - b2) * 
									.a * (log(ifelse(x > 0, x, NA)) - b3))/.den^2)))
	}
	attr(.value, "gradient") <- .grad
	.value
}


lin.model <- function(x, beta, w = 1, type)
{
	##################################################
	## called for the linear model fit in lin.pom   ##
	##                                              ##
	## header added by Matthew Mitchell on 3/18/03  ##
	## actually I don't this is used? 3/28/03 MM    ##
	##################################################
	
	if(type == "lin") {
		.value <- beta[1] + beta[2] * x
		.grad <- array(.value, c(length(.value), 2), list(NULL, 
						c("b1", "b2")))
		.grad[, "b1"] <- w
		.grad[, "b2"] <- w * x
	}
	
	if(type == "quad") {
		.value <- beta[1] + beta[2] * x + beta[3] * x * x
		.grad <- array(.value, c(length(.value), 3), list(NULL, 
						c("b1", "b2", "b3")))
		.grad[, "b1"] <- w
		.grad[, "b2"] <- w * x
		.grad[, "b3"] <- w * x * x
	}
	attr(.value, "gradient") <- .grad
	.value
}

