###################################################
### chunk number 1: loadLibrary
###################################################
library(calib)


###################################################
### chunk number 2: loadData
###################################################
data(HPLC)
data(ELISA)


###################################################
### chunk number 3: nameData
###################################################
conc.hplc <- HPLC[,1]
resp.hplc <- HPLC[,2]
conc.elisa <- ELISA[,1]
resp.elisa <- ELISA[,2]


###################################################
### chunk number 4: plotDataAndLinFit
###################################################
## Plot of the data with a std OLS fit
par(mfrow=c(1,2))
plot(conc.hplc,resp.hplc,
	xlab = "Concentration (ng/ml)",
	ylab = "Response",
	main = "HPLC data")
linmodel <- lm(resp.hplc~conc.hplc)
# The predicted response
linPredResp <- fitted(linmodel)
# Linear regression fit
lines(conc.hplc,linPredResp)

## Plot of the data with a std FPL fit
#par(mar = c(4,4,2,2))
plot(log(conc.elisa),resp.elisa,
	xlab = "log(Concentration (ng/ml))",
	ylab = "Response",
	main = "ELISA data")
fplmodel <- calib.fit(conc.elisa,resp.elisa,type="log.fpl")
# The predicted response
fplPredResp <- fplmodel@fitted.values
# fpl regression fit
lines(log(conc.elisa),fplPredResp)


###################################################
### chunk number 5: plotLinRes
###################################################
par(mfrow=c(1,2))
## Residuals from linear fit
linres <- residuals(linmodel)/summary(linmodel)$sigma
plot(linPredResp,linres,
	 xlab = "Predicted Value of Mean (LS)",
	 ylab = "Standardized Residuals",
	 main = "HPLC data",
	 ylim = c(-5,5))
abline(h=0)

## Residuals from fpl fit
fplres <- fplmodel@residuals/fplmodel@sigma
plot(fplPredResp,fplres,
	xlab = "Predicted Value of Mean (FPL)",
	ylab = "Standardized Residuals",
	main = "ELISA data",
	ylim = c(-5,5))
abline(h=0)


###################################################
### chunk number 6: logHPLCplot
###################################################
#par(mar = c(3.5,3.5,1.5,1.5))
plot(conc.hplc,log(resp.hplc),
	 xlab = "Concentration (ng/ml)",
	 ylab = "log(Response)")


###################################################
### chunk number 7: sqrtHPLCplot
###################################################
#par(mar = c(3.5,3.5,1.5,1.5))
plot(conc.hplc,sqrt(resp.hplc),
	 xlab = "Concentration (ng/ml)",
	 ylab = "sqrt(Response)")


###################################################
### chunk number 8: logAbsLogPred
###################################################
par(mfrow=c(1,2))
#par(mar = c(3.5,3.5,1.5,1.5))
plot(log(linPredResp),log(abs(linres)),
	xlab = "Log(LS predicted values)",
	ylab = "Log(absolute LS residuals)",
	main = "HPLC data",
	ylim = c(-6,2))
linresmodel <- lm(log(abs(linres))~log(linPredResp))
lines(log(linPredResp),fitted(linresmodel))

#par(mar = c(3.5,3.5,1.5,1.5))
plot(log(fplPredResp),log(abs(fplres)),
	xlab = "Log(FPL predicted values)",
	ylab = "Log(absolute FPL residuals)",
	main = "ELISA data",
	ylim = c(-6,2))
fplresmodel <- lm(log(abs(fplres))~log(fplPredResp))
lines(log(fplPredResp),fitted(fplresmodel))


###################################################
### chunk number 9: calib_fit
###################################################
cal.fpl <- calib.fit(conc.elisa,resp.elisa,type="log.fpl")
cal.lin.pom <- calib.fit(conc.hplc,resp.hplc,type="lin.pom")
cal.fpl.pom <- calib.fit(conc.elisa,resp.elisa,type="log.fpl.pom")

linpom.fit <- cal.lin.pom@fitted.values
fplpom.fit <- cal.fpl.pom@fitted.values

sig.lin <- cal.lin.pom@sigma
sig.fpl <- cal.fpl.pom@sigma

theta.lin <- cal.lin.pom@theta
theta.fpl <- cal.fpl.pom@theta

linpom.res <- cal.lin.pom@residuals*(1/((linpom.fit^theta.lin)*sig.lin))
fplpom.res <- cal.fpl.pom@residuals*(1/((fplpom.fit^theta.fpl)*sig.fpl))


###################################################
### chunk number 10: vferesid
###################################################
par(mfrow=c(1,2))
plot(linpom.fit,linpom.res,
	xlab = "Fitted Values (LS)",
	ylab = "Standardized Residuals",
	main = "HPLC data")
	
plot(fplpom.fit,fplpom.res,
	xlab = "Fitted Values (FPL)",
	ylab = "Standardized Residuals",
	main = "ELISA data")


###################################################
### chunk number 11: linCI
###################################################
ciu <- fitted(linmodel) + summary(linmodel)$sigma*qt(.975,linmodel$df)
cil <- fitted(linmodel) - summary(linmodel)$sigma*qt(.975,linmodel$df)


###################################################
### chunk number 12: HPLCcompareCI
###################################################
par(mfrow=c(1,2))
plot(HPLC, main = "HPLC data", sub="Without POM fit",col="blue",pch=16)
lines(conc.hplc,fitted(linmodel),col="lightblue")
lines(conc.hplc,ciu,col="lightblue",lty=2)
lines(conc.hplc,cil,col="grey",lty=2)

plot(cal.lin.pom,print=FALSE,main = "HPLC data", sub = "With POM fit",xlab = "Concentration", ylab = "Response")


###################################################
### chunk number 13: ELISAcompareCI
###################################################
par(mfrow=c(1,2))
plot(cal.fpl,print=FALSE,main = "ELISA data", sub = "Without POM fit",xlab = "Concentration", ylab = "Response")
plot(cal.fpl.pom,print=FALSE,main = "ELISA data", sub = "With POM fit",xlab = "Concentration", ylab = "Response")


###################################################
### chunk number 14: calibCalc
###################################################
calib.lin <- calib(cal.lin.pom,resp.hplc)
calib.fpl <- calib(cal.fpl.pom,resp.elisa)


###################################################
### chunk number 15: HPLCcalibPlot
###################################################
plot(calib.lin,main="HPLC data")


###################################################
### chunk number 16: ELISAcalibPlot
###################################################
plot(calib.fpl,main="ELISA data")


###################################################
### chunk number 17: loadLibrary
###################################################
library(calib)


###################################################
### chunk number 18: loadData
###################################################
data(HPLC)
data(ELISA)


###################################################
### chunk number 19: nameData
###################################################
conc.hplc <- HPLC[,1]
resp.hplc <- HPLC[,2]
conc.elisa <- ELISA[,1]
resp.elisa <- ELISA[,2]


###################################################
### chunk number 20: plottingData
###################################################
par(mfrow=c(1,2))
#par(mar = c(3.5,3.5,1.5,1.5))
plot(conc.hplc,resp.hplc,
	xlab = "Concentration (ng/ml)",
	ylab = "Response",
	main = "HPLC data")
linmodel <- lm(resp.hplc~conc.hplc)
# The predicted response
linPredResp <- fitted(linmodel)
# Linear regression fit
lines(conc.hplc,linPredResp)

## Plot of the data with a std FPL fit
#par(mar = c(3.5,3.5,1.5,1.5))
plot(log(conc.elisa),resp.elisa,
	xlab = "log(Concentration (ng/ml))",
	ylab = "Response",
	main = "ELISA data")
fplmodel <- calib.fit(conc.elisa,resp.elisa,type="log.fpl")
# The predicted response
fplPredResp <- fplmodel@fitted.values
# fpl regression fit
lines(log(conc.elisa),fplPredResp)


###################################################
### chunk number 21: plottingLinRes
###################################################
par(mfrow=c(1,2))
#par(mar = c(3.5,3.5,1.5,1.5))
## Residuals from linear fit
linres <- residuals(linmodel)/summary(linmodel)$sigma
plot(linPredResp,linres,
	 xlab = "Predicted Value of Mean (LS)",
	 ylab = "Standardized Residuals",
	 main = "HPLC data",
	 ylim = c(-5,5))
abline(h=0)

## Residuals from fpl fit
fplres <- fplmodel@residuals/fplmodel@sigma
#par(mar = c(3.5,3.5,1.5,1.5))
plot(fplPredResp,fplres,
	xlab = "Predicted Value of Mean (FPL)",
	ylab = "Standardized Residuals",
	main = "ELISA data",
	ylim = c(-5,5))
abline(h=0)


###################################################
### chunk number 22: plottingLogHPLC
###################################################
plot(conc.hplc,log(resp.hplc),
	 xlab = "Concentration (ng/ml)",
	 ylab = "log(Response)")


###################################################
### chunk number 23: plottingSqrtHPLC
###################################################
plot(conc.hplc,sqrt(resp.hplc),
	 xlab = "Concentration (ng/ml)",
	 ylab = "sqrt(Response)")


###################################################
### chunk number 24: plottingLogAbsLog
###################################################
par(mfrow=c(1,2))
#par(mar = c(3.5,3.5,1.5,1.5))
plot(log(linPredResp),log(abs(linres)),
	xlab = "Log(LS predicted values)",
	ylab = "Log(absolute LS residuals)",
	main = "HPLC data",
	ylim = c(-6,2))
linresmodel <- lm(log(abs(linres))~log(linPredResp))
lines(log(linPredResp),fitted(linresmodel))

#par(mar = c(3.5,3.5,1.5,1.5))
plot(log(fplPredResp),log(abs(fplres)),
	xlab = "Log(FPL predicted values)",
	ylab = "Log(absolute FPL residuals)",
	main = "ELISA data",
	ylim = c(-6,2))
fplresmodel <- lm(log(abs(fplres))~log(fplPredResp))
lines(log(fplPredResp),fitted(fplresmodel))


###################################################
### chunk number 25: calib_fitFuns
###################################################
cal.fpl <- calib.fit(conc.elisa,resp.elisa,type="log.fpl")
cal.lin.pom <- calib.fit(conc.hplc,resp.hplc,type="lin.pom")
cal.fpl.pom <- calib.fit(conc.elisa,resp.elisa,type="log.fpl.pom")

linpom.fit <- cal.lin.pom@fitted.values
fplpom.fit <- cal.fpl.pom@fitted.values

sig.lin <- cal.lin.pom@sigma
sig.fpl <- cal.fpl.pom@sigma

theta.lin <- cal.lin.pom@theta
theta.fpl <- cal.fpl.pom@theta

linpom.res <- cal.lin.pom@residuals*(1/((linpom.fit^theta.lin)*sig.lin))
fplpom.res <- cal.fpl.pom@residuals*(1/((fplpom.fit^theta.fpl)*sig.fpl))


###################################################
### chunk number 26: plottingVferesid
###################################################
par(mfrow=c(1,2))
#par(mar = c(3.5,3.5,1.5,1.5))
plot(linpom.fit,linpom.res,
	xlab = "Fitted Values (LS)",
	ylab = "Standardized Residuals",
	main = "HPLC data")

#par(mar = c(3.5,3.5,1.5,1.5))	
plot(fplpom.fit,fplpom.res,
	xlab = "Fitted Values (FPL)",
	ylab = "Standardized Residuals",
	main = "ELISA data")


###################################################
### chunk number 27: linCIFuns
###################################################
ciu <- fitted(linmodel) + summary(linmodel)$sigma*qt(.975,linmodel$df)
cil <- fitted(linmodel) - summary(linmodel)$sigma*qt(.975,linmodel$df)


###################################################
### chunk number 28: plottingHPLCcompareCI
###################################################
par(mfrow=c(1,2))
#par(mar = c(3.5,3.5,1.5,1.5))
plot(HPLC, main = "HPLC data", sub="Without POM fit",col="blue",pch=16)
lines(conc.hplc,fitted(linmodel),col="lightblue")
lines(conc.hplc,ciu,col="lightblue",lty=2)
lines(conc.hplc,cil,col="grey",lty=2)

#par(mar = c(3.5,3.5,1.5,1.5))
plot(cal.lin.pom,print=FALSE,main = "HPLC data", sub = "With POM fit",xlab = "Concentration", ylab = "Response")


###################################################
### chunk number 29: plottingELISAcompareCI
###################################################
par(mfrow=c(1,2))
#par(mar = c(3.5,3.5,1.5,1.5))
plot(cal.fpl,print=FALSE,main = "ELISA data", sub = "Without POM fit",xlab = "Concentration", ylab = "Response")

#par(mar = c(3.5,3.5,1.5,1.5))
plot(cal.fpl.pom,print=FALSE,main = "ELISA data", sub = "With POM fit",xlab = "Concentration", ylab = "Response")


###################################################
### chunk number 30: calibCalcFun
###################################################
calib.lin <- calib(cal.lin.pom,resp.hplc)
calib.fpl <- calib(cal.fpl.pom,resp.elisa)


###################################################
### chunk number 31: plottingHPLCcalibPlot
###################################################
plot(calib.lin,main="HPLC data")


###################################################
### chunk number 32: plottingELISAcalibPlot
###################################################
plot(calib.fpl,main="ELISA data")


