#############################################################
## This plots the output from calib.fpl.pom                ##
## modified by Matthew Mitchell on 3/27/03                 ##
##                                                         ##
## input:                                                  ##
##  calout: the output data frame from calib.fpl.pom       ##
##  int.meth: indicates which intervals are shown on       ##
##            the plot.  Originally the default was "i"    ##
##            for inverse intervals but this was changed   ##
##            by MM on 3/27 to "w" for the Wald intervals  ##
##            because when m is a vector with unequal      ##
##            elements, the inverse intervals are not      ##
##            computed corectly.                           ##
##  ylim, xlim, ylab, xlab, join give plotting info        ##
##  mdcline: default is T but it is changed in the body    ##
##           of this program when the dilutions are not    ##
##           all 1 or when the number of reps from the     ##
##           samples is different from those used for      ##
##           the standard curve                            ##
##  text.size, oma, cex are graphing parameters            ##
##                                                         ##
## the original comments are indicated with # or ###       ##
## and comments added by MM are indicated with a ##        ##
## modified on 9/1/04                                      ##
#############################################################


plot.calib <- function(calib.out, int.meth = "w", ylim, xlim = c(min(0,
    calib.out@labels$plot.ind), 1.1 * max(calib.out@labels$
    plot.ind)), ylab, xlab, join = calib.out@labels$
    cont.pl, mdcline = T, txt.size = 1, addcal.mat, ...,
    oma = c(2, 2, 2, 2), cex = 1, extrap=calib.out@extrap,
    main, las = 2, cex.axis = .7)
{

  ## Not sure about the use of this line
#calout$calib <- calout[ ,1]  ## MWM 9/1/04 (I renamed this variable)
  
#browser()
    old.par <- par("oma", "cex")
    on.exit(par(old.par))
## oma is the vector c(bottom, left, top, right).
## This gives the size of the outer margins in lines of text.
    par(oma = oma, cex = cex) #outer margin text lines and point size
    max.x <- calib.out@max.x * (max(1/calib.out@dilution,na.rm=TRUE))
    if(missing(ylim))
        ylim <- c(0, 1.05 * max.x)
    cal.plot <- calib.out@Estimated.x
    too.low <- calib.out@labels$too.low
    too.high <- calib.out@labels$too.high
if(!extrap){
    cal.plot <- ifelse(too.high == 1, max.x * 1.01, cal.plot)
    cal.plot <- ifelse(too.low == 1, 0, cal.plot)}

## line added by MM on 3/27/03
    mdcline <- ifelse(all(calib.out@dilution == 1) & calib.out@repeq == T,T,F) 
    
    if(int.meth == "i") {
        xlow <- calib.out@inver.low
        xup <- calib.out@inver.up
    }
    else {
        xlow <- calib.out@wald.low
        xup <- calib.out@wald.up
    }
    cal.low <- ifelse(is.finite(xlow), xlow, max.x)
    cal.up <- ifelse(!is.finite(xup), max(ylim), xup)
    ylim[2] <- 1.05*max(cal.up)
    x <- calib.out@labels$plot.ind
    x.ax <- calib.out@labels$xlabs
    if(missing(ylab)) {
        if(calib.out@labels$dose.units != "")
            ylab <- paste("Calibrated", calib.out@labels$dose.name,
                          "concentration (", 
                calib.out@labels$dose.units, 
                ")")
        else ylab <- paste("Calibrated", calib.out@labels$dose.name,
                           "concentration")
    }
    if(missing(xlab)) {
        if(calib.out@labels$xunits == "")
            xlab <- calib.out@labels$name
        else xlab <- paste(calib.out@labels$name, " (", 
                calib.out@labels$xunits, ")", 
                sep = "")
    }
    if(missing(main))
	    main <- "Plot of Calibrated Concentrations by Sample Number"
	else
		main <- main
    plot(x, cal.plot, type = "n", xlab = xlab, ylab = ylab, pch = 5,
        ylim = ylim, xlim = xlim, xaxt = "n", main = main, ...)  
    ## fix done in calib so that these points() make sense (moc 5/8/95)
    grid(nx=NA,ny=NULL)

## Since the oor of range are set to 0 or 1.01max.x we don't want
## confidence intervals to show up either
## These loope were added by MM on 3/27/03
    if(!extrap){
    cal.low <- ifelse(too.high == 1|too.low == 1,cal.plot,cal.low)
    cal.up <-  ifelse(too.high == 1|too.low == 1,cal.plot,cal.up)}
    
    segments(x, cal.low, x, cal.up,col="lightblue")

#    indColor <- ifelse(cal.low>attributes(calout)$mdc[1],T,F)
#    segments(x[indColor], cal.low[indColor], x[indColor], cal.up[indColor],
#            col="lightblue",lwd=1.5)

    points(x[too.low == 0 & too.high == 0], cal.plot[too.low == 0 & 
        too.high == 0], pch = 16, col="blue")
    points(x[too.low == 1 | too.high == 1], cal.plot[too.low == 1 | 
        too.high == 1], pch = 4)
    axis(1, at = x, labels = x.ax, cex = 1 * txt.size, las = las, cex.axis = cex.axis)
    if(join == T)
        lines(x[order(x)], cal.plot[order(x)])
    points(x[cal.up > max.x], cal.up[cal.up > max.x], pch = "^")
    abline(h = 0, lty = 1, col="gray")

    if(length(calib.out@truth)!=0) {
        truth <- calib.out@truth
        min.diff <- function(b)
        {
            len.b <- length(b)
            bmat <- cbind(b[1:(len.b - 1)], b[2:len.b])
            diffs <- bmat[, 1] - bmat[, 2]
            mdiff <- min(abs(diffs[diffs != 0]))
            mdiff
        }
        slen <- min.diff(x)/4
        segments(x - slen, truth, x + slen, truth)
    }
    if(mdcline) {
        mdcval <- calib.out@mdc[1]
        abline(h = mdcval, lty = 2, col="violetred1")
        mtext("MDC", side = 4, line=1, at = mdcval, las=1, col="violetred1")
    }
    if(sum(too.low + too.high,na.rm=TRUE) > 0)
        text(x[1], ylim[2], "x = response out of range", adj = 
            0, cex = txt.size)
    if(!missing(addcal.mat))
        points(addcal.mat[, 1], addpts.mat[, 2], pch = 0)
    invisible()
}

plot.fpl.pom <- function(fpl.out, conf = 0.9, m = fpl.out@m, main = "", xlab = "",
     ylab = "", xlim = c(0, max(fpl.out@x)), ylim = c(0, 1.1 * max(fpl.out@y, na.rm = T)),
     logplot = F, print.txt = T, click.txt = F, start.txt, increm, dig = c(4, 3), txt.size = 1,
     heading = T, lof.print = F, pred.lim = F, ...)
{
######################################################################
## This plots objects of class fpl.pom                              ##
## It produces a plot with the data points and the fitted curve.    ##
## The bands represent point-wise confidence intervals.             ##
## Summary statistics are added to the plot.                        ##
##                                                                  ##
## inputs variables                                                 ##
##  fpl.out: output data set from the fpl or fpl.pom function       ##
##  conf: confidence level for point-wise intervals for y           ##
##  m: number of reps (number of ys observed for a given x)         ##
##  main: title of the graph, no title is given by default          ##
##  xlab: the label for the x-axis, is blank by default             ##
##  ylab: the label for the y-axis, is blank by default             ##
##  xlim: sets the upper and lower limits of the x-axis             ##
##  ylim: sets the upper and lower limits of the y-axis             ##
##  logplot: determines of x is plotted on logscale, if             ##
##           parm=1 (see fpl, fpl.pom) then logplot=F, if           ##
##           parm=2 then logplot=T                                  ##
##  print.txt: is T if you want the stats printed on the graph      ##
##  click.txt: is T if you click to determine where it goes         ##
##             if parm=2 this is T                                  ##
##  start.txt: values are set by where you click                    ##
##  increm: helps determine the position of the text                ##
##  dig:                                                            ##
##  txt.size:                                                       ##
##  heading: if true if header is printed                           ##
##  lof.print: is T is lof is printed                               ##
##  pred.lim: F by default, conf limits set, T prediction limits    ##                                         
##                                                                  ##
## modified by Matthew Mitchell on 3/14/03                          ##
## comments added by Matthew Mitchell                               ##
## are indicated with ##, original comments with ###                ##
######################################################################
  
    x <- fpl.out@x
    y <- fpl.out@y
    if(fpl.out@parm == 2)
        logplot <- T
    
    if(print.txt) {
    if(missing(increm)) increm <- c(0.1 * max(x,    
###  increm <- c(0.1 * max(x, na.rm = T), 0.025 * max(y, na.rm = T))
### Min modified for fit to PC 
            na.rm = T), 0.04 * max(y, na.rm = T))   #  }
        if(missing(start.txt)) {
            if(fpl.out@coefficients[1] > fpl.out@coefficients[2])
                start.txt <- c(0.5 * max(x, na.rm = T), 
                  0.8 * max(y, na.rm = T))
            else start.txt <- c(0.5 * max(x, na.rm = T), 
                  0.5 * max(y, na.rm = T))
          }
    }
### if logplot=T, set up lower x-axis limit for logarithimic plot
    logind <- ifelse(logplot, "x", "")
    if(logplot) {
        xlim[1] <- ifelse(xlim[1] == 0, ifelse(min(x) > 0, min(x) - 0.0001,
         min(x[x != 0]) - 0.0001), xlim[1])
    if(xlim[1]<0) xlim[1]<-0.000001
## confusing when logplot is T click.txt was set equal to T so how can this be executed?
        if((print.txt & !click.txt)) {
            increm[1] <- log(increm[1])
        }
        on.exit(par(xaxt = "s"))
    }

############################
### plot the data points ###
############################
    plot(x, y, ylim = ylim, xlim = xlim, xlab = xlab, ylab = ylab, 
        main = main, type = "n", log = logind, ...)
    points(x, y, pch = 16, cex=.8, col="blue")   

### get the necessary information from fpl.out
    b <- fpl.out@coefficients
    n <- length(x)
    sigma <- fpl.out@sigma
    sg.inv <- fpl.out@cov.unscaled
    se.bhat <- fpl.out@se.coefficients

#####################################################################
### set up value of x's and calculate predicted values for y then ###
###  draw the smooth curve                                     ###
#####################################################################
    minmax <- c(min(b[1], b[2]), max(b[1], b[2]))
    yp1 <- c(minmax[1], 0.2 * minmax[2] + 0.8 * minmax[1], 0.4 * 
        minmax[2] + 0.6 * minmax[1], 0.6 * minmax[2] + 0.4 * 
        minmax[1], 0.8 * minmax[2] + 0.2 * minmax[1])
    if(logplot) {
        xp1 <- (fpl.inverse(fpl.out, yp1))[-1]
        if(min(x[x > 0]) < xp1[1])
            xp1 <- c(min(x[x > 0]), xp1)
        if(max(x) > 500) {
            xp <- c(log(exp(seq(xp1[1], min(xp1[2], 500), 
                len = 200))), seq(min(xp1[2], 500), max(
                xp1[3], 500), len = 50), seq(max(xp1[3],
                500), max(max(x), 500), len = 50))
        }
        else {
            xp <- c(log(exp(seq(xp1[1], xp1[2], len = 200))
                ), seq(xp1[2], xp1[3], len = 50), seq(
                xp1[3], max(x), len = 50))
        }
    }
    else {
        xp1 <- fpl.inverse(fpl.out, yp1)
        xp <- c(seq(0, xp1[3], len = 150), seq(xp1[3], max(x), 
            len = 150))
    }
    yp <- fpl.model(xp, b, parm = fpl.out@parm)
    lines(xp, as.vector(yp), col="lightblue")    

#######################################################################
### calculate the confidence limits based on the GLS-FPL fit        ###
### then draw the lines on the graph.  if pred.lim=T, calculate     ###
### prediction limits, otherwise calculate calibration limits.      ###
#######################################################################
    qn.term2 <- diag(attributes(yp)$gradient %*% 
    			sg.inv %*% t(attributes(yp)$gradient))
    qn.term1.c <- (as.vector(yp)^2)^fpl.out@theta/m
    qn.term1.p <- (as.vector(yp)^2)^fpl.out@theta
    qn.term2 <- as.vector(qn.term2)
    if(pred.lim)
        qn <- sigma * sqrt(qn.term1.p + qn.term2)
    else qn <- sigma * sqrt(qn.term1.c + qn.term2)
    n <- length(x)
    tcrit <- qt(1 - (1 - conf)/2, fpl.out@df.residual)
    ucl <- as.vector(yp) + tcrit * qn
    lcl <- as.vector(yp) - tcrit * qn
    lines(xp, ucl, lty = 2, col="lightblue") #col="gray")
    lines(xp, lcl, lty = 2, col="gray") 
    grid()


par(new=T) #, xaxs="d")
plot(0:1, 0:1, xlab="", ylab="", type="n", axes=F)
if(b[2] - b[1]){start.txt <- c(0,.9)}
else{start.txt <- c(.5,.9)}

increm <- c(.15,.0375)
#############################################
### text for the model summary statistics ###
#############################################
    if(print.txt) {

### print fpl equation
        if(heading) {
            if(fpl.out@parm == 1) {
                if(fpl.out@vfe.method == "OLS")
                  text(start.txt[1], start.txt[2] + 1.2 * increm[2], 
                    #"FPL-OLS: y=(b1-b2)/(1+(x/b3)^b4)+b2 + \n sigma*error",
                    expression(FPL[OLS]: y==frac(beta[1]-beta[2],1+(frac(x,beta[3]))^beta[4])
                        +beta[2]+sigma*epsilon), adj = 0, cex = txt.size)
                else text(start.txt[1], start.txt[2] + 1.2 * increm[2], 
###                    "FPL-POM: y=(b1-b2)/(1+(x/b3)^b4)+b2 + \n (mu^theta)*sigma*error",
## typo here              expression(FPL[OLS]: y==frac(beta[1]-beta[2],1+(frac(x,beta[3]))^beta[4])
## changed to below by MM 3/6/03  +beta[2]+mu^theta*sigma*epsilon), adj = 0, cex = txt.size)
                         expression(FPL[POM]: y==frac(beta[1]-beta[2],1+(frac(x,beta[3]))^beta[4])
                        +beta[2]+mu^theta*sigma*epsilon), adj = 0, cex = txt.size)
                          
            }
            else {
                if(fpl.out@vfe.method == "OLS")
                  text(start.txt[1], start.txt[2] + 1.2 * increm[2], 
###                    "FPL-OLS: y=(b1-b2)/(1+exp(b4(log(x)-b3)))+b2 + \n sigma*error",
## wrong function              expression(FPL[OLS]: y==frac(beta[1]-beta[2],1+(frac(log(x),beta[3]))^beta[4])
## changed to below by MWM 3/6/03   +beta[2]+sigma*epsilon), adj = 0, cex = txt.size)
                        expression(FPL[OLS]: y==frac(beta[1]-beta[2],1+ exp(beta[4]*(log(x)- beta[3])))
                        +beta[2]+ sigma*epsilon), adj = 0, cex = txt.size)
                else text(start.txt[1], start.txt[2] + 1.2 * increm[2], 
###                    "FPL-POM: y=(b1-b2)/(1+exp(b4(log(x)-b3)))+b2 + \n (mu^theta)*sigma*error",
## wrong function          expression(FPL[OLS]: y==frac(beta[1]-beta[2],1+(frac(log(x),beta[3]))^beta[4])
## changed to below by MM 3/6/03              +beta[2]+mu^theta*sigma*epsilon), adj = 0, cex = txt.size)
                          expression(FPL[POM]: y==frac(beta[1]-beta[2],1+ exp(beta[4]*(log(x)- beta[3])))
                        +beta[2]+mu^theta*sigma*epsilon), adj = 0, cex = txt.size)
                          
  
            }
        }
click.txt <- F
## how can the next line execute? - we just set it equal to F
        if(click.txt) {
            increm <- array(c(0, 0))
            cat("\n Click Mouse for beginning of first line of table \n")
            start.txt2 <- locator(1)
            increm[2] <- 0.6 * (start.txt1$y - start.txt2$y)
            col1 <- start.txt2$x
        }
        else col1 <- start.txt[1] + 0.1 * increm[1]
        text(col1, start.txt[2] - increm[2], "Parameter", adj = 0, cex = 0.8 * txt.size)
        if(click.txt) {
            cat("\n now click for placement of second column \n")
            col2 <- locator(1)$x
        }
        else col2 <- col1 + 1.1 * increm[1]
        text(col2, start.txt[2] - increm[2], "Estimate", adj = 0, cex = 0.8 * txt.size)
        if(click.txt) {
            cat("\n now click for placement of third column \n")
            col3 <- locator(1)$x
        }
        else col3 <- col2 + 1.1 * increm[1]
        text(col3, start.txt[2] - increm[2], "Std. Error", adj = 0, cex = 0.8 * txt.size)
        if(length(dig) == 1)
            dig <- rep(dig, 2)  
### print parameter estimates and standard errors
        
        text(col1, start.txt[2] - 2 * increm[2], expression(beta[1]), adj = 0, cex = 0.8 * txt.size)
        text(col2, start.txt[2] - 2 * increm[2],
            paste(format(round(b[1], dig[1]))), adj = 0, cex = 0.8 * txt.size)
        text(col3, start.txt[2] - 2 * increm[2],
            paste(format(round(se.bhat[1], dig[1]))), adj = 0, cex = 0.8 * txt.size)
        text(col1, start.txt[2] - 3 * increm[2], expression(beta[2]), adj = 0, cex = 0.8 * txt.size)
        text(col2, start.txt[2] - 3 * increm[2],
            paste(format(round(b[2], dig[1]))), adj = 0, cex = 0.8 * txt.size)
        text(col3, start.txt[2] - 3 * increm[2],
            paste(format(round(se.bhat[2], dig[1]))), adj = 0, cex = 0.8 * txt.size)
        if(fpl.out@parm == 1)
            text(col1, start.txt[2] - 4 * increm[2], expression(beta[3]==B50), adj = 0, cex = 0.8 * txt.size)
        else text(col1, start.txt[2] - 4 * increm[2], 
                expression(beta[3]==log(B50)), adj = 0, cex = 0.8 * txt.size)
        text(col2, start.txt[2] - 4 * increm[2],
            paste(format(round(b[3], dig[1]))), adj = 0, cex = 0.8 * txt.size)
        text(col3, start.txt[2] - 4 * increm[2],
            paste(format(round(se.bhat[3], dig[1]))), adj = 0, cex = 0.8 * txt.size)
        text(col1, start.txt[2] - 5 * increm[2], expression(beta[4]), adj = 0, cex = 0.8 * txt.size)
        text(col2, start.txt[2] - 5 * increm[2], 
            paste(format(round(b[4], dig[1]))), adj = 0, cex = 0.8 * txt.size)
        text(col3, start.txt[2] - 5 * increm[2], paste(format(
            round(se.bhat[4], dig[1]))), adj = 0, cex = 0.8 * txt.size)

### print mdc and rdl
        if(!is.null(fpl.out@mdc)) {
            text(col1, start.txt[2] - 6.5 * increm[2], 
                paste("MDC (m=", fpl.out@m, ") = ", format(round(fpl.out@mdc, dig[1])),
                 sep = ""), adj = 0, cex = 0.7 * txt.size)
            text((col2 + 1.1 * col3)/2, start.txt[2] - 6.5 * increm[2],
                paste("sigma = ", format(round(sigma, dig[1])),
                 sep = ""), adj = 0, cex = 0.7 * txt.size)
            text(col1, start.txt[2] - 7.5 * increm[2], 
                paste("RDL (m=", fpl.out@m, ") = ", format(round(fpl.out@rdl, dig[1])),
                 sep = ""), adj = 0, cex = 0.7 * txt.size)
            if(is.numeric(fpl.out@loq)) {
                text(col1, start.txt[2] - 8.5 * increm[2],
                 paste("LOQ (m=", fpl.out@m, ", CV=", fpl.out@cv, ") = ",
                  format(round(fpl.out@loq, dig[1])), sep = ""),
                  adj = 0, cex = 0.7 * txt.size)
            }
            else {
                text(col1, start.txt[2] - 8.5 * increm[2],
                 paste("LOQ (m=", fpl.out@m, ", CV=", fpl.out@cv, ") = NA",
                  sep = ""), adj = 0, cex = 0.7 * txt.size)
            }
            last.line <- 8.5
        }
        else last.line <- 5
        if(lof.print) {
            if(length(x) > length(unique(x))) {
                if(is.null(fpl.out@lof.test
                  ))
                  lof.out <- lof.test(fpl.out)
                else lof.out <- fpl.out@lof.test
                text(col1, start.txt[2] - 10 * increm[2
                  ], "LOF: Source", adj = 0, cex = 0.75 *
                  txt.size)
                text(col2, start.txt[2] - 10 * increm[2
                  ], "SS", adj = 0, cex = 0.75 * 
                  txt.size)
                text((col2 + col3)/2, start.txt[2] - 10 *
                  increm[2], "df", adj = 0, cex = 0.75 * 
                  txt.size)
                text(col3, start.txt[2] - 10 * increm[2
                  ], "MS", adj = 0, cex = 0.75 * 
                  txt.size)
                text(col1, start.txt[2] - 11 * increm[2
                  ], "Error", adj = 0, cex = 0.7 * 
                  txt.size)
                text(col1, start.txt[2] - 12 * increm[2
                  ], "Pure Error", adj = 0, cex = 0.7 * 
                  txt.size)
                text(col1, start.txt[2] - 13 * increm[2
                  ], "Lack of Fit", adj = 0, cex = 0.7 * 
                  txt.size)
                text(col2, start.txt[2] - 11 * increm[2
                  ], paste(format(round(lof.out@sse, 
                  dig[2]))), adj = 0, cex = 0.7 * 
                  txt.size)
                text(col2, start.txt[2] - 12 * increm[2
                  ], paste(format(round(lof.out@
                  pure.error, dig[2]))), adj = 0, cex
                   = 0.7 * txt.size)
                text(col2, start.txt[2] - 13 * increm[2
                  ], paste(format(round(lof.out@lofss, 
                  dig[2]))), adj = 0, cex = 0.7 * 
                  txt.size)
                text((col2 + col3)/2, start.txt[2] - 11 *
                  increm[2], paste(lof.out@df.sse), adj
                   = 0, cex = 0.7 * txt.size)
                text((col2 + col3)/2, start.txt[2] - 12 *
                  increm[2], paste(lof.out@
                  df.pure.error), adj = 0, cex = 0.7 * 
                  txt.size)
                text((col2 + col3)/2, start.txt[2] - 13 *
                  increm[2], paste(lof.out@df.lof), adj
                   = 0, cex = 0.7 * txt.size)
                text(col3, start.txt[2] - 11 * increm[2
                  ], paste(format(round(lof.out@sse/
                  lof.out@df.sse, dig[2]))), adj = 0, 
                  cex = 0.7 * txt.size)
                text(col3, start.txt[2] - 12 * increm[2
                  ], paste(format(round(lof.out@
                  pure.error/lof.out@df.pure.error, dig[
                  2]))), adj = 0, cex = 0.7 * txt.size)
                text(col3, start.txt[2] - 13 * increm[2
                  ], paste(format(round(lof.out@lofss/
                  lof.out@df.lof, dig[2]))), adj = 0, 
                  cex = 0.7 * txt.size)
                text(col1, start.txt[2] - 14.5 * increm[
                  2], paste("F =", format(round(lof.out@
                  Fstat, 2))), adj = 0, cex = 0.75 * 
                  txt.size)
                text((col2 + col3)/2, start.txt[2] - 
                  14.5 * increm[2], paste("P-value =", 
                  format(round(lof.out@p.value, 3))), 
                  adj = 0, cex = 0.75 * txt.size)
                last.line <- 14.5
            }
        }
### print POM variance function parameter if relevant
        if(fpl.out@vfe.method != "OLS")
            text(col1, start.txt[2] - (last.line + 1.5) * 
                increm[2], paste("POM Theta =", format(
                round(fpl.out@theta, 2))), adj = 0, cex
                 = 0.75 * txt.size)
    }
### TDS, 990105; the next lines copies the completed graph to the clipboard
### for pasting into Excel:
###    guiExecuteBuiltIn("$$Standard$Copy", wait = T)
    invisible()
}

plot.thpl.pom <- function(thpl.out, conf = 0.9, m = thpl.out@m, main = "", xlab = "",
     ylab = "", xlim = c(0, max(thpl.out@x)), ylim = c(0, 1.1 * max(thpl.out@y, na.rm = T)),
     logplot = F, print.txt = T, 
	 click.txt = F,  
	 start.txt, increm, dig = c(4, 3), txt.size = 1,
     heading = T, lof.print = F, bios=F, pred.lim = F, ...)
{
######################################################################
## This plots objects of class thpl.pom                             ##
## It produces a plot with the data points and the fitted curve.    ##
## The bands represent point-wise confidence intervals.             ##
## Summary statistics are added to the plot.                        ##
##                                                                  ##
## inputs variables                                                 ##
##  thpl.out: output data set from the thpl or thpl.pom function    ##
##  conf: confidence level for point-wise intervals for y           ##
##  m: number of reps (number of ys observed for a given x)         ##
##  main: title of the graph, no title is given by default          ##
##  xlab: the label for the x-axis, is blank by default             ##
##  ylab: the label for the y-axis, is blank by default             ##
##  xlim: sets the upper and lower limits of the x-axis             ##
##  ylim: sets the upper and lower limits of the y-axis             ##
##  logplot: determines of x is plotted on logscale                 ##
##  print.txt: is T if you want the stats printed on the graph      ##
##  click.txt: is T if you click to determine where it goes         ##
##  start.txt: values are set by where you click                    ##
##  increm: helps determine the position of the text                ##
##  dig:                                                            ##
##  txt.size:                                                       ##
##  heading: if true if header is printed                           ##
##  lof.print: is T is lof is printed                               ##
##  pred.lim: F by default, conf limits set, T prediction limits    ##
##  bios: if T prints model with f, glucose, F_inf, F_0, K_d        ##
##                                                                  ##
## written by Matthew Mitchell on 5/20/03                           ##
## comments added by Matthew Mitchell                               ##
## are indicated with ##, original comments with ###                ##
######################################################################
  
    x <- thpl.out@x
    y <- thpl.out@y
  
    if(print.txt) {
    if(missing(increm)) increm <- c(0.1 * max(x,na.rm = T), 0.04 * max(y, na.rm = T))  
        if(missing(start.txt)) {
            if(thpl.out@coefficients[1] > thpl.out@coefficients[2])
                start.txt <- c(0.5 * max(x, na.rm = T), 
                  0.8 * max(y, na.rm = T))
            else start.txt <- c(0.5 * max(x, na.rm = T), 
                  0.5 * max(y, na.rm = T))
          }
    }

    logind <- ifelse(logplot, "x", "")
    if(logplot) {
        xlim[1] <- ifelse(xlim[1] == 0, ifelse(min(x) > 0, min(x) - 0.0001,
         min(x[x != 0]) - 0.0001), xlim[1])

        if(print.txt & !click.txt) {
            increm[1] <- log(increm[1])
        }
        on.exit(par(xaxt = "s"))
    }

############################
### plot the data points ###
############################
    plot(x, y, ylim = ylim, xlim = xlim, xlab = xlab, ylab = ylab, 
        main = main, type = "n", log = logind, ...)
    points(x, y, pch = 16, cex=.8, col="blue")   

### get the necessary information from thpl.out
    b <- thpl.out@coefficients
    n <- length(x)
    sigma <- thpl.out@sigma
    sg.inv <- thpl.out@cov.unscaled
    se.bhat <- thpl.out@se.coefficients

#####################################################################
### set up value of x's and calculate predicted values for y then ###
###  draw the smooth curve                                     ###
#####################################################################
    minmax <- c(min(b[1], b[2]), max(b[1], b[2]))
    yp1 <- c(minmax[1], 0.2 * minmax[2] + 0.8 * minmax[1], 0.4 * 
        minmax[2] + 0.6 * minmax[1], 0.6 * minmax[2] + 0.4 * 
        minmax[1], 0.8 * minmax[2] + 0.2 * minmax[1])
    if(logplot) {
        xp1 <- (thpl.inverse(thpl.out, yp1))[-1]
        if(min(x[x > 0]) < xp1[1])
            xp1 <- c(min(x[x > 0]), xp1)
        if(max(x) > 500) {
            xp <- c(log(exp(seq(xp1[1], min(xp1[2], 500), 
                len = 200))), seq(min(xp1[2], 500), max(
                xp1[3], 500), len = 50), seq(max(xp1[3],
                500), max(max(x), 500), len = 50))
        }
        else {
            xp <- c(log(exp(seq(xp1[1], xp1[2], len = 200))
                ), seq(xp1[2], xp1[3], len = 50), seq(
                xp1[3], max(x), len = 50))
        }
    }
    else {
        xp1 <- thpl.inverse(thpl.out, yp1)
        xp <- c(seq(0, xp1[3], len = 150), seq(xp1[3], max(x), 
            len = 150))
    }
    yp <- thpl.model(xp, b)
    lines(xp, as.vector(yp), col="lightblue")    

#######################################################################
### calculate the confidence limits based on the GLS-THPL fit        ###
### then draw the lines on the graph.  if pred.lim=T, calculate     ###
### prediction limits, otherwise calculate calibration limits.      ###
#######################################################################
    qn.term2 <- diag(attributes(yp)$gradient %*% sg.inv %*% t(
        attributes(yp)$gradient))
    qn.term1.c <- as.vector(yp)^(2 * thpl.out@theta)/m
    qn.term1.p <- as.vector(yp)^(2 * thpl.out@theta)
    qn.term2 <- as.vector(qn.term2)
    if(pred.lim)
        qn <- sigma * sqrt(qn.term1.p + qn.term2)
    else qn <- sigma * sqrt(qn.term1.c + qn.term2)
    n <- length(x)
    tcrit <- qt(1 - (1 - conf)/2, thpl.out@df.residual)
    ucl <- as.vector(yp) + tcrit * qn
    lcl <- as.vector(yp) - tcrit * qn
    lines(xp, ucl, lty = 2, col="lightblue") #col="gray")
    lines(xp, lcl, lty = 2, col="gray") 
    grid()


par(new=T) #, xaxs="d")
plot(0:1, 0:1, xlab="", ylab="", type="n", axes=F)
if(b[2] - b[1] > 0){start.txt <- c(0,.9)}
else{start.txt <- c(.5,.9)}
increm <- c(.15,.0375)
#############################################
### text for the model summary statistics ###
#############################################
    if(print.txt) {
		## Does not seem necessary anymore DVS
#        if(click.txt) {
#            start.txt <- array(c(0, 0))
#           ## This doesn't display until after you click!  MM 5/22/03
#           ## cat("\n Click Mouse on where you would like text to start \n"   
#           ##     )
#            start.txt1 <- locator(1)
#            start.txt[1] <- start.txt1$x
#            start.txt[2] <- start.txt1$y
#        }

### print thpl equation
        if(heading) {
                if(!bios){
                if(thpl.out@vfe.method == "OLS")
                  text(start.txt[1], start.txt[2] + 1.2 * increm[2], 
                        expression(THPL[OLS]: y==frac(beta[1]-beta[2],1+frac(x,beta[3]))
                        +beta[2]+sigma*epsilon), adj = 0, cex = txt.size)
                else text(start.txt[1], start.txt[2] + 1.2 * increm[2], 
                         expression(THPL[POM]: y==frac(beta[1]-beta[2],1+frac(x,beta[3]))
                       +beta[2]+mu^theta*sigma*epsilon), adj = 0, cex = txt.size)}
                if(bios){
                if(thpl.out@vfe.method == "OLS")
                  text(start.txt[1], start.txt[2] + 1.2 * increm[2], 
                        expression(THPL[OLS]: f==frac(F[0]-F[inf],1+frac(glc,K[d]))
                        +F[inf]+sigma*epsilon), adj = 0, cex = txt.size)
                else text(start.txt[1], start.txt[2] + 1.2 * increm[2], 
                         expression(THPL[POM]: f==frac(F[0]-F[inf],1+frac(glc,K[d]))
                       +F[inf]+mu^theta*sigma*epsilon), adj = 0, cex = txt.size)}
              }
        
click.txt <- F

        if(click.txt) {
            increm <- array(c(0, 0))
            cat("\n Click Mouse for beginning of first line of table \n")
            start.txt2 <- locator(1)
            increm[2] <- 0.6 * (start.txt1$y - start.txt2$y)
            col1 <- start.txt2$x
        }
        else col1 <- start.txt[1] + 0.1 * increm[1]
        text(col1, start.txt[2] - increm[2], "Parameter", adj = 0, cex = 0.8 * txt.size)
        if(click.txt) {
            cat("\n now click for placement of second column \n")
            col2 <- locator(1)$x
        }
        else col2 <- col1 + 1.1 * increm[1]
        text(col2, start.txt[2] - increm[2], "Estimate", adj = 0, cex = 0.8 * txt.size)
        if(click.txt) {
            cat("\n now click for placement of third column \n")
            col3 <- locator(1)$x
        }
        else col3 <- col2 + 1.1 * increm[1]
        text(col3, start.txt[2] - increm[2], "Std. Error", adj = 0, cex = 0.8 * txt.size)
        if(length(dig) == 1)
            dig <- rep(dig, 2)  
### print parameter estimates and standard errors
        if(!bios){
        text(col1, start.txt[2] - 2 * increm[2], expression(beta[1]), adj = 0, cex = 0.8 * txt.size)
        text(col2, start.txt[2] - 2 * increm[2],
            paste(format(round(b[1], dig[1]))), adj = 0, cex = 0.8 * txt.size)
        text(col3, start.txt[2] - 2 * increm[2],
            paste(format(round(se.bhat[1], dig[1]))), adj = 0, cex = 0.8 * txt.size)
        text(col1, start.txt[2] - 3 * increm[2], expression(beta[2]), adj = 0, cex = 0.8 * txt.size)
        text(col2, start.txt[2] - 3 * increm[2],
            paste(format(round(b[2], dig[1]))), adj = 0, cex = 0.8 * txt.size)
        text(col3, start.txt[2] - 3 * increm[2],
            paste(format(round(se.bhat[2], dig[1]))), adj = 0, cex = 0.8 * txt.size)
        text(col1, start.txt[2] - 4 * increm[2], expression(beta[3]==B50), adj = 0, cex = 0.8 * txt.size)
        text(col2, start.txt[2] - 4 * increm[2],
            paste(format(round(b[3], dig[1]))), adj = 0, cex = 0.8 * txt.size)
        text(col3, start.txt[2] - 4 * increm[2],
            paste(format(round(se.bhat[3], dig[1]))), adj = 0, cex = 0.8 * txt.size)}

        if(bios){
        text(col1, start.txt[2] - 2 * increm[2], expression(F[0]), adj = 0, cex = 0.8 * txt.size)
        text(col2, start.txt[2] - 2 * increm[2],
            paste(format(round(b[1], dig[1]))), adj = 0, cex = 0.8 * txt.size)
        text(col3, start.txt[2] - 2 * increm[2],
            paste(format(round(se.bhat[1], dig[1]))), adj = 0, cex = 0.8 * txt.size)
        text(col1, start.txt[2] - 3 * increm[2], expression(F[inf]), adj = 0, cex = 0.8 * txt.size)
        text(col2, start.txt[2] - 3 * increm[2],
            paste(format(round(b[2], dig[1]))), adj = 0, cex = 0.8 * txt.size)
        text(col3, start.txt[2] - 3 * increm[2],
            paste(format(round(se.bhat[2], dig[1]))), adj = 0, cex = 0.8 * txt.size)
        text(col1, start.txt[2] - 4 * increm[2], expression(K[d]), adj = 0, cex = 0.8 * txt.size)
        text(col2, start.txt[2] - 4 * increm[2],
            paste(format(round(b[3], dig[1]))), adj = 0, cex = 0.8 * txt.size)
        text(col3, start.txt[2] - 4 * increm[2],
            paste(format(round(se.bhat[3], dig[1]))), adj = 0, cex = 0.8 * txt.size)}
       

### print mdc and rdl
        if(!is.null(thpl.out@mdc)) {
            text(col1, start.txt[2] - 6.5 * increm[2], 
                paste("MDC (m=", thpl.out@m, ") = ", format(round(thpl.out@mdc, dig[1])),
                 sep = ""), adj = 0, cex = 0.7 * txt.size)
            text((col2 + 1.1 * col3)/2, start.txt[2] - 6.5 * increm[2],
                paste("sigma = ", format(round(sigma, dig[1])),
                 sep = ""), adj = 0, cex = 0.7 * txt.size)
            text(col1, start.txt[2] - 7.5 * increm[2], 
                paste("RDL (m=", thpl.out@m, ") = ", format(round(thpl.out@rdl, dig[1])),
                 sep = ""), adj = 0, cex = 0.7 * txt.size)
            if(is.numeric(thpl.out@loq)) {
                text(col1, start.txt[2] - 8.5 * increm[2],
                 paste("LOQ (m=", thpl.out@m, ", CV=", thpl.out@cv, ") = ",
                  format(round(thpl.out@loq, dig[1])), sep = ""),
                  adj = 0, cex = 0.7 * txt.size)
            }
            else {
                text(col1, start.txt[2] - 8.5 * increm[2],
                 paste("LOQ (m=", thpl.out@m, ", CV=", thpl.out@cv, ") = NA",
                  sep = ""), adj = 0, cex = 0.7 * txt.size)
            }
            last.line <- 8.5
        }
        else last.line <- 5
        if(lof.print) {
            if(length(x) > length(unique(x))) {
                if(length(thpl.out@lof.test@Fstat)==0)
                  lof.out <- lof.test(thpl.out)
                else lof.out <- thpl.out@lof.test
                text(col1, start.txt[2] - 10 * increm[2
                  ], "LOF: Source", adj = 0, cex = 0.75 *
                  txt.size)
                text(col2, start.txt[2] - 10 * increm[2
                  ], "SS", adj = 0, cex = 0.75 * 
                  txt.size)
                text((0.35*col2 + 0.65*col3), start.txt[2] - 10 *
                  increm[2], "df", adj = 0, cex = 0.75 * 
                  txt.size)
                text(col3, start.txt[2] - 10 * increm[2
                  ], "MS", adj = 0, cex = 0.75 * 
                  txt.size)
                text(col1, start.txt[2] - 11 * increm[2
                  ], "Error", adj = 0, cex = 0.7 * 
                  txt.size)
                text(col1, start.txt[2] - 12 * increm[2
                  ], "Pure Error", adj = 0, cex = 0.7 * 
                  txt.size)
                text(col1, start.txt[2] - 13 * increm[2
                  ], "Lack of Fit", adj = 0, cex = 0.7 * 
                  txt.size)
                text(col2, start.txt[2] - 11 * increm[2
                  ], paste(format(round(lof.out@sse, 
                  dig[2]))), adj = 0, cex = 0.7 * 
                  txt.size)
                text(col2, start.txt[2] - 12 * increm[2
                  ], paste(format(round(lof.out@
                  pure.error, dig[2]))), adj = 0, cex
                   = 0.7 * txt.size)
                text(col2, start.txt[2] - 13 * increm[2
                  ], paste(format(round(lof.out@lofss, 
                  dig[2]))), adj = 0, cex = 0.7 * 
                  txt.size)
                text((0.35*col2 + 0.65*col3), start.txt[2] - 11 *
                  increm[2], paste(lof.out@df.sse), adj
                   = 0, cex = 0.7 * txt.size)
                text((0.35*col2 + 0.65*col3), start.txt[2] - 12 *
                  increm[2], paste(lof.out@
                  df.pure.error), adj = 0, cex = 0.7 * 
                  txt.size)
                text((0.35*col2 + 0.65*col3), start.txt[2] - 13 *
                  increm[2], paste(lof.out@df.lof), adj
                   = 0, cex = 0.7 * txt.size)
                text(col3, start.txt[2] - 11 * increm[2
                  ], paste(format(round(lof.out@sse/
                  lof.out@df.sse, dig[2]))), adj = 0, 
                  cex = 0.7 * txt.size)
                text(col3, start.txt[2] - 12 * increm[2
                  ], paste(format(round(lof.out@
                  pure.error/lof.out@df.pure.error, dig[
                  2]))), adj = 0, cex = 0.7 * txt.size)
                text(col3, start.txt[2] - 13 * increm[2
                  ], paste(format(round(lof.out@lofss/
                  lof.out@df.lof, dig[2]))), adj = 0, 
                  cex = 0.7 * txt.size)
                text(col1, start.txt[2] - 14.5 * increm[
                  2], paste("F =", format(round(lof.out@
                  Fstat, 2))), adj = 0, cex = 0.75 * 
                  txt.size)
                text((col2 + col3)/2, start.txt[2] - 
                  14.5 * increm[2], paste("P-value =", 
                  format(round(lof.out@p.value, 3))), 
                  adj = 0, cex = 0.75 * txt.size)
                last.line <- 14.5
            }
        }
### print POM variance function parameter if relevant
        if(thpl.out@vfe.method != "OLS")
            text(col1, start.txt[2] - (last.line + 1.5) * 
                increm[2], paste("POM Theta =", format(
                round(thpl.out@theta, 2))), adj = 0, cex
                 = 0.75 * txt.size)
    }
    invisible()
}

#####################################################
## This plots output from the linear POM model     ##
##                                                 ##
## MM 9/26/03 made major changes                   ##
#####################################################

plot.lin.pom <- function(lin.out, conf = 0.9, m = lin.out@m, main = "", xlab = "", ylab
	 = "", xlim = c(0, max(lin.out@x)), ylim = c(0, 1.1 * max(
	lin.out@y, na.rm = T)), print.txt = T, start.txt, increm = 0.03,
	digits = c(2, 2, 4, 2, 2), txt.size = 0.8, lof.print = F, logplot=F,
	pred.lim = F, sub = "",...)  
{

## digits refers to the number of digits reported on the graph for paramter estimates
## digits[1]: digits of intercept
## digits[2]: digits of slope (and quadratic term if there is one)
## digits[3]: digits for MDC, RDL, LOQ
## digits[4]: digits used for reporting theta
## digits[5]: digits used for reporting sigma
## previously default was c(1,1,1,2,2) but this often gave 0s so I changed it

  
	x <- lin.out@x
	y <- lin.out@y	# assay performance characteristics
	mdc <- lin.out@mdc
	rdl <- lin.out@rdl
	loqcv <- lin.out@cv
	loq <- lin.out@loq
	theta <- lin.out@theta
	df.den <- lin.out@df.residual	# calc pred limits
	n <- length(x)
	xmax <- max(lin.out@x)
	xp <- seq(0, xmax, length = 100)
	if(lin.out@mmod == "lin") {
		yp <- lin.out@coefficients[1] + xp * lin.out@coefficients[2]
		xmat <- cbind(1, xp)}
	
	if(lin.out@mmod == "quad") {
		yp <- lin.out@coefficients[1] + xp * lin.out@coefficients[2] + xp * xp * 
			lin.out@coefficients[3]
		xmat <- cbind(1, xp, xp * xp)
	}
	if(pred.lim)
		mq <- 1
	else mq <- m
	sg.inv <- lin.out@cov.unscaled
	conf <- 1 - (1 - conf)/2
	tcrit <- qt(conf, df.den)
	pred.limit <- function(pom, xp, m, tcrit, upper = T)
	{
		if(missing(m))
			m <- pom@m
		if(attr(pom, "mmod") == "lin") {
			yp <- pom@coefficients[1] + xp * pom@coefficients[2]
			xmat <- cbind(1, xp)}
	
		if(pom@mmod == "quad") {
			yp <- pom@coefficients[1] + xp * pom@coefficients[2] + xp * xp * 
				pom@coefficients[3]
			xmat <- cbind(1, xp, xp * xp)
		}
		qn.term1 <- yp^(2 * pom@theta)/m
		qn.term2 <- xmat %*% pom@cov.unscaled %*% t(xmat)
		qn.term2 <- diag(qn.term2)
		qn.p <- pom@sigma * sqrt(qn.term1 + qn.term2)
		if(upper)
			cl <- yp + tcrit * qn.p
		else cl <- yp - tcrit * qn.p
		cl
	}
	ucl <- pred.limit(lin.out, xp, m = mq, tcrit, upper = T)
	lcl <- pred.limit(lin.out, xp, m = mq, tcrit, upper = F)	
	## plot cal curve
	y <- y[x <= xmax]
	x <- x[x <= xmax]
	ymax <- max(y, na.rm = T)
	if(missing(ylim))
		ylim <- c(0, 1.1 * max(y, na.rm = T))
	plot(x, y, xlim = xlim, ylim = ylim, main = main, xlab = xlab, 
		ylab = ylab, type = "n", sub = sub)#, ...) This addtional arg giving an error
        ## leaving it out unless they serve some other function DVS 10/23/06
	points(x, y, pch = 16, col="blue")	## 9/26 changed pch from 5 to 16 and made color blue
	## draw regression line and calibration limits
	if(lin.out@mmod == "lin")
		abline(lin.out@coefficients[1], lin.out@coefficients[2], col="lightblue") ## MM 9/26 changed from black to blue
	if(lin.out@mmod == "quad")
		lines(xp[yp <= max(ylim)], yp[yp <= max(ylim)], col="lightblue") ## changed from black 
	lines(xp[ucl <= max(ylim)], ucl[ucl <= max(ylim)], lty = 2, col="lightblue") ## changed from black
	lines(xp[lcl <= max(ylim)], lcl[lcl <= max(ylim)], lty = 2, col="grey") ## changed from black
	if(print.txt) {
		if(missing(start.txt)) {
			if(lin.out@coefficients[2] > 0)
				start.txt <- c(0, max(ylim))
			else start.txt <- c(mean(xlim), max(ylim))  ## added 11/21 MM
		}
		ind <- start.txt[2] - diff(ylim) * increm * seq(0, 14)
## MM: I completely revamped the printout on the plot 9/26/03
       if(lin.out@mmod == "lin" & lin.out@var.model=="OLS") {
                    text(0, ind[1], expression(Linear (OLS): y == a + bx + sigma*epsilon),
                         adj=0, cex = 1.2*txt.size)}
       if(lin.out@mmod == "lin" & lin.out@var.model=="POM") {
                    text(0, ind[1], expression(Linear (POM): y == a + bx + mu^theta*sigma*epsilon),
                         adj=0, cex = 1.2*txt.size)}
       if(lin.out@mmod == "quad" & lin.out@var.model=="OLS") {
		    text(0, ind[1], expression(Quadratic (OLS): y == a + bx + cx^2 + sigma*epsilon),
                         adj=0, cex=1.2*txt.size)}
       if(lin.out@mmod == "quad" & lin.out@var.model=="POM") {
		    text(0, ind[1], expression(Quadratic (POM): y == a + bx + cx^2 + mu^theta*sigma*epsilon),
                         adj=0, cex=1.2*txt.size)}
       text(0, ind[3], "Parameter", adj=0, cex=txt.size)
       text(0.2*xlim[2], ind[3], "Estimate",  adj=0, cex=txt.size)               
       text(0.4*xlim[2], ind[3], "Std. Error",  adj=0, cex=txt.size)
              
                if(lin.out@mmod == "lin"){
                  text(0, ind[4], expression(a (intercept)), adj=0, cex=txt.size)
                  text(0.2*xlim[2], ind[4], format(round(lin.out@coefficients[1],digits[1])), adj=0, cex=txt.size)
                  text(0.4*xlim[2], ind[4], format(round(lin.out@se.coefficients[1],digits[1])), adj=0, cex=txt.size)
                  
                  text(0, ind[5], expression(b (slope)), adj=0, cex=txt.size)
                  text(0.2*xlim[2], ind[5], format(round(lin.out@coefficients[2],digits[2])), adj=0, cex=txt.size)
                  text(0.4*xlim[2], ind[5], format(round(lin.out@se.coefficients[2],digits[2])), adj=0, cex=txt.size)
                  
                  text(0.2*xlim[2], ind[7], paste("sigma =", format(round(lin.out@sigma, digits[5]))), adj=0, cex=txt.size)
                  text(0, ind[7], paste("Rsq =", format(round(lin.out@Rsq, digits[5]))), adj=0, cex=txt.size)
                }
         
                if(lin.out@var.model=="POM"){
                  text(0.4*xlim[2], ind[7], paste("theta =", format(round(lin.out@theta, digits[4]))), adj=0, cex=txt.size)}

                text(0, ind[9], paste("MDC [m=", m, "] =", format(round(
                     mdc, digits[3]))), adj = 0, cex = txt.size)
                text(0.3*xlim[2], ind[9], paste("RDL [m=", m, "] =", format(round(
                     rdl, digits[3]))), adj = 0, cex = txt.size)
                text(0, ind[10], paste("LOQ [m=", m, ",CV=", loqcv, 
                     "] =", format(round(loq, digits[3]))), adj = 0, 
                     cex = txt.size)
              
         if(lof.print){
         lofl <- lin.out@lof.test
		  if(length(lofl@Fstat)>0){
		      text(0, ind[12], "Lack of Fit Test:", adj=0, cex=txt.size)
                      text(0, ind[13], paste("F =", format(round(lofl@Fstat, 3))), adj=0, cex=txt.size)
                      text(0.2*xlim[2], ind[13], paste("p-value =",  format(round(lofl@p.value, 3))),
                           adj=0, cex=txt.size)}}
                

                if(lin.out@mmod == "quad"){
                  text(0, ind[4], expression(a (intercept)), adj=0, cex=txt.size)
                  text(0.2*xlim[2], ind[4], format(round(lin.out@coefficients[1],digits[1])), adj=0, cex=txt.size)
                  text(0.4*xlim[2], ind[4], format(round(lin.out@se.coefficients[1],digits[1])), adj=0, cex=txt.size)
                  
                  text(0, ind[5], expression(b (lin)), adj=0, cex=txt.size)
                  text(0.2*xlim[2], ind[5], format(round(lin.out@coefficients[2],digits[2])), adj=0, cex=txt.size)
                  text(0.4*xlim[2], ind[5], format(round(lin.out@se.coefficients[2],digits[2])), adj=0, cex=txt.size)
                  
                  text(0, ind[6], expression(c (quad)), adj=0, cex=txt.size)
                  text(0.2*xlim[2], ind[6], format(round(lin.out@coefficients[3],digits[2])), adj=0, cex=txt.size)
                  text(0.4*xlim[2], ind[6], format(round(lin.out@se.coefficients[3],digits[2])), adj=0, cex=txt.size)
                  
                  text(0.2*xlim[2], ind[8], paste("sigma =", format(round(lin.out@sigma, digits[5]))), adj=0, cex=txt.size)
                  text(0, ind[8], paste("Rsq =", format(round(lin.out@Rsq, digits[5]))), adj=0, cex=txt.size) 
         
                  #if(lin.out@var.model=="POM"){
                    #text(0.4*xlim[2], ind[8], paste("theta =", format(round(lin.out@theta, digits[4]))), adj=0, cex=txt.size)}
                  
                  #text(0, ind[10], paste("MDC [m=", m, "] =", format(round(
                  #     mdc, digits[3]))), adj = 0, cex = txt.size)
                  #text(0.3*xlim[2], ind[10], paste("RDL [m=", m, "] =", format(round(
                  #     rdl, digits[3]))), adj = 0, cex = txt.size)
	 #text(0, ind[11], paste("LOQ [m=", m, ",CV=", loqcv, 
		#	"] =", format(round(loq, digits[3]))), adj = 0, 
        #                cex = txt.size)
         if(lof.print){
           lofl <- lin.out@lof.test
           if(length(lofl@Fstat)>0){
             #text(0, ind[13], "Lack of Fit Test:", adj=0, cex=txt.size)
             #text(0, ind[14], paste("F =", format(round(lofl@Fstat, 3))), adj=0, cex=txt.size)
             #text(0.2*xlim[2], ind[14], paste("p-value =",  format(round(lofl@p.value, 3))),
             #     adj=0, cex=txt.size)
                  }}
       }

         
              }
	invisible()
         }

diagplot <- function(pom.out, ptype = "all", ..., cex = 1, oma = c(1, 1, 3, 1))
{
###########################################################################
## diagnostic plot for the pom fit                                       ##
##                                                                       ##
## It seems like what is called the studentized residuals (which         ##
##  are also called jackknifed residuals) are                            ##
##  actually standard residuals:  standres = e_i/(s sqrt(1 - h_ii))      ##
##  while studres = (y_i - haty_(i))/sqrt(var(numerator)) =              ##
##  = (standres)*sqrt((n-p-standres^2)/(n-p-1))                          ##
##                                                                       ##
## comments added by Matthew Mitchell on 3/18/03 are indicated with ##   ##
## the original comments have ###                                        ##
###########################################################################


  
### if standard=T, plots standardized residuals versus predicted values for FPL or FPL-POM
### else plots regular residuals versus fitted (absolute residuals if absl=T)
	old.par <- par("mfrow", "cex", "oma")
	on.exit(par(old.par))
	if(ptype == "all")
		par(mfrow = c(2, 2))
	par(cex = cex)
	par(oma = oma)
	b <- pom.out@coefficients
	if(pom.out@vfe.method != "OLS")
		theta <- pom.out@theta
	else theta <- 0
	res <- pom.out@residuals
	if(pom.out@type == "fpl")
		fits <- fpl.model(pom.out@x, b[1], b[2], b[3], b[4], 
			parm = pom.out@parm)
	if(pom.out@type == "thpl")
		fits <- thpl.model(pom.out@x, b[1],b[2],b[3])
	if(pom.out@type == "lin")
		fits <- lin.model(pom.out@x, b, mmod = pom.out@mmod)
	xmat <- fits@gradient
	xprx.inv <- solve(t(xmat) %*% xmat)
	hats <- diag(xmat %*% xprx.inv %*% t(xmat))
	ress <- res/(((as.vector(fits)^(theta)) * pom.out@sigma) * sqrt(
		1 - hats))
	if(ptype == "all" | ptype == "raw") {
		plot(as.vector(fits), res, xlab = "Predicted Value", 
			ylab = "Residual", pch = 4, ...)
		abline(h = 0)
	}
	if(ptype == "all" | ptype == "stud") {
		plot(as.vector(fits), ress, xlab = "Predicted Value", 
			ylab = "Studentized Residual", pch = 4, ...)
		abline(h = 0)
	}
	if(ptype == "all" | ptype == "abs") {
		resa <- res/(((as.vector(fits)^(theta)) * pom.out@sigma))
		resa <- log(abs(resa))
		fitsa <- log(as.vector(fits))
		plot(fitsa, resa, ylab = "log (Absolute Residual)", 
			xlab = "log (Predicted Value)", pch = 4, ...)
## doesn't work in R, rreg is an SPlus command
## MM 4/1/03            abline(rreg(fitsa, resa), lty = 2)
	}
	if(ptype == "all" | ptype == "cube") {
		resc <- exp((1/3) * log(ress * ress))
		fitsc <- log(as.vector(fits))
		plot(fitsc, resc, ylab = 
			"Cube Root Squared Student Residual", xlab ## shortened title to fit MM 4/1/03
			 = "log (Predicted Value)", pch = 4, ...)
	}
	invisible()
}

precprof <- function(calib.fit.out,m=calib.fit.out@m,cv=calib.fit.out@cv,
            vlen=500,mit=1000,tolder=0.001,ylim,xlab="dose",ylab="CV",
            txt.size=1){#,dos=calib.fit.out@dos){
  if (calib.fit.out@type=="fpl"){
    precprof.fpl.pom(calib.fit.out, m = m, cv = cv, vlen = vlen, 
    mit = mit, toler = toler, ylim=ylim, xlab = xlab, ylab = ylab, 
    txt.size = txt.size)
  }
  else {if(calib.fit.out@type=="thpl"){
    precprof.thpl.pom(calib.fit.out, m = m, cv = cv, vlen = vlen, 
    mit = mit, toler = toler, ylim=ylim, xlab = xlab, ylab = ylab, 
    txt.size = txt.size)
  }
  else {if(calib.fit.out@type=="lin"){
    precprof.lin.pom(calib.fit.out, m = m, cv = cv, vlen = vlen, 
    ylim=ylim, xlab = xlab, ylab = ylab, txt.size = txt.size)
  }
      }
      }
}


precprof.fpl.pom <- function(fpl.out, m = fpl.out@m, cv = fpl.out@cv, vlen = 500, mit = 
	1000, toler = 0.001, ylim, xlab = "dose", ylab = "CV", txt.size
	 = 1,...)#, dos = fpl.out@dos, ...)
{

#########################################################
## This creates a plot of the dose vs. cv for output   ##
##  objects from fpl fit                               ##
## At the top of the graph the lowest cv acceptable,   ##
##  the working range of the assay, and the            ##
##  number of reps is printed                          ##
##                                                     ##
## not sure why loq from fpl is not used?              ##
##                                                     ##
## header added by Matthew Mitchell on 3/18/03         ##
#########################################################
  
	b <- fpl.out@coefficients
	cov.un <- fpl.out@cov.unscaled
	theta <- fpl.out@theta
	b50 <- ifelse(fpl.out@parm == 1, b[3], exp(b[3]))
	xpstart <- min(c(0.0005, min(fpl.out@x[fpl.out@x > 0])))
	xp <- c(seq(xpstart, b50, length = round(vlen/2, 0)), seq(b50, 
		max(fpl.out@x), length = round(vlen/2, 0)))
	yp <- as.vector(fpl.model(xp, b, parm = fpl.out@parm))
	dh.dy <- xp * ((b[2] - b[1])/(b[4] * (b[1] - yp) * (yp - b[2])))
	dh.db1 <- xp/(b[4] * (b[1] - yp))
	dh.db2 <- xp/(b[4] * (yp - b[2]))
	if(fpl.out@parm == 1)
		dh.db3 <- xp/b[3]
	else dh.db3 <- xp
	dh.db4 <- ( - xp/(b[4] * b[4])) * log((b[1] - yp)/(yp - b[2]))
	sigma2 <- fpl.out@sigma * fpl.out@sigma
	var.xnot.hat <- (((dh.dy * dh.dy) * sigma2 * (yp^(2 * theta)))/
		m + sigma2 * (dh.db1 * (dh.db1 * cov.un[1, 1] + dh.db2 * 
		cov.un[2, 1] + dh.db3 * cov.un[3, 1] + dh.db4 * cov.un[
		4, 1]) + dh.db2 * (dh.db1 * cov.un[1, 2] + dh.db2 * 
		cov.un[2, 2] + dh.db3 * cov.un[3, 2] + dh.db4 * cov.un[
		4, 2]) + dh.db3 * (dh.db1 * cov.un[1, 3] + dh.db2 * 
		cov.un[2, 3] + dh.db3 * cov.un[3, 3] + dh.db4 * cov.un[
		4, 3]) + dh.db4 * (dh.db1 * cov.un[1, 4] + dh.db2 * 
		cov.un[2, 4] + dh.db3 * cov.un[3, 4] + dh.db4 * cov.un[
		4, 4])))
	sd <- sqrt(var.xnot.hat)
	plot.y <- sd/xp
	if(missing(ylim)) {
		ylim <- c(0, 2 * cv)
		if(ylim[2] < min(plot.y))
			ylim <- c(0, 1.5 * (min(plot.y)))
	}
	plot.x <- xp[plot.y < max(ylim)]
	plot.y <- plot.y[plot.y < max(ylim)]
	
	if((length(plot.x)==0)|length(plot.y)==0)
		stop("Error was produced in precprof estimated CV larger CV cutoff")
	
	plot(plot.x, plot.y, type = "l", xlab = xlab, ylab = ylab, ylim
		 = ylim, ...)
	abline(h = cv, lty = 3)
	xmin <- xp[sd/xp == min(sd/xp)]
	if(cv > min(sd/xp[xp != 0])) {
		xl <- xp[xp <= xmin]
		xu <- xp[xp >= xmin]

			sdl <- sd[xp <= xmin]
			sdu <- sd[xp >= xmin]
			loqsl <- approx(cv - (sdl/xl), xl, 0,ties="ordered")$y
			loqsu <- approx(cv - (sdu/xu), xu, 0,ties="ordered")$y
			lims <- c(loqsl, ifelse(is.na(loqsu), max(xp), 
				loqsu))
       
		text(min(plot.x) + 0.1 * (max(plot.x) - min(plot.x)), 
			0.9 * max(ylim), paste("CV acceptable =", cv, 
			"\n\n\nWorking range of assay: ", format(round(
			min(lims), 1)), " to ", format(round(max(lims), 
			1)), "\n\n\n", 
			"Number of replicates based on: ", format(round(
			m, 0))), adj = 0, cex = txt.size)
		lines(rep(lims[1], 2), c(0, cv), lty = 7)
		lines(rep(lims[2], 2), c(0, cv), lty = 7)
              }
	else {
		text(min(plot.x) + 0.1 * (max(plot.x) - min(plot.x)), 
			0.9 * max(ylim), paste(
			"Acceptable precision not acheived in assay range for m=",
			m, ",  CV=", cv, sep = ""), adj = 0, cex = 
			txt.size)
	}
	invisible()
}

precprof.thpl.pom <- function(thpl.out, m = thpl.out@m, cv = thpl.out@cv, vlen = 500, mit = 
	1000, toler = 0.001, ylim, xlab = "dose", ylab = "CV", txt.size
	 = 1,...)#, dos = fpl.out@dos, ...)
{

#########################################################
## This creates a plot of the dose vs. cv for output   ##
##  objects from fpl fit                               ##
## At the top of the graph the lowest cv acceptable,   ##
##  the working range of the assay, and the            ##
##  number of reps is printed                          ##
##                                                     ##
## not sure why loq from fpl is not used?              ##
##                                                     ##
## header added by Matthew Mitchell on 3/18/03         ##
#########################################################
  
	b <- thpl.out@coefficients
	cov.un <- thpl.out@cov.unscaled
	theta <- thpl.out@theta
	b50 <- b[3]
	xpstart <- min(c(0.0005, min(thpl.out@x[thpl.out@x > 0])))
	xp <- c(seq(xpstart, b50, length = round(vlen/2, 0)), seq(b50, 
		max(thpl.out@x), length = round(vlen/2, 0)))
	yp <- as.vector(thpl.model(xp, b))
	dh.dy <- xp * ((b[2] - b[1])/(b[1] - yp) * (yp - b[2]))
	dh.db1 <- xp/(b[1] - yp)
	dh.db2 <- xp/(yp - b[2])
	dh.db3 <- xp/b[3]
	sigma2 <- thpl.out@sigma^2
	var.xnot.hat <- (((dh.dy * dh.dy) * sigma2 * (yp^(2 * theta)))/
		m + sigma2 * (dh.db1 * (dh.db1 * cov.un[1, 1] + dh.db2 * 
		cov.un[2, 1] + dh.db3 * cov.un[3, 1]) + dh.db2 * (dh.db1 * cov.un[1, 2] + dh.db2 * 
		cov.un[2, 2] + dh.db3 * cov.un[3, 2]) + dh.db3 * (dh.db1 * cov.un[1, 3] + dh.db2 * 
		cov.un[2, 3] + dh.db3 * cov.un[3, 3])))
	sd <- sqrt(var.xnot.hat)
	plot.y <- sd/xp
	if(missing(ylim)) {
		ylim <- c(0, 2 * cv)
		if(ylim[2] < min(plot.y))
			ylim <- c(0, 1.5 * (min(plot.y)))
	}
	plot.x <- xp[plot.y < max(ylim)]
	plot.y <- plot.y[plot.y < max(ylim)]
	
	if((length(plot.x)==0)|length(plot.y)==0)
		stop("Error was produced in precprof estimated CV larger CV cutoff")
		
	plot(plot.x, plot.y, type = "l", xlab = xlab, ylab = ylab, ylim
		 = ylim, ...)
	abline(h = cv, lty = 3)
	xmin <- xp[sd/xp == min(sd/xp)]
	if(cv > min(sd/xp[xp != 0])) {
		xl <- xp[xp <= xmin]
		xu <- xp[xp >= xmin]
		sdl <- sd[xp <= xmin]
		sdu <- sd[xp >= xmin]
		loqsl <- approx(cv - (sdl/xl), xl, 0,ties="ordered")$y
		loqsu <- approx(cv - (sdu/xu), xu, 0,ties="ordered")$y
		lims <- c(loqsl, ifelse(is.na(loqsu), max(xp), 
			loqsu))
		text(min(plot.x) + 0.1 * (max(plot.x) - min(plot.x)), 
			0.9 * max(ylim), paste("CV acceptable =", cv, 
			"\n\n\nWorking range of assay: ", format(round(
			min(lims), 1)), " to ", format(round(max(lims), 
			1)), "\n\n\n", 
			"Number of replicates based on: ", format(round(
			m, 0))), adj = 0, cex = txt.size)
		lines(rep(lims[1], 2), c(0, cv), lty = 7)
		lines(rep(lims[2], 2), c(0, cv), lty = 7)
              }
	else {
		text(min(plot.x) + 0.1 * (max(plot.x) - min(plot.x)), 
			0.9 * max(ylim), paste(
			"Acceptable precision not acheived in assay range for m=",
			m, ",  CV=", cv, sep = ""), adj = 0, cex = 
			txt.size)
	}
	invisible()
}

precprof.lin.pom <- function(lin.out, m = lin.out@m, cv = lin.out@cv, vlen = 500, ylim, 
	xlab = "dose", ylab = "CV", txt.size = 1, ...)
{
#########################################################
## This creates a plot of the dose vs. cv for output   ##
##  objects from lin fit                               ##
## At the top of the graph the lowest cv acceptable,   ##
##  the working range of the assay, and the            ##
##  number of reps is printed                          ##
##                                                     ##
## not sure why loq from fpl is not used?              ##
##                                                     ##
## header added by Matthew Mitchell on 3/18/03         ##
## original comments hae # or  ###                     ##
#########################################################
  
	beta <- lin.out@coefficients
	cov.un <- lin.out@cov.unscaled
	sigma <- lin.out@sigma
	theta <- lin.out@theta
	x <- lin.out@x
	x.range <- c(ifelse(min(lin.out@x) == 0, min(lin.out@x) + 0.005,
		min(lin.out@x)), max(lin.out@x))
	xp <- exp(seq(log(min(min(lin.out@x[lin.out@x > 0]), 1e-010)), 
		log(max(lin.out@x)), length = vlen))
	sigma2 <- sigma * sigma	
	### calculate the standard deviation for xps
	if(lin.out@mmod == "lin") {
		yp <- beta[1] + beta[2] * xp
		var.xnot.hat <- (sigma2) * ((yp^(2 * theta))/(beta[2] * 
			beta[2] * m) + (cov.un[1, 1] + 2 * cov.un[1, 2] *
			xp + cov.un[2, 2] * xp * xp)/(beta[2]^2))
		lin.ind <- 1
	}
	if(lin.out@mmod == "quad") {
		yp <- beta[1] + beta[2] * xp + beta[3] * xp * xp
		dh.dyp <- 1/(xp + (beta[2]/(2 * beta[1])))
		dh.da <- beta[2]/(2 * beta[1] * beta[1]) + ((yp - beta[
			3])/beta[1]) * (1/(xp + (beta[2]/(2 * beta[1]))
			)) - (xp + (beta[2]/(2 * beta[1])))/(2 * beta[1
			] * beta[1])
		dh.db <- -1/(2 * beta[1]) + (beta[2]/(2 * beta[1])) * (
			1/(xp + (beta[2]/(2 * beta[1]))))
		dh.dc <- -1/(xp + (beta[2]/(2 * beta[1])))
		var.xnot.hat <- ((dh.dyp * dh.dyp) * sigma2 * (yp^(2 * 
			lin.out@theta)))/m + sigma2 * (dh.da * (dh.da * 
			cov.un[1, 1] + dh.db * cov.un[1, 2] + dh.dc * 
			cov.un[1, 3]) + dh.db * (dh.da * cov.un[1, 2] + 
			dh.db * cov.un[2, 2] + dh.dc * cov.un[2, 3]) + 
			dh.dc * (dh.da * cov.un[1, 3] + dh.db * cov.un[
			2, 3] + dh.dc * cov.un[3, 3]))
		lin.ind <- 0
	}
	sd <- sqrt(var.xnot.hat)
	if(missing(ylim))
		ylim <- c(0, 2 * cv)
	plot.y <- sd/xp
	plot.x <- xp[plot.y < max(ylim)]
	plot.y <- plot.y[plot.y < max(ylim)]
	
	if((length(plot.x)==0)|length(plot.y)==0)
		stop("Error was produced in precprof estimated CV larger CV cutoff")
		
	plot(plot.x, plot.y, type = "l", xlab = xlab, ylab = ylab, ylim
		 = ylim, ...)	#  plot(xp, sd/xp, 
#       type = "l", 
#       main= main, 
#       ylim=ylim, 
#       xlab = xlab, 
#       ylab = ylab, 
#       cex=1)	
	abline(h = cv, lty = 3)	
	### mtext("    Acceptable\n   precision", side = 4, outer = T,	  at = cv, cex=cex)
### xmin is the inflection point (if one exists), else it is the maximum of x
	est.cv <- sd/xp
	xmin <- xp[sd/xp == min(est.cv)]	
	### xl and xu are the starting endpoint values for the search algorithm
	xl <- c(min(xp), xmin)
	xu <- c(xmin, max(x.range))
	if(xu[1] < xl[1]) {
		tempx <- xl[1]
		xl[1] <- xu[1]
		xu[1] <- tempx
	}
	svec <- c(cov.un[, 1], cov.un[, 2])	
	### if reasonable find the endpoints for the x values which have CV lower
###   than acceptable level (cv)
	if(cv > min(est.cv[xp != 0])) {
		xl <- xp[xp <= xmin]
		xu <- xp[xp >= xmin]
		sdl <- sd[xp <= xmin]
		sdu <- sd[xp >= xmin]
		loqsl <- approx(cv - (sdl/xl), xl, 0)$y
		if(est.cv[length(est.cv)] < cv)
			loqsu <- xu[length(xu)]
		else loqsu <- approx(cv - (sdu/xu), xu, 0)$y
		lims <- c(loqsl, ifelse(is.na(loqsu), max(xp), loqsu))	
	#    }
#      else
#	{
#	  if (!loaded) dyn.load("/local/lib/calib/sun4/src/precprof.pom.o")
#	  snum <-2
#	  uprange <- NULL
#	  if (est.cv[length(est.cv)]<cv)
#	    {
#	      xl <- xl[1]
#	      xu <- xu[1]
#	      snum <-1
#	      uprange <- max(xp)
#	    }
#	  
#	  bp.out <-.C("pomprec", 
#		      xl=as.double(xl),xu=as.double(xu), 
#		      as.double(cv),as.double(beta), 
#		      as.double(svec),as.integer(maxit), 
#		      as.double(toler),as.double(theta),as.double(sigma), 
#		      as.integer(snum),as.integer(m),as.integer(lin.ind))
#	  lims <- c(bp.out$xl,uprange)
#	}
### print out the results on the graph
		text(min(plot.x) + 0.1 * (max(plot.x) - min(plot.x)), 
			0.9 * max(ylim), paste(
			"Acceptable precision (CV) =", cv, 
			"\nWorking range of assay:  \n", format(round(
			min(lims), 1)), " to ", format(round(max(lims), 
			1)), "\n", "Based on", format(round(m, 0)), 
			"replicate(s)"), adj = 0, cex = txt.size)
		lines(rep(lims[1], 2), c(0, cv), lty = 7)
		lines(rep(lims[2], 2), c(0, cv), lty = 7)
	}
	else text(min(plot.x) + 0.1 * (max(plot.x) - min(plot.x)), 0.9 * 
			max(ylim), paste("Acceptable precision (CV) = ",
			cv, "\n", 
			"Acceptable precision not achieved in assay range",
			sep = ""), adj = 0, cex = txt.size)
	invisible()
}


