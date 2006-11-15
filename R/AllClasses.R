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

setClass("lof.test",
         representation(Fstat="numeric",
                        p.value="numeric",
                        lofss="numeric",
                        df.lof="numeric",
                        pure.error="numeric",
                        df.pure.error="numeric",
                        sse="numeric",
                        df.sse="intOrNum"))

setClass("calib.fit",
         representation(coefficients="numeric",
                        se.coefficients="numeric",
                        sigma="numeric",
                        cov.unscaled="matrix",
                        theta="numeric",
                        df.residual="intOrNum",
                        fitted.values="numeric",
                        residuals="numeric",
                        vfe.method="character",
                        kused="numeric",
                        status="character",
                        x="numOrMat",
                        y="numOrMat",
                        parm="numeric",
                        m="numeric",
                        cv="numeric",
                        mdc="numOrNULL",
                        rdl="numOrNULL",
                        loq="numOrNULL",
                        cf="numeric",
                        #dos="logical", 04-17-07 Removing from calib
                        ## does not seem to be doing anything of use
                        gradient="matrix",
                        lof.test="lof.test",
                        var.model="character",
                        conf.level="numeric",
                        mmod="character",
                        SST="numeric",
                        SSE="numeric",
                        SSR="numeric",
                        Rsq="numeric",
                        type="character",
						rdlwarn="character"
                        ))

setClass("Calib",
         representation(Estimated.x="numeric",
                        PredStdErr="numeric",
                        inver.low="numOrNULL",
                        inver.up="numOrNULL",
                        wald.low="numOrNULL",
                        wald.up="numOrNULL",
                        waldl.low="numOrNULL",
                        waldl.up="numOrNULL",
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
                        samp.names="charOrNULL",
                        lof.test="lof.test",
                        calib.fit="calib.fit"
                        )
         )
