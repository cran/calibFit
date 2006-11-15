## calibplotMethods.R
## 
## Plot methods for calib.fit and calib objects
##
## Author: samarov
###############################################################################

## ============================================================================
## Plot methods for calib objects
## =======================================================================
setMethod("plot",signature(x="Calib",y="missing"),
          function(x,...) plot.calib(x,...))
## ============================================================================
## Plot mehods for calib.fit objects
## =======================================================================
setMethod("plot",signature(x="calib.fit",y="missing"),
	function(x,type,...){
    	if (missing(type)){
        	if (x@type=="fpl"){
            	plot.fpl.pom(x,...)
              }
              else {if (x@type=="thpl"){
                plot.thpl.pom(x,...)
              }
              else {if (x@type=="lin"){
                plot.lin.pom(x,...)
              }
              }
              }
		}
        else {if (type=="precprof"){
        	if (x@type=="fpl")
		    	precprof.fpl.pom(x,...)
		        else {if (x@type=="lin")
                	precprof.lin.pom(x,...)
				else {if (x@type=="thpl")
                	precprof.thpl.pom(x,...)
                }
                }
		}
        else {if (type=="diagplot"){
        	diagplot(x,...)
            }
		}
		}
	})
