plotPost = function( paramSampleVec , credMass=0.95 , compVal=NULL ,
           HDItextPlace=0.7 , ROPE=NULL , yaxt=NULL , ylab=NULL ,
           xlab=NULL , cex.lab=NULL , cex=NULL , xlim=NULL , main=NULL ,
           col=NULL , border=NULL , showMode=F , showCurve=F , ... ) {
    # Override defaults of hist function, if not specified by user:
    # (additional arguments "..." are passed to the hist function)
    if ( is.null(xlab) ) xlab="Parameter"
    if ( is.null(cex.lab) ) cex.lab=1.5
    if ( is.null(cex) ) cex=1.4
    if ( is.null(xlim) ) xlim=range( c( compVal , paramSampleVec ) )
    if ( is.null(main) ) main=""
    if ( is.null(yaxt) ) yaxt="n"
    if ( is.null(ylab) ) ylab=""
    if ( is.null(col) ) col="grey"
    if ( is.null(border) ) border="white"
    # Plot histogram.
    if ( !showCurve ) {
      par(xpd=NA)
      histinfo = hist( paramSampleVec , xlab=xlab , yaxt=yaxt , ylab=ylab ,
                       freq=F , border=border , col=col ,
                       xlim=xlim , main=main , cex=cex , cex.lab=cex.lab ,
                       ... )
    }
    if ( showCurve ) {
      histinfo = hist( paramSampleVec , plot=F )
      densCurve = density( paramSampleVec , adjust=2 )
      plot( densCurve , type="l" , lwd=5 , col=col , 
            xlim=xlim , xlab=xlab , yaxt=yaxt , ylab=ylab ,
            main=main , cex=cex , cex.lab=cex.lab , ... )
    }
    # Display mean or mode:
    if ( showMode==F ) {
        meanParam = mean( paramSampleVec )
        text( meanParam , .9*max(histinfo$density) ,
              bquote(mean==.(signif(meanParam,3))) , adj=c(.5,0) , cex=cex )
    } else {
        dres = density( paramSampleVec )
        modeParam = dres$x[which.max(dres$y)]
        text( modeParam , .9*max(histinfo$density) ,
              bquote(mode==.(signif(modeParam,3))) , adj=c(.5,0) , cex=cex )
    }
    # Display the comparison value.
    if ( !is.null( compVal ) ) {
      cvCol = "darkgreen"
      pcgtCompVal = round( 100 * sum( paramSampleVec > compVal )
                            / length( paramSampleVec )  , 1 )
       pcltCompVal = 100 - pcgtCompVal
       lines( c(compVal,compVal) , c(.5*max(histinfo$density),0) ,
              lty="dashed" , lwd=2 , col=cvCol )
       text( compVal , .5*max(histinfo$density) ,
             bquote( .(pcltCompVal)*"% <= " *
                     .(signif(compVal,3)) * " < "*.(pcgtCompVal)*"%" ) ,
             adj=c(pcltCompVal/100,-0.2) , cex=cex , col=cvCol )
    }
    # Display the ROPE.
    if ( !is.null( ROPE ) ) {
      ropeCol = "darkred"
       pcInROPE = ( sum( paramSampleVec > ROPE[1] & paramSampleVec < ROPE[2] )
                            / length( paramSampleVec ) )
       ROPEtextHt = .35*max(histinfo$density)
       lines( c(ROPE[1],ROPE[1]) , c(ROPEtextHt,0) , lty="dotted" , lwd=2 , col=ropeCol )
       lines( c(ROPE[2],ROPE[2]) , c(ROPEtextHt,0) , lty="dotted" , lwd=2 , col=ropeCol)
       text( mean(ROPE) , ROPEtextHt ,
             bquote( .(round(100*pcInROPE))*"% in ROPE" ) ,
             adj=c(.5,-0.2) , cex=1 , col=ropeCol )
    }
    # Display the HDI.
    source("HDIofMCMC.R")
    HDI = HDIofMCMC( paramSampleVec , credMass )
    lines( HDI , c(0,0) , lwd=4 )
    text( mean(HDI) , 0 , bquote(.(100*credMass) * "% HDI" ) ,
          adj=c(.5,-1.9) , cex=cex )
    text( HDI[1] , 0 , bquote(.(signif(HDI[1],3))) ,
          adj=c(HDItextPlace,-0.5) , cex=cex )
    text( HDI[2] , 0 , bquote(.(signif(HDI[2],3))) ,
          adj=c(1.0-HDItextPlace,-0.5) , cex=cex )
    par(xpd=F)
    return( histinfo )
}
