# For the graphical displays at the end of the program to run, you must have the
# auxiliary programs plotPost.R and HDIofMCMC.R in the same folder as this
# program, and set R's working directory to this folder. These programs are
# available from http://www.indiana.edu/~kruschke/DoingBayesianDataAnalysis/

graphics.off()
rm(list=ls(all=TRUE))
fileNameRoot="SplitPlotJags" # for constructing output filenames
if ( .Platform$OS.type != "windows" ) { 
  windows <- function( ... ) X11( ... ) 
}
require(rjags)
# Programmed by John Kruschke, Dec 2011 - Jan 2012, in the style of
# Kruschke, J. K. (2011). Doing Bayesian data analysis:
# A Tutorial with R and BUGS. Academic Press / Elsevier.

#------------------------------------------------------------------------------
# THE MODEL.

# Split-plot design, with single between-subject factor and 
# single within-subject factor.
# Factor s is nested within factor a.
# Factor b has s repeated across its levels.

modelstring = "
model {
  for ( i in 1:Ntotal ) {
    y[i] ~ dnorm( mu[i] , tau )
    mu[i] <- base + a[aLvl[i]] + s[sLvl[i]] + b[bLvl[i]] + axb[aLvl[i],bLvl[i]]
  }
  # Prior:
  tau <- pow( sigma , -2 )
  sigma ~ dunif(0,1000) # change to be appropriate for scale of data!
  #
  base ~ dnorm(500,1.0E-6) # change to be appropriate for scale of data!
  #
  for ( j in 1:NaLvl ) { a[j] ~ dnorm( 0.0 , aTau ) }
  aTau <- 1 / pow( aSD , 2 )
  aSD <- abs( aSDunabs ) + .1
  aSDunabs ~ dt( 0 , 1.0E-7 , 2 )
  #
  for ( j in 1:NbLvl ) { b[j] ~ dnorm( 0.0 , bTau ) }
  bTau <- 1 / pow( bSD , 2 )
  bSD <- abs( bSDunabs ) + .1
  bSDunabs ~ dt( 0 , 1.0E-7 , 2 )
  #
  for ( j in 1:NsLvl ) { s[j] ~ dnorm( 0.0 , sTau ) }
  sTau <- 1 / pow( sSD , 2 )
  sSD <- abs( sSDunabs ) + .1
  sSDunabs ~ dt( 0 , 1.0E-7 , 2 )
  #
  for ( j in 1:NaLvl ) { for ( k in 1:NbLvl ) {
    axb[j,k] ~ dnorm( 0.0 , axbTau )
  } }
  axbTau <- 1 / pow( axbSD , 2 )
  axbSD <- abs( axbSDunabs ) + .1
  axbSDunabs ~ dt( 0 , 1.0E-7 , 2 )

  # Conversion to sum-to-zero effects happens outside JAGS in R
}
" # close quote for modelstring.
writeLines(modelstring,con="model.txt") # Write model to a file.

#------------------------------------------------------------------------------
# THE DATA.

# Specify data source:
dataSource = c( "MaxwellDelaneyCh12data" )[1]

# Load the data:
if ( dataSource == "MaxwellDelaneyCh12data" ) {
  fileNameRoot = paste( fileNameRoot , dataSource , sep="" )
  datarecord = read.csv( "MaxwellDelaneyCh12data.csv" )
  # Factor s is nested within factor a.
  # Factor b has s repeated across its levels.
  # Model assumes that the level codes are consecutive integers starting at 1.
  # Rename variables to generic factor names:  
  y =  as.numeric(datarecord$Y)
  sName = "Subj" # must match column name in datarecord
  sLvl = as.numeric(datarecord[,sName])
  aName = "Age" # must match column name in datarecord
  aLvl = as.numeric(datarecord[,aName])
  aLvlNames = levels(datarecord[,aName])
  bName = "Angle" # must match column name in datarecord
  bLvl = as.numeric(datarecord[,bName])
  bLvlNames = levels(datarecord[,bName])
  # The remainder is standard for any data:
  Ntotal = length(y)
  NsLvl = length(unique(sLvl))
  NaLvl = length(unique(aLvl))
  NbLvl = length(unique(bLvl))
  # Specify desired contrasts:
  aContrastList = list( A2vA1=c(-1,1) )
  bContrastList = list( B2vB1=c(-1,1,0) ,  B3vB2=c(0,-1,1) )
  axbContrastList = list( A1B1vA2B1=matrix(c(1,0,0,-1,0,0), 
                                           byrow=T,nrow=NaLvl,ncol=NbLvl),
                          DiffQuadTrend=matrix(c(-1,2,-1,1,-2,1),
                                               byrow=T,nrow=NaLvl,ncol=NbLvl))
  # Create aLvlOFsLvl , the 'a' level of a given 's' level. 
  # It is used for conversion to sum-to-zero effects.
  aLvlOFsLvl = rep( 0 , NsLvl )
  for ( sIdx in 1:NsLvl ) {
    aLvlOFsLvl[sIdx] = unique( aLvl[ sLvl==sIdx ] ) # assumes there is unique aLVL
  }
}

# Specify the data in a form that is compatible with model, as a list:
#ySDorig = sd(y)
#yMorig = mean(y)
#zy = ( y - yMorig ) / ySDorig
dataList = list(
  y = y , # zy ,
  aLvl = aLvl ,
  bLvl = bLvl ,
  sLvl = sLvl ,
  Ntotal = Ntotal ,
  NsLvl = NsLvl ,
  NaLvl = NaLvl ,
  NbLvl = NbLvl
)

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

# Use data to determine reasonable initial values.
# Compute ab cell means, marginal, and grand mean:
abMean = matrix( aggregate( dataList$y , list( aLvl , bLvl ) , mean )[,"x"] , 
                 nrow = NaLvl , ncol=NbLvl , byrow=F )
aMean = rowMeans(abMean)
bMean = colMeans(abMean)
gMean = mean(abMean)
# Convert into a, b, axb effects:
aEff = aMean - gMean
bEff = bMean - gMean
axbEff = abMean - ( gMean + outer( aEff , bEff , "+" ) )
# Compute subject effect 
sMean = aggregate( dataList$y , list( sLvl ) , mean )[,"x"]
sEff = 0*sMean
for ( sIdx in 1:NsLvl ) {
  sEff[sIdx] = sMean[sIdx] - aMean[aLvlOFsLvl[sIdx]]
}

initsList = list(
            base = gMean ,
            a = aEff ,
            b = bEff ,
            s = sEff ,
            axb = axbEff ,
            sigma = sd(dataList$y)/2 , # lazy
            aSDunabs = sd(aEff) ,
            bSDunabs = sd(bEff) ,
            axbSDunabs = sd(as.vector(axbEff))
        )

#------------------------------------------------------------------------------
# RUN THE CHAINS

parameters = c( "base" ,  "a" ,  "b" ,  "s" , "axb" ,
                "sigma" , "aSD" , "bSD" , "sSD" , "axbSD" )
adaptSteps = 1000 
burnInSteps = 2000 
nChains = 10 
numSavedSteps=100000 # Requires a long chain for stable results. More is better! 
thinSteps=1
nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )
# Create, initialize, and adapt the model:
jagsModel = jags.model( "model.txt" , data=dataList , inits=initsList , 
                        n.chains=nChains , n.adapt=adaptSteps )
# Burn-in:
cat( "Burning in the MCMC chain...\n" )
update( jagsModel , n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
codaSamples = coda.samples( jagsModel , variable.names=parameters , 
                            n.iter=nIter , thin=thinSteps )
# resulting codaSamples object has these indices: 
#   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS
checkConvergence = F  #  sum-to-zero is not imposed yet.
if ( checkConvergence ) {
  show( summary( codaSamples ) )
  plot( codaSamples , ask=T )
  autocorr.plot( codaSamples , ask=T )
}

source("plotPost.R")

# Convert coda-object codaSamples to matrix object for easier handling.
# But note that this concatenates the different chains into one long chain.
# Result is mcmcChain[ stepIdx , paramIdx ]
mcmcChain = as.matrix( codaSamples )
chainLength = NROW(mcmcChain)

# Extract parameter values:
# "base" ,  "a" ,  "b" ,  "s" , "axb" , "sigma" , "aSD" , "bSD" , "axbSD"
base = mcmcChain[,"base"]
sigma = mcmcChain[,"sigma"]
aSD = mcmcChain[,"aSD"]
bSD = mcmcChain[,"bSD"]
sSD = mcmcChain[,"sSD"]
axbSD = mcmcChain[,"axbSD"]
a = array( 0 , dim=c( NaLvl , chainLength ) )
for ( j in 1:NaLvl ) { a[j,] = mcmcChain[, paste("a[",j,"]",sep="") ] }
b = array( 0 , dim=c( NbLvl , chainLength ) )
for ( j in 1:NbLvl ) { b[j,] = mcmcChain[, paste("b[",j,"]",sep="") ] }
s = array( 0 , dim=c( NsLvl , chainLength ) )
for ( j in 1:NsLvl ) { s[j,] = mcmcChain[, paste("s[",j,"]",sep="") ] }
axb = array( 0 , dim=c( NaLvl , NbLvl , chainLength ) )
for ( j in 1:NaLvl ) { 
  for ( k in 1:NbLvl ) {
    axb[j,k,] = mcmcChain[, paste("axb[",j,",",k,"]",sep="") ] 
  }
}

# Convert parameters to xum-to-zero values:
# Define placeholders for STZ values:
baseSTZ = base ; aSTZ=a ; bSTZ=b ; sSTZ=s ; axbSTZ=axb
# Define placeholder for predicted Y at each step of chain:
predY = rep(0,Ntotal)
# Now go through each step of chain and compute STZ values:
cat("\nConverting MCMC to sum-to-zero (STZ) values.\n")
cat("(If you find a way to speed this up, let me know!)\n")
cat("Percent done:     ")
displayIdx = round(seq(chainLength/50,chainLength,by=chainLength/50))
for ( mcIdx in 1:chainLength ) {
  if ( any(mcIdx==displayIdx) ) { 
    cat("\b\b\b\b. ",round(100*mcIdx/chainLength)) 
  }
  # cConstruct predicted values:
  for ( drIdx in 1:length(predY) ) {
    predY[drIdx] = ( base[mcIdx] 
                    + a[aLvl[drIdx],mcIdx]
                    + b[bLvl[drIdx],mcIdx]
                    + s[sLvl[drIdx],mcIdx]
                    + axb[aLvl[drIdx],bLvl[drIdx],mcIdx] )
  }
  # Compute ab cell means, marginal, and grand mean:
  abMean = matrix( aggregate( predY , list( aLvl , bLvl ) , mean )[,"x"] , 
                   nrow = NaLvl , ncol=NbLvl , byrow=F )
  aMean = rowMeans(abMean)
  bMean = colMeans(abMean)
  gMean = mean(abMean)
  # Convert into a, b, axb effects:
  aEff = aMean - gMean
  bEff = bMean - gMean
  axbEff = abMean - ( gMean + outer( aEff , bEff , "+" ) )
  # Compute subject effect 
  sMean = aggregate( dataList$y , list( sLvl ) , mean )[,"x"]
  sEff = 0*sMean
  for ( sIdx in 1:NsLvl ) {
    sEff[sIdx] = sMean[sIdx] - aMean[aLvlOFsLvl[sIdx]]
  }
  # Now store results:
  baseSTZ[mcIdx] = gMean
  aSTZ[,mcIdx] = aEff
  bSTZ[,mcIdx] = bEff
  sSTZ[,mcIdx] = sEff
  axbSTZ[,,mcIdx] = axbEff
}

# Now make plots of results
aColor = "pink"
bColor = "lawngreen"
axbColor = "lavender"

# Display the SD's
windows()
layout( matrix(c(1:6),nrow=2,byrow=T) )
par( mar=c(3,1,2.5,0) , mgp=c(2,0.7,0) )
histInfo = plotPost( aSD , xlab="aSD" , main=paste(aName,"SD") , 
                     breaks=30 , col=aColor , showMode=T )
histInfo = plotPost( bSD , xlab="bSD" , main=paste(bName,"SD") , 
                     breaks=30 , col=bColor , showMode=T )
histInfo = plotPost( axbSD , xlab="axbSD" , main=paste(aName,"x",bName,"SD") , 
                     breaks=30 , col=axbColor , showMode=T )
histInfo = plotPost( sigma , xlab="sigma" , main="Within Cell SD" , 
                     breaks=30 , col="skyblue" , showMode=T )
histInfo = plotPost( sSD , xlab="sSD" , main=paste(sName,"[within",aName,"] SD") , 
                     breaks=30 , col="skyblue" , showMode=T )
plot(-1,-1,xlim=c(0,2),ylim=c(0,2),type="n",xaxt="n",xlab="",yaxt="n",ylab="")
text(1,1,"Caution: Even long MCMC chains \nproduce unstable estimates of SD")
savePlot( file=paste(fileNameRoot,"SDmode.eps",sep="") , type="eps" )
savePlot( file=paste(fileNameRoot,"SDmode.jpg",sep="") , type="jpg" )
windows()
layout( matrix(c(1:6),nrow=2,byrow=T) )
par( mar=c(3,1,2.5,0) , mgp=c(2,0.7,0) )
histInfo = plotPost( aSD , xlab="aSD" , main=paste(aName,"SD") , 
                     breaks=30 , col=aColor , showMode=F )
histInfo = plotPost( bSD , xlab="bSD" , main=paste(bName,"SD") , 
                     breaks=30 , col=bColor , showMode=F )
histInfo = plotPost( axbSD , xlab="axbSD" , main=paste(aName,"x",bName,"SD") , 
                     breaks=30 , col=axbColor , showMode=F )
histInfo = plotPost( sigma , xlab="sigma" , main="Within Cell SD" , 
                     breaks=30 , col="skyblue" , showMode=F )
histInfo = plotPost( sSD , xlab="sSD" , main=paste(sName,"[within",aName,"] SD") , 
                     breaks=30 , col="skyblue" , showMode=F )
plot(-1,-1,xlim=c(0,2),ylim=c(0,2),type="n",xaxt="n",xlab="",yaxt="n",ylab="")
text(1,1,"Caution: Even long MCMC chains \nproduce unstable estimates of SD")
savePlot( file=paste(fileNameRoot,"SDmean.eps",sep="") , type="eps" )
savePlot( file=paste(fileNameRoot,"SDmean.jpg",sep="") , type="jpg" )


# Display a, b, axb values:
windows((NaLvl+1)*2.75,(NbLvl+1)*2.0)
layoutMat = matrix( 0 , nrow=(NaLvl+1) , ncol=(NbLvl+1) )
layoutMat[1,1] = 1 # base
lastCellIdx = layoutMat[1,1]
layoutMat[1+1:NaLvl,1] = 1:NaLvl + lastCellIdx # a
lastCellIdx = layoutMat[NaLvl+1,1]
layoutMat[1,1+1:NbLvl] = 1:NbLvl + lastCellIdx # b
lastCellIdx = layoutMat[1,NbLvl+1]
layoutMat[2:(NaLvl+1),2:(NbLvl+1)] = matrix( lastCellIdx + 1:(NaLvl*NbLvl) ,
    nrow=NaLvl , ncol=NbLvl , byrow=F )
layout( layoutMat )
par( mar=c(4,0.5,2.5,0.5) , mgp=c(2,0.7,0) )
histinfo = plotPost( baseSTZ , xlab=expression(base) , main="Baseline" ,
                     breaks=30 , col="skyblue"  )
xLim=range( c(aSTZ) )
for ( j in 1:NaLvl ) {
  histinfo = plotPost( aSTZ[j,] , xlab=bquote(a[.(j)]) , 
                       main=aLvlNames[j] , 
                       breaks=30 , col=aColor , xlim=xLim )
}
xLim=range( c(bSTZ) )
for ( j in 1:NbLvl ) {
  histinfo = plotPost( bSTZ[j,] , xlab=bquote(b[.(j)]) , 
                       main=bLvlNames[j] , 
                       breaks=30 , col=bColor , xlim=xLim )
}
xLim=range( c(axbSTZ) )
for ( k in 1:NbLvl ) {
  for ( j in 1:NaLvl ) {
    histinfo = plotPost( axbSTZ[j,k,] , xlab=bquote(axb[.(j)*","*.(k)]) , 
                         main=paste(aLvlNames[j],"x",bLvlNames[k]) , 
                         breaks=30 , col=axbColor , xlim=xLim )
  }
}
savePlot( file=paste(fileNameRoot,"ab.eps",sep="") , type="eps" )
savePlot( file=paste(fileNameRoot,"ab.jpg",sep="") , type="jpg" )

# Display contrast analyses
nContrasts = length( aContrastList )
if ( nContrasts > 0 ) {
   nPlotPerRow = 5
   nPlotRow = ceiling(nContrasts/nPlotPerRow)
   nPlotCol = ceiling(nContrasts/nPlotRow)
   windows(3.75*nPlotCol,2.5*nPlotRow)
   layout( matrix(1:(nPlotRow*nPlotCol),nrow=nPlotRow,ncol=nPlotCol,byrow=T) )
   par( mar=c(4,0.5,2.5,0.5) , mgp=c(2,0.7,0) )
   for ( cIdx in 1:nContrasts ) {
       contrast = matrix( aContrastList[[cIdx]],nrow=1) # make it a row matrix
       incIdx = contrast!=0
       histInfo = plotPost( contrast %*% aSTZ , compVal=0 , breaks=30 ,
                xlab=paste( round(contrast[incIdx],2) , aLvlNames[incIdx] ,
                            c(rep("+",sum(incIdx)-1),"") , collapse=" " ) ,
                cex.lab = 1.0 ,
                main=paste( "A Contrast:", names(aContrastList)[cIdx] ) ,
                            col=aColor )
   }
   savePlot( file=paste(fileNameRoot,"aContrasts.eps",sep="") , type="eps" )
   savePlot( file=paste(fileNameRoot,"aContrasts.jpg",sep="") , type="jpg" )
}
#
nContrasts = length( bContrastList )
if ( nContrasts > 0 ) {
   nPlotPerRow = 5
   nPlotRow = ceiling(nContrasts/nPlotPerRow)
   nPlotCol = ceiling(nContrasts/nPlotRow)
   windows(3.75*nPlotCol,2.5*nPlotRow)
   layout( matrix(1:(nPlotRow*nPlotCol),nrow=nPlotRow,ncol=nPlotCol,byrow=T) )
   par( mar=c(4,0.5,2.5,0.5) , mgp=c(2,0.7,0) )
   for ( cIdx in 1:nContrasts ) {
       contrast = matrix( bContrastList[[cIdx]],nrow=1) # make it a row matrix
       incIdx = contrast!=0
       histInfo = plotPost( contrast %*% bSTZ , compVal=0 , breaks=30 ,
                xlab=paste( round(contrast[incIdx],2) , bLvlNames[incIdx] ,
                            c(rep("+",sum(incIdx)-1),"") , collapse=" " ) ,
                cex.lab = 1.0 ,
                main=paste( "B Contrast:", names(bContrastList)[cIdx] ) ,
                            col=bColor )
   }
   savePlot( file=paste(fileNameRoot,"bContrasts.eps",sep="") , type="eps" )
   savePlot( file=paste(fileNameRoot,"bContrasts.jpg",sep="") , type="jpg" )
}
#
nContrasts = length( axbContrastList )
if ( nContrasts > 0 ) {
   nPlotPerRow = 5
   nPlotRow = ceiling(nContrasts/nPlotPerRow)
   nPlotCol = ceiling(nContrasts/nPlotRow)
   windows(3.75*nPlotCol,2.5*nPlotRow)
   layout( matrix(1:(nPlotRow*nPlotCol),nrow=nPlotRow,ncol=nPlotCol,byrow=T) )
   par( mar=c(4,0.5,2.5,0.5) , mgp=c(2,0.7,0) )
   for ( cIdx in 1:nContrasts ) {
       contrast = axbContrastList[[cIdx]]
       contrastArr = array( rep(contrast,chainLength) ,
                            dim=c(NROW(contrast),NCOL(contrast),chainLength) )
       contrastLab = ""
       for ( aIdx in 1:NaLvl ) {
         for ( bIdx in 1:NbLvl ) {
           if ( contrast[aIdx,bIdx] != 0 ) {
             contrastLab = paste( contrastLab , "+" ,
                                  signif(contrast[aIdx,bIdx],2) ,
                                  aLvlNames[aIdx] , bLvlNames[bIdx] )
           }
         }
       }
       histInfo = plotPost( apply( contrastArr * axbSTZ , 3 , sum ) ,
                compVal=0 , breaks=30 , xlab=contrastLab , cex.lab = 0.75 ,
                main=paste( names(axbContrastList)[cIdx] ) , col=axbColor )
   }
   savePlot( file=paste(fileNameRoot,"axbContrasts.eps",sep="") , type="eps" )
   savePlot( file=paste(fileNameRoot,"axbContrasts.jpg",sep="") , type="jpg" )
}

#==============================================================================
# Do NHST ANOVA:
theData = data.frame( y=dataList$y , aLvl=dataList$aLvl , bLvl=dataList$bLvl ,
                      sLvl=dataList$sLvl  )
# ???
#==============================================================================
