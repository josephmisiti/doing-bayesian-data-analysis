setwd('~/learning/bayes_book/ch4_bayes_rule/heads')
# Theta is the vector of candidate values for the parameter theta.
# nThetaVals is the number of candidate theta values.
# To produce the examples in the book, set nThetaVals to either 3 or 63.
nThetaVals <- 20  #times we're going to flip the coin
# Now make the vector of theta values:
# 99 unfair lufte-coins, one fair coin
# Lets make leftward skewed pmf
Theta <- seq(.01, .99, by=0.01)

# pTheta is the vector of prior probabilities on the theta values.
# this is your prior belief. We assign a distribution to our belief of the
# likelihood of different kinds of coins - about how fair or unfair they are
# Lets change this.
# 99 unfair ones that give Heads 90% of time, 1 fair one
pTheta <- dnorm(Theta, mean=0.9, sd=0.1)

pTheta <- pTheta / sum( pTheta )  # Makes sure that beliefs sum to 1.

# Specify the data. To produce the examples in the book, use either
# Get 50 heads, and then 50 tails
for(i in 1:nThetaVals){
heads_flipped <- i
Data <- c(rep(1, heads_flipped), rep(0, nThetaVals - heads_flipped))
nHeads <- sum( Data == 1 )
nTails <- sum( Data == 0 )

# Compute the likelihood of the data for each value of theta:
pDataGivenTheta <- Theta^nHeads * (1-Theta)^nTails

# Compute the posterior:
# This is total probability rule?
pData <- sum( pDataGivenTheta * pTheta )  
pThetaGivenData <- pDataGivenTheta * pTheta / pData   # This is Bayes' rule!
# Do you really need the pData? Or can you just 'normalize all the probs'
# No, you don't!
# you can calc it this way too by normalizing, which is in effect exactly what
# the total probability does
pThetaGivenData1 <- (pDataGivenTheta * pTheta) / sum(pDataGivenTheta * pTheta)
all(pThetaGivenData == pThetaGivenData1)  # Damn right it's TRUE


# which gives you [0.71, 0.30, 0.001]
# which means that if you got 3 heads out of 12 flips, you got damn near no
# chance that you got the coin that comes up 75% heads

# Plot the results.
# windows(7,10) # create window of specified width,height inches.
# x11(7,10) # create window of specified width,height inches.

filename <- paste('heads_', as.character(i), '.jpeg', sep='')
jpeg(filename = filename, width = 480, height = 481, units = "px", 
     pointsize = 12) 
 # quality = 75, bg = "white", res = NA, ..., type = c("cairo", "Xlib", 
 # "quartz"), antialias)
layout( matrix( c( 1,2,3 ) ,nrow=3 ,ncol=1 ,byrow=FALSE ) ) # 3x1 panels
par(mar=c(3,3,1,0))         # number of margin lines: bottom,left,top,right
par(mgp=c(2,1,0))           # which margin lines to use for labels
par(mai=c(0.5,0.5,0.3,0.1)) # margin size in inches: bottom,left,top,right

# Plot the prior:
plot( Theta , pTheta , type="h" , lwd=3 , main="Prior" ,
      xlim=c(0,1) , xlab=bquote(theta) ,
      ylim=c(0,1.1*max(pThetaGivenData)) , ylab=bquote(p(theta)) ,
      cex.axis=1.2 , cex.lab=1.5 , cex.main=1.5 )

# Plot the likelihood:
plot( Theta , pDataGivenTheta , type="h" , lwd=3 , main="Likelihood" ,
      xlim=c(0,1) , xlab=bquote(theta) ,
      ylim=c(0,1.1*max(pDataGivenTheta)) , ylab=bquote(paste("p(D|",theta,")")),
      cex.axis=1.2 , cex.lab=1.5 , cex.main=1.5 )
text( .55 , .85*max(pDataGivenTheta) , cex=2.0 ,
      bquote( "D=" * .(nHeads) * "H," * .(nTails) * "T" ) , adj=c(0,.5) )

# Plot the posterior:
plot( Theta , pThetaGivenData , type="h" , lwd=3 , main="Posterior" ,
      xlim=c(0,1) , xlab=bquote(theta) ,
      ylim=c(0,1.1*max(pThetaGivenData)) , ylab=bquote(paste("p(",theta,"|D)")),
      cex.axis=1.2 , cex.lab=1.5 , cex.main=1.5 )
text( .55 , .85*max(pThetaGivenData) , cex=2.0 ,
      bquote( "p(D)=" * .(signif(pData,3)) ) , adj=c(0,.5) )

# Save the plot as an EPS file.
dev.off()
}
#test comment
