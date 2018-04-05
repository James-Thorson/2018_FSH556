
setwd( "C:/Users/James.Thorson/Desktop/Project_git/2018_FSH556/Week 2 -- mixed-effects/Lab 2/")

png( file="Laplace_example.png", width=10, height=7.5, res=200, units="in")
  par( mfrow=c(2,2), mgp=c(2,0.5,0), tck=-0.02, xaxs="i", yaxs="i", mar=c(2,2,2,1), oma=c(2,2,0,0) )
  for(i in 1:4){
    if(i==1) X = seq(-5,12,length=1000)
    if(i==2) X = seq(-5,25,length=1000)
    if(i==3) X = seq(0,50,length=1000)
    if(i==4) X = seq(0,100,length=1000)
    DF = c(4,8,16,32)[i]
    Y1 = dchisq( X, df=DF )
    Y2 = dnorm( X, mean=DF-2, sd=sqrt(2*(DF-2)) )
    Y2 = Y2/max(Y2)*max(Y1)
    matplot( x=X, y=cbind(Y1,Y2), type="l", lwd=2, col=c("black","red"), ylim=c(0,max(c(Y1,Y2))*1.2), xlab="", ylab="", main=paste("Degrees of freedom = ",DF) )
    abline( v=1, lty="dotted" )
    Integral1 = sum(Y1) * mean(diff(X))
    Integral2 = sum(Y2) * mean(diff(X))
    legend( "topleft", col="black", legend=paste0("Integral= ",formatC(Integral1,format="f",digits=3)), bty="n" )
    legend( "topright", col="red", legend=paste0("Approximation= ",formatC(Integral1,format="f",digits=3)), bty="n" )
    mtext( side=1:2, text=c("Value","Density"), outer=TRUE, cex=1.5 )
  }
dev.off()

