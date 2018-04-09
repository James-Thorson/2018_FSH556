
setwd( "C:/Users/James.Thorson/Desktop/Project_git/2018_FSH556/Week 3 -- Temporal Models/Lecture" )
Use_REML = FALSE
set.seed(2)

######################
# Simulate data
######################

nt = 100
x0 = 3
sigmaP = 0.2
sigmaM = 0.2
alpha = 0.04

# Simulate predictors
x_t = y_t = rep(NA, nt)
x_t[1] = x0
for( t in 2:nt ){
  x_t[t] = x_t[t-1] + rnorm( 1, mean=alpha, sd=sigmaP )
}
for( t in 1:nt ){
  y_t[t] = x_t[t] + rnorm( 1, mean=0, sd=sigmaM )
}

######################
# Linear model in R
######################

Lm = lm( y_t ~ I(1:nt) )
ypred_t = predict(Lm, se=TRUE)

for( z in 1:3 ){
  png( file=paste0("lm_",z,".png"), width=8, height=5, res=200, units="in" )
    par( mar=c(3,3,1,1), mgp=c(2,0.5,0), tck=-0.02 )
    plot( x=1:nt, y=y_t, col="blue", cex=1.2, xlab="Time", ylab="Value", pch=20 )
    if(z>=2){
      points( x=1:nt, y=ypred_t$fit, col="red" )
      for( t in 1:nt) lines( x=rep(t,2), y=ypred_t$fit[t]+c(-1.96,1.96)*ypred_t$se.fit[t], col="red" )
    }
    if(z>=3) lines( x=1:nt, y=x_t, col="black", lwd=2 )
  dev.off()
}


######################
# Loess smoother in R
######################

Loess = loess( y_t ~ I(1:nt) )
ypred_t = predict(Loess, se=TRUE)

for( z in 1:3 ){
  png( file=paste0("loess_",z,".png"), width=8, height=5, res=200, units="in" )
    par( mar=c(3,3,1,1), mgp=c(2,0.5,0), tck=-0.02 )
    plot( x=1:nt, y=y_t, col="blue", cex=1.2, xlab="Time", ylab="Value", pch=20 )
    if(z>=2){
      points( x=1:nt, y=ypred_t$fit, col="red" )
      for( t in 1:nt) lines( x=rep(t,2), y=ypred_t$fit[t]+c(-1.96,1.96)*ypred_t$se.fit[t], col="red" )
    }
    if(z>=3) lines( x=1:nt, y=x_t, col="black", lwd=2 )
  dev.off()
}


######################
# Smoother in R
######################

library(mgcv)

Gam = gam( y_t ~ s(I(1:nt)) )
ypred_t = predict(Gam, se=TRUE)

for( z in 1:3 ){
  png( file=paste0("gam_",z,".png"), width=8, height=5, res=200, units="in" )
    par( mar=c(3,3,1,1), mgp=c(2,0.5,0), tck=-0.02 )
    plot( x=1:nt, y=y_t, col="blue", cex=1.2, xlab="Time", ylab="Value", pch=20 )
    if(z>=2){
      points( x=1:nt, y=ypred_t$fit, col="red" )
      for( t in 1:nt) lines( x=rep(t,2), y=ypred_t$fit[t]+c(-1.96,1.96)*ypred_t$se.fit[t], col="red" )
    }
    if(z>=3) lines( x=1:nt, y=x_t, col="black", lwd=2 )
  dev.off()
}

######################
# Run in TMB
######################

library(TMB)

# Compile model
Version = "kalman_filter"
compile( paste0(Version,".cpp") )

# Build inputs
Data = list( "nt"=nt, "y_t"=y_t )
Parameters = list( "x0"=0, "log_sigmaP"=1, "log_sigmaM"=1, "alpha"=0, "x_t"=rep(0,nt) )
Random = c("x_t")
if( Use_REML==TRUE ) Random = union( Random, c("x0","alpha") )

# Build object
dyn.load( dynlib("kalman_filter") )
Obj = MakeADFun(data=Data, parameters=Parameters, random=Random)  #

# Prove that function and gradient calls work
Obj$fn( Obj$par )
Obj$gr( Obj$par )

# Optimize
start_time = Sys.time()
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, control=list("trace"=1) )
  Opt[["final_gradient"]] = Obj$gr( Opt$par )
  Opt[["total_time"]] = Sys.time() - start_time

# Get reporting and SEs
Report = Obj$report()
  SD = sdreport( Obj )

for( z in 1:3 ){
  png( file=paste0("kalman_",z,".png"), width=8, height=5, res=200, units="in" )
    par( mar=c(3,3,1,1), mgp=c(2,0.5,0), tck=-0.02 )
    plot( x=1:nt, y=y_t, col="blue", cex=1.2, xlab="Time", ylab="Value", pch=20 )
    if(z>=2){
      points( x=1:nt, y=as.list(SD,"Estimate")$x_t, col="red" )
      for( t in 1:nt) lines( x=rep(t,2), y=as.list(SD,"Estimate")$x_t[t]+c(-1.96,1.96)*as.list(SD,"Std. Error")$x_t[t], col="red" )
    }
    if(z>=3) lines( x=1:nt, y=x_t, col="black", lwd=2 )
  dev.off()
}

