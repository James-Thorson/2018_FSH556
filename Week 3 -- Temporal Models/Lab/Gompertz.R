
setwd( "C:/Users/James.Thorson/Desktop/Project_git/2018_FSH556/Week 3 -- Temporal Models/Lab" )
Use_REML = FALSE
set.seed(2)

#####################
# Explore Gompertz model
#####################

beta= 0.2
alpha = 1
d_equil = exp(alpha/beta)

d_1 = seq(0,d_equil*2,length=1e4)
d_2 = d_1 * exp(alpha) * exp( - beta*log(d_1))

# Dynamics
png( file="Gompertz_dynamics.png", width=8, height=4, res=200, units="in")
  par( mfrow=c(1,2), mar=c(3,3.5,2,0), mgp=c(1.75,0.25,0), tck=-0.02)
  plot( x=d_1, y=d_2, type="l", lwd=3, xlab=expression(Biomass[t]), ylab=expression(Biomass[t+1]), main="Production")
  abline( a=0, b=1, lty="dotted")
  # Log-dynamics
  plot( x=log(d_1[-1]), y=log(d_2[-1]/d_1[-1]), type="l", lwd=3, xlab=expression(log(Biomass[t])), ylab=expression(log(Biomass[t+1]/Biomass[t])), main="log-Biomass ratio" )
  abline( a=1, b=0, lty="dotted")
dev.off()

######################
# Simulate data
######################

nt = 100
log_d0 = 3
sigmaP = 0.5
sigmaM = 0.5
alpha = 1
beta = 0.1

# Simulate predictors
log_d_t = log_b_t = rep(NA, nt)
log_d_t[1] = log_d0
for( t in 2:nt ){
  log_d_t[t] = alpha + (1-beta)*log_d_t[t-1] + rnorm( 1, mean=0, sd=sigmaP )
}
for( t in 1:nt ){
  log_b_t[t] = log_d_t[t] + rnorm( 1, mean=0, sd=sigmaM )
}



######################
# Run in TMB
######################

library(TMB)

# Compile model
Version = "gompertz"
compile( paste0(Version,".cpp") )

# Build inputs
Data = list( "nt"=nt, "log_b_t"=log_b_t )
Parameters = list( "log_d0"=0, "log_sigmaP"=1, "log_sigmaM"=1, "alpha"=0, "rho"=0, "log_d_t"=rep(0,Data$nt) )
Random = c("log_d_t")
if( Use_REML==TRUE ) Random = union( Random, c("log_d0","alpha","rho") )

# Build object
dyn.load( dynlib("gompertz") )
Obj = MakeADFun(data=Data, parameters=Parameters, random=Random, DLL="gompertz")  #

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
  Opt$SD = sdreport( Obj )

for( z in 1:3 ){
  png( file=paste0("gompertz_",z,".png"), width=8, height=5, res=200, units="in" )
    par( mar=c(3,3,1,1), mgp=c(2,0.5,0), tck=-0.02 )
    plot( x=1:Data$nt, y=Data$log_b_t, col="blue", cex=1.2, xlab="Time", ylab="Value", pch=20 )
    if(z>=2){
      points( x=1:Data$nt, y=as.list(Opt$SD,"Estimate")$log_d_t, col="red" )
      for( t in 1:Data$nt) lines( x=rep(t,2), y=as.list(Opt$SD,"Estimate")$log_d_t[t]+c(-1.96,1.96)*as.list(Opt$SD,"Std. Error")$log_d_t[t], col="red" )
    }
    if(z>=3) lines( x=1:Data$nt, y=log_d_t, col="black", lwd=2 )
  dev.off()
}

######################
# Download real data and run again
######################

# devtools::install_github("james-thorson/FishData")

# Download data for Alaska pollock
CPUE = FishData::download_catch_rates( survey="Eastern_Bering_Sea", species_set="Gadus chalcogrammus", error_tol=0.01, localdir=paste0(getwd(),"/") )
B_t = tapply( CPUE[,'Wt'], INDEX=CPUE[,'Year'], FUN=mean )

# Run Gompertz model again
Data = list( "nt"=length(B_t), "log_b_t"=log(B_t) )
Parameters = list( "log_d0"=0, "log_sigmaP"=1, "log_sigmaM"=1, "alpha"=0, "rho"=0, "log_d_t"=rep(0,Data$nt) )
Obj = MakeADFun(data=Data, parameters=Parameters, random=Random, DLL="gompertz")  #
Opt = TMBhelper::Optimize( obj=Obj )

# Get reporting and SEs
Report = Obj$report()

for( z in 1:2 ){
  png( file=paste0("pollock_",z,".png"), width=8, height=5, res=200, units="in" )
    par( mar=c(3,3,1,1), mgp=c(2,0.5,0), tck=-0.02 )
    plot( x=1:Data$nt, y=Data$log_b_t, col="blue", cex=1.2, xlab="Time", ylab="Value", pch=20 )
    if(z>=2){
      points( x=1:Data$nt, y=as.list(Opt$SD,"Estimate")$log_d_t, col="red" )
      for( t in 1:Data$nt) lines( x=rep(t,2), y=as.list(Opt$SD,"Estimate")$log_d_t[t]+c(-1.96,1.96)*as.list(Opt$SD,"Std. Error")$log_d_t[t], col="red" )
    }
  dev.off()
}

write.csv( cbind(Opt$diagnostics, summary(Opt$SD,"fixed")[,'Std. Error']), file=paste0(getwd(),"/pollock.csv") )

